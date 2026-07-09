#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00_calibrate_init_eir.R
# ------------------------------------------------------------------------------
# Careful-workflow stage 0: solve for the init_EIR that yields the target
# *seasonal-mean* PfPR2-10 (microscopy). Two-step procedure:
#
#   1. Analytical anchor: malariaEquilibrium-based inversion via the
#      production helper `calibrate_eir_from_pfpr()`. This is fast (~5 s)
#      but ignores seasonality, so it overshoots the target under the
#      seasonal simulator.
#   2. Simulation-based refinement: 4-point log-spaced coarse grid plus
#      log-linear interpolation, then up to 3 refinement iterations. Each
#      grid point runs a short warmup (default 2 years) through
#      `msimGD_run_truth()` (the production helper) and measures the realised
#      annual-mean PfPR over the final year. Total runtime ~3-5 minutes.
#
# This matches the production audit-cell pattern: analytical inversion as a
# starting point, then a real simulation-based calibration on top.
#
# Output (under output/careful/):
#   calibrated_init_eir.rds  list(init_EIR, target, realised, search_log)
#   calibration_log.csv      all (init_EIR, realised_pfpr) pairs visited
#   calibration_search_diagnostics.csv / .txt
#                            bracketing and monotonicity diagnostics for the
#                            candidate-EIR search
#   calibration_node_summary.csv
#                            final validation-run node-level PfPR, EIR,
#                            adult mosquitoes, and clinical incidence
#   calibration_global_summary.csv
#                            population-weighted/global final validation
#                            summaries and incidence reference comparison
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/Busia_landscape/careful/00_calibrate_init_eir.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

load_local_malariasimulationGD <- function(example_dir) {
  pkg_root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(pkg_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (requireNamespace("malariasimulationGD", quietly = TRUE)) {
    suppressPackageStartupMessages(library(malariasimulationGD))
    return(invisible(TRUE))
  }
  stop("malariasimulationGD is not installed and pkgload is unavailable.",
       call. = FALSE)
}

load_local_malariasimulationGD(example_dir)
suppressPackageStartupMessages(library(malariaEquilibrium))

source(file.path(example_dir, "lib", "movement_mu.R"))
source(file.path(example_dir, "lib", "calibrate_eir.R"))
source(file.path(example_dir, "lib", "msimGD_truth_generation.R"))
source(file.path(example_dir, "lib", "synthetic_covariate.R"))
source(file.path(example_dir, "lib", "simulation_calibration.R"))

source(file.path(example_dir, "config", "busia_landscape.R"))
source(file.path(example_dir, "config", "movement.R"))
source(file.path(example_dir, "config", "seasonality.R"))
source(file.path(example_dir, "config", "covariate.R"))
source(file.path(example_dir, "config", "homing_drive.R"))

out_dir <- Sys.getenv(
  "BUSIA_CAREFUL_OUT_DIR",
  unset = file.path(example_dir, "output", "careful")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Calibration knobs ------------------------------------------------------

TARGET_PFPR    <- busia_default_target_pfpr()
TARGET_AGE_MIN_YEARS <- 2
TARGET_AGE_MAX_YEARS <- 10
TARGET_AGE_MIN_DAYS  <- as.integer(round(365 * TARGET_AGE_MIN_YEARS))
TARGET_AGE_MAX_DAYS  <- as.integer(round(365 * TARGET_AGE_MAX_YEARS))

# Production mosquito-biology lifecycle
theta <- list(
  qE = 1 / 3, nE = 2L,
  qL = 1 / 7, nL = 3L,
  qP = 1,     nP = 2L,
  muE = 0.05, muL = 0.15, muP = 0.05,
  muF = 0.132, muM = 0.132,
  beta = 16,
  nEIP = 5L,
  nu = 1 / (4 / 24)
)
saveRDS(theta, file.path(out_dir, "theta.rds"))

# Short-warmup length used for each calibration probe (years). Production
# uses simulation runs that are similar in scale.
SEARCH_WARMUP_DAYS <- as.integer(Sys.getenv("BUSIA_SEARCH_WARMUP_DAYS", "3650"))
if (is.na(SEARCH_WARMUP_DAYS) || SEARCH_WARMUP_DAYS < 1L) {
  stop("BUSIA_SEARCH_WARMUP_DAYS must be a positive integer.", call. = FALSE)
}
MEASUREMENT_WINDOW_DAYS <- min(365L, SEARCH_WARMUP_DAYS)
MAX_REFINE_ITER <- as.integer(Sys.getenv("BUSIA_CALIBRATION_MAX_REFINE", "3"))
if (is.na(MAX_REFINE_ITER) || MAX_REFINE_ITER < 1L) {
  stop("BUSIA_CALIBRATION_MAX_REFINE must be a positive integer.", call. = FALSE)
}
MONOTONICITY_ABS_TOL <- as.numeric(
  Sys.getenv("BUSIA_MONOTONICITY_ABS_TOL", "0.05")
)
if (!is.finite(MONOTONICITY_ABS_TOL) || MONOTONICITY_ABS_TOL < 0) {
  stop("BUSIA_MONOTONICITY_ABS_TOL must be a finite non-negative number.",
       call. = FALSE)
}
BUSIA_REQUIRE_CALIBRATION_DIAGNOSTICS <- tolower(Sys.getenv(
  "BUSIA_REQUIRE_CALIBRATION_DIAGNOSTICS",
  unset = "false"
)) %in% c("true", "t", "1", "yes", "y")
REFERENCE_INCIDENCE_PER_PY <- 3.2
RNG_SEED <- 20260514L

required_col <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) < 1L) {
    stop(sprintf(
      "Could not find %s column. Available columns: %s",
      label, paste(names(df), collapse = ", ")
    ), call. = FALSE)
  }
  hit[[1L]]
}

extract_data_list <- function(res) {
  if (is.list(res) && all(vapply(res, is.data.frame, logical(1)))) {
    res
  } else if (!is.null(res$data)) {
    res$data
  } else if (!is.null(res$sim$data)) {
    res$sim$data
  } else {
    stop("Could not locate per-node data.frames in msimGD_run_truth output.",
         call. = FALSE)
  }
}

search_diagnostics <- function(search_log, target_pfpr, monotonicity_abs_tol) {
  if (!is.data.frame(search_log) || nrow(search_log) < 1L) {
    stop("cal_sim$search_log must be a non-empty data.frame.", call. = FALSE)
  }
  eir_col <- required_col(search_log, c("init_EIR"), "candidate init_EIR")
  pfpr_col <- required_col(search_log, c("realised_pfpr"), "realised PfPR")
  raw <- data.frame(
    init_EIR = as.numeric(search_log[[eir_col]]),
    realised_pfpr = as.numeric(search_log[[pfpr_col]])
  )
  if (any(!is.finite(raw$init_EIR)) || any(raw$init_EIR <= 0) ||
      any(!is.finite(raw$realised_pfpr))) {
    stop("Calibration search log contains non-finite or non-positive candidate values.",
         call. = FALSE)
  }
  agg <- aggregate(raw$realised_pfpr,
                   by = list(init_EIR = raw$init_EIR), FUN = mean)
  names(agg)[2] <- "realised_pfpr"
  agg <- agg[order(agg$init_EIR), , drop = FALSE]

  below <- agg[agg$realised_pfpr <= target_pfpr, , drop = FALSE]
  above <- agg[agg$realised_pfpr >= target_pfpr, , drop = FALSE]
  below_row <- if (nrow(below) > 0L) below[which.max(below$realised_pfpr), ] else NULL
  above_row <- if (nrow(above) > 0L) above[which.min(above$realised_pfpr), ] else NULL
  bracketed <- nrow(below) > 0L && nrow(above) > 0L

  deltas <- diff(agg$realised_pfpr)
  violation_drops <- pmax(-deltas - monotonicity_abs_tol, 0)
  violation_idx <- which(deltas < -monotonicity_abs_tol)
  spearman <- if (nrow(agg) >= 2L) {
    suppressWarnings(stats::cor(log(agg$init_EIR), agg$realised_pfpr,
                                method = "spearman"))
  } else {
    NA_real_
  }

  list(
    aggregated_search = agg,
    summary = data.frame(
      target_pfpr = target_pfpr,
      n_search_rows = nrow(search_log),
      n_unique_init_EIR = nrow(agg),
      min_init_EIR = min(agg$init_EIR),
      max_init_EIR = max(agg$init_EIR),
      min_realised_pfpr = min(agg$realised_pfpr),
      max_realised_pfpr = max(agg$realised_pfpr),
      target_bracketed = bracketed,
      below_init_EIR = if (is.null(below_row)) NA_real_ else below_row$init_EIR,
      below_realised_pfpr = if (is.null(below_row)) NA_real_ else below_row$realised_pfpr,
      above_init_EIR = if (is.null(above_row)) NA_real_ else above_row$init_EIR,
      above_realised_pfpr = if (is.null(above_row)) NA_real_ else above_row$realised_pfpr,
      spearman_log_eir_pfpr = spearman,
      monotonicity_abs_tol = monotonicity_abs_tol,
      adjacent_monotonicity_violations = length(violation_idx),
      max_adjacent_decline = if (length(deltas) > 0L) max(pmax(-deltas, 0)) else 0,
      total_material_decline_excess = sum(violation_drops),
      approximately_non_decreasing = length(violation_idx) == 0L,
      stringsAsFactors = FALSE
    )
  )
}

write_search_diagnostics <- function(diag, out_dir) {
  summary_df <- diag$summary
  utils::write.csv(summary_df,
                   file.path(out_dir, "calibration_search_diagnostics.csv"),
                   row.names = FALSE)
  lines <- c(
    "Busia calibration search diagnostics",
    sprintf("Target PfPR2-10 microscopy: %.3f", summary_df$target_pfpr),
    sprintf("Explored init_EIR range: %.6f to %.6f",
            summary_df$min_init_EIR, summary_df$max_init_EIR),
    sprintf("Realised PfPR range: %.6f to %.6f",
            summary_df$min_realised_pfpr, summary_df$max_realised_pfpr),
    sprintf("Target bracketed by tested simulations: %s",
            summary_df$target_bracketed),
    sprintf("Candidate below target: init_EIR=%s, PfPR=%s",
            ifelse(is.na(summary_df$below_init_EIR), "NA",
                   sprintf("%.6f", summary_df$below_init_EIR)),
            ifelse(is.na(summary_df$below_realised_pfpr), "NA",
                   sprintf("%.6f", summary_df$below_realised_pfpr))),
    sprintf("Candidate above target: init_EIR=%s, PfPR=%s",
            ifelse(is.na(summary_df$above_init_EIR), "NA",
                   sprintf("%.6f", summary_df$above_init_EIR)),
            ifelse(is.na(summary_df$above_realised_pfpr), "NA",
                   sprintf("%.6f", summary_df$above_realised_pfpr))),
    sprintf("Spearman cor(log(init_EIR), realised PfPR): %s",
            ifelse(is.na(summary_df$spearman_log_eir_pfpr), "NA",
                   sprintf("%.6f", summary_df$spearman_log_eir_pfpr))),
    sprintf("Adjacent monotonicity violations > %.3f: %d",
            summary_df$monotonicity_abs_tol,
            summary_df$adjacent_monotonicity_violations),
    sprintf("Maximum adjacent decline: %.6f",
            summary_df$max_adjacent_decline),
    sprintf("Approximately non-decreasing: %s",
            summary_df$approximately_non_decreasing)
  )
  writeLines(lines, file.path(out_dir, "calibration_search_diagnostics.txt"))
}

validate_search_diagnostics <- function(diag, strict) {
  summary_df <- diag$summary
  problems <- character(0)
  if (!isTRUE(summary_df$target_bracketed)) {
    problems <- c(problems, "target PfPR is not bracketed by tested simulations")
  }
  if (!isTRUE(summary_df$approximately_non_decreasing)) {
    problems <- c(
      problems,
      sprintf(
        "search has %d adjacent monotonicity violation(s) above tolerance %.3f",
        summary_df$adjacent_monotonicity_violations,
        summary_df$monotonicity_abs_tol
      )
    )
  }
  if (length(problems) > 0L) {
    msg <- paste(
      "Calibration search diagnostics flagged:",
      paste(problems, collapse = "; ")
    )
    if (isTRUE(strict)) {
      stop(msg, call. = FALSE)
    }
    warning(msg, call. = FALSE)
  }
}

summarise_validation_run <- function(data_list, nodes, measurement_window_days,
                                     age_band_days, clinical_age_days,
                                     reference_incidence_per_py,
                                     realised_pfpr) {
  if (length(data_list) != nrow(nodes)) {
    stop(sprintf("Validation run returned %d node outputs, expected %d.",
                 length(data_list), nrow(nodes)), call. = FALSE)
  }
  prev_n_col <- sprintf("n_detect_lm_%d_%d",
                        age_band_days[[1]], age_band_days[[2]])
  prev_d_col <- sprintf("n_age_%d_%d",
                        age_band_days[[1]], age_band_days[[2]])
  clinical_col <- sprintf("n_inc_clinical_%d_%d",
                          clinical_age_days[[1]], clinical_age_days[[2]])
  clinical_age_col <- sprintf("n_age_%d_%d",
                              clinical_age_days[[1]], clinical_age_days[[2]])
  required <- c("timestep", prev_n_col, prev_d_col, clinical_col,
                clinical_age_col)
  optional <- c("EIR_gamb", "total_M_gamb")

  node_rows <- vector("list", length(data_list))
  unavailable <- character(0)
  for (i in seq_along(data_list)) {
    df <- data_list[[i]]
    missing <- setdiff(required, names(df))
    if (length(missing) > 0L) {
      stop(sprintf("Validation output for node %d is missing required columns: %s",
                   i, paste(missing, collapse = ", ")), call. = FALSE)
    }
    max_t <- max(df$timestep)
    start_t <- max_t - measurement_window_days + 1L
    win <- df[df$timestep >= start_t & df$timestep <= max_t, , drop = FALSE]
    if (nrow(win) != measurement_window_days) {
      stop(sprintf(
        "Validation window for node %d expected %d days but found %d.",
        i, measurement_window_days, nrow(win)
      ), call. = FALSE)
    }
    for (col in optional) {
      if (!(col %in% names(df))) {
        unavailable <- c(unavailable, sprintf("node %d: %s unavailable", i, col))
      }
    }
    prev_den <- sum(win[[prev_d_col]])
    clinical_person_days <- sum(win[[clinical_age_col]])
    clinical_cases <- sum(win[[clinical_col]])
    node_name <- if ("node_name" %in% names(nodes)) {
      as.character(nodes$node_name[[i]])
    } else if ("name" %in% names(nodes)) {
      as.character(nodes$name[[i]])
    } else {
      NA_character_
    }
    node_rows[[i]] <- data.frame(
      node_id = as.integer(nodes$node[[i]]),
      node_name = node_name,
      human_population = as.integer(nodes$NH_per_node[[i]]),
      population_weight = as.numeric(nodes$NH_per_node[[i]]) / sum(nodes$NH_per_node),
      pfpr2_10_microscopy = sum(win[[prev_n_col]]) / pmax(prev_den, 1),
      mean_eir = if ("EIR_gamb" %in% names(win)) mean(win$EIR_gamb) else NA_real_,
      mean_adult_female_mosquitoes = if ("total_M_gamb" %in% names(win)) {
        mean(win$total_M_gamb)
      } else {
        NA_real_
      },
      clinical_incident_cases = clinical_cases,
      person_years_at_risk = clinical_person_days / 365,
      clinical_incidence_per_person_year =
        clinical_cases / pmax(clinical_person_days / 365, .Machine$double.eps),
      diagnostics = if (is.na(node_name)) "node_name unavailable in nodes table" else "",
      stringsAsFactors = FALSE
    )
  }
  node_summary <- do.call(rbind, node_rows)
  global_pfpr <- stats::weighted.mean(
    node_summary$pfpr2_10_microscopy,
    node_summary$population_weight
  )
  global_cases <- sum(node_summary$clinical_incident_cases)
  global_py <- sum(node_summary$person_years_at_risk)
  global_summary <- data.frame(
    scope = "global",
    n_nodes = nrow(node_summary),
    measurement_window_days = measurement_window_days,
    pfpr2_10_microscopy = global_pfpr,
    calibration_realised_pfpr = realised_pfpr,
    abs_difference_from_calibration_realised_pfpr = abs(global_pfpr - realised_pfpr),
    clinical_incident_cases = global_cases,
    person_years_at_risk = global_py,
    clinical_incidence_per_person_year = global_cases / pmax(global_py, .Machine$double.eps),
    reference_incidence_per_person_year = reference_incidence_per_py,
    difference_from_reference =
      global_cases / pmax(global_py, .Machine$double.eps) - reference_incidence_per_py,
    ratio_to_reference =
      global_cases / pmax(global_py, .Machine$double.eps) / reference_incidence_per_py,
    incidence_age_definition = sprintf(
      "all ages rendered as %d-%d days; person-years from realised daily n_age_%d_%d",
      clinical_age_days[[1]], clinical_age_days[[2]],
      clinical_age_days[[1]], clinical_age_days[[2]]
    ),
    unavailable_output_fields = paste(unique(unavailable), collapse = "; "),
    stringsAsFactors = FALSE
  )
  list(node_summary = node_summary, global_summary = global_summary)
}

cat(sprintf("[%s] careful/00: simulation-based calibration to PfPR%d-%d = %.3f\n",
            format(Sys.time(), "%F %T"),
            TARGET_AGE_MIN_YEARS, TARGET_AGE_MAX_YEARS, TARGET_PFPR))

# ============================================================================
# Step 1 — analytical anchor (non-seasonal, fast)
# ============================================================================

land <- build_busia_landscape()
validate_busia_landscape(land)
NH_per_node <- as.integer(land$nodes$NH_per_node)
validate_busia_demography_for_NH(NH_per_node)
NH_mean     <- as.integer(round(mean(NH_per_node)))
seas        <- seven_node_seasonality()

cat(sprintf("Busia population input: total = %d | mean = %.1f | sd = %.1f\n",
            sum(NH_per_node), mean(NH_per_node), stats::sd(NH_per_node)))

baseline_p <- msimGD_build_calibration_parameters(
  NH = NH_mean,
  theta = theta,
  overrides = list(
    individual_mosquitoes = FALSE,
    native_mosquito_backend = TRUE,
    model_seasonality = TRUE,
    g0 = seas$g0, g = seas$g, h = seas$h,
    rainfall_floor = seas$rainfall_floor,
    progress_bar = FALSE
  )
)
baseline_p <- apply_busia_demography(baseline_p, total_population = NH_mean)
t0 <- Sys.time()
cal_analytical <- calibrate_eir_from_pfpr(
  target_pfpr      = TARGET_PFPR,
  parameters       = baseline_p,
  age_min          = TARGET_AGE_MIN_YEARS,
  age_max          = TARGET_AGE_MAX_YEARS,
  prevalence_col   = "pos_M",
  eir_range        = c(0.01, 500),
  tol              = 1e-6
)
init_EIR_analytical <- cal_analytical$init_EIR %||%
  cal_analytical$eir_root %||% cal_analytical$init_EIR_calibrated %||%
  stop("Could not extract analytical init_EIR.")
analytical_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Step 1 (analytical anchor): init_EIR = %.4f  [%.1f s]\n",
            init_EIR_analytical, analytical_seconds))

# ============================================================================
# Step 2 — simulation-based refinement under seasonality
# ============================================================================

mv_settings <- seven_node_movement_settings()
cov         <- build_seven_node_covariate(land$D)
cube        <- build_seven_node_drive_cube()

allowed_mat <- busia_allowed_matrix(land$n_nodes)
setup <- list(D = land$D, allowed = allowed_mat)

contact_surface <- list(
  type = "contact_surface",
  contact_multiplier = stats::setNames(
    as.numeric(cov$contact_multiplier),
    as.character(seq_len(land$n_nodes))
  )
)

parameter_modifier <- function(parameters, node_index, warmup_days) {
  parameters <- apply_busia_demography(
    parameters,
    total_population = parameters$human_population
  )
  parameters$model_seasonality <- TRUE
  parameters$g0 <- seas$g0; parameters$g <- seas$g; parameters$h <- seas$h
  parameters$rainfall_floor <- seas$rainfall_floor
  malariasimulationGD::apply_node_contact_surface(
    parameters = parameters, contact_surface = contact_surface,
    node_index = as.integer(node_index)
  )
}

cat(sprintf("Step 2 (simulation-based refinement, ~%.0f-min budget):\n",
            (5 * SEARCH_WARMUP_DAYS / 730) * 0.5))
t0 <- Sys.time()
cal_sim <- simulation_calibrate_init_eir(
  target_pfpr        = TARGET_PFPR,
  setup              = setup,
  cube               = cube,
  NH                 = NH_per_node,
  mu                 = mv_settings$mu,
  p_move             = mv_settings$p_move,
  theta              = theta,
  parameter_modifier = parameter_modifier,
  warmup_days        = SEARCH_WARMUP_DAYS,
  measurement_window_days = MEASUREMENT_WINDOW_DAYS,
  age_band_days      = c(TARGET_AGE_MIN_DAYS, TARGET_AGE_MAX_DAYS),
  starting_init_EIR  = init_EIR_analytical,
  tolerance_rel      = 0.05,
  max_refine_iter    = MAX_REFINE_ITER,
  seed               = RNG_SEED,
  verbose            = TRUE
)
sim_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Step 2 done in %.1f s (%.1f min).\n",
            sim_seconds, sim_seconds / 60))
cat(sprintf("Calibrated init_EIR = %.4f  | realised PfPR = %.4f  (target %.4f)\n",
            cal_sim$init_EIR, cal_sim$realised_pfpr, TARGET_PFPR))

search_diag <- search_diagnostics(
  search_log = cal_sim$search_log,
  target_pfpr = TARGET_PFPR,
  monotonicity_abs_tol = MONOTONICITY_ABS_TOL
)
write_search_diagnostics(search_diag, out_dir)
validate_search_diagnostics(
  search_diag,
  strict = BUSIA_REQUIRE_CALIBRATION_DIAGNOSTICS
)

# ============================================================================
# Step 3 — final no-release validation and incidence reporting
# ============================================================================

busia_demography <- busia_demography_spec(total_population = sum(NH_per_node))
CLINICAL_INCIDENCE_AGE_MIN_DAYS <- 0L
CLINICAL_INCIDENCE_AGE_MAX_DAYS <- max(busia_demography$age_max_days)
INCIDENCE_AGE_DEFINITION <- sprintf(
  "all ages, rendered as %d-%d days; person-years from realised daily n_age_%d_%d over trailing %d days",
  CLINICAL_INCIDENCE_AGE_MIN_DAYS,
  CLINICAL_INCIDENCE_AGE_MAX_DAYS,
  CLINICAL_INCIDENCE_AGE_MIN_DAYS,
  CLINICAL_INCIDENCE_AGE_MAX_DAYS,
  MEASUREMENT_WINDOW_DAYS
)

cat(sprintf("Step 3 (final validation run at calibrated init_EIR = %.4f):\n",
            cal_sim$init_EIR))
t0 <- Sys.time()
validation_res <- msimGD_run_truth(
  setup        = setup,
  cube         = cube,
  NF           = NULL,
  NH           = NH_per_node,
  tmax         = SEARCH_WARMUP_DAYS,
  mu           = mv_settings$mu,
  p_move       = mv_settings$p_move,
  release      = NULL,
  theta        = theta,
  init_EIR     = cal_sim$init_EIR,
  prevalence_rendering_min_age = TARGET_AGE_MIN_DAYS,
  prevalence_rendering_max_age = TARGET_AGE_MAX_DAYS,
  clinical_incidence_min_age   = CLINICAL_INCIDENCE_AGE_MIN_DAYS,
  clinical_incidence_max_age   = CLINICAL_INCIDENCE_AGE_MAX_DAYS,
  warmup_days  = 0L,
  parameter_modifier = parameter_modifier,
  seed = RNG_SEED
)
validation_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
validation_data <- extract_data_list(validation_res)
validation_summary <- summarise_validation_run(
  data_list = validation_data,
  nodes = land$nodes,
  measurement_window_days = MEASUREMENT_WINDOW_DAYS,
  age_band_days = c(TARGET_AGE_MIN_DAYS, TARGET_AGE_MAX_DAYS),
  clinical_age_days = c(CLINICAL_INCIDENCE_AGE_MIN_DAYS,
                        CLINICAL_INCIDENCE_AGE_MAX_DAYS),
  reference_incidence_per_py = REFERENCE_INCIDENCE_PER_PY,
  realised_pfpr = cal_sim$realised_pfpr
)
node_summary <- validation_summary$node_summary
global_summary <- validation_summary$global_summary
global_summary$validation_seconds <- validation_seconds

utils::write.csv(node_summary,
                 file.path(out_dir, "calibration_node_summary.csv"),
                 row.names = FALSE)
utils::write.csv(global_summary,
                 file.path(out_dir, "calibration_global_summary.csv"),
                 row.names = FALSE)

if (global_summary$abs_difference_from_calibration_realised_pfpr > 0.05) {
  warning(sprintf(
    paste(
      "Final validation population-weighted node PfPR %.4f differs from",
      "cal_sim$realised_pfpr %.4f by %.4f."
    ),
    global_summary$pfpr2_10_microscopy,
    cal_sim$realised_pfpr,
    global_summary$abs_difference_from_calibration_realised_pfpr
  ), call. = FALSE)
}

cat(sprintf("Step 3 done in %.1f s (%.1f min).\n",
            validation_seconds, validation_seconds / 60))
cat(sprintf("Global final-year PfPR2-10: %.3f (target %.3f)\n",
            global_summary$pfpr2_10_microscopy, TARGET_PFPR))
cat(sprintf("Global clinical incidence: %.3f episodes/person-year\n",
            global_summary$clinical_incidence_per_person_year))
cat(sprintf("Reference incidence: %.3f episodes/person-year\n",
            REFERENCE_INCIDENCE_PER_PY))

# --- Persist ----------------------------------------------------------------

saveRDS(
  list(
    init_EIR             = cal_sim$init_EIR,
    realised_pfpr        = cal_sim$realised_pfpr,
    target_prevalence    = TARGET_PFPR,
    age_min_years        = TARGET_AGE_MIN_YEARS,
    age_max_years        = TARGET_AGE_MAX_YEARS,
    age_min_days         = TARGET_AGE_MIN_DAYS,
    age_max_days         = TARGET_AGE_MAX_DAYS,
    NH_used_for_analytical = NH_mean,
    NH_per_node           = NH_per_node,
    population_total      = sum(NH_per_node),
    population_mean       = mean(NH_per_node),
    population_sd         = stats::sd(NH_per_node),
    busia_age_structure   = read_busia_age_structure(),
    busia_demography      = busia_demography,
    theta                = theta,
    init_EIR_analytical  = init_EIR_analytical,
    analytical_pfpr      = TARGET_PFPR,    # by construction the analytical hit exactly
    search_warmup_days   = SEARCH_WARMUP_DAYS,
    measurement_window_days = MEASUREMENT_WINDOW_DAYS,
    analytical_seconds   = analytical_seconds,
    simulation_seconds   = sim_seconds,
    validation_seconds   = validation_seconds,
    search_diagnostics   = search_diag$summary,
    node_summary         = node_summary,
    global_summary       = global_summary,
    reference_incidence_per_py = REFERENCE_INCIDENCE_PER_PY,
    incidence_age_definition = INCIDENCE_AGE_DEFINITION
  ),
  file.path(out_dir, "calibrated_init_eir.rds")
)
utils::write.csv(cal_sim$search_log,
                 file.path(out_dir, "calibration_log.csv"), row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(out_dir, "calibrated_init_eir.rds")))
cat(sprintf("Saved: %s\n", file.path(out_dir, "calibration_log.csv")))
cat(sprintf("Saved: %s\n", file.path(out_dir, "calibration_search_diagnostics.csv")))
cat(sprintf("Saved: %s\n", file.path(out_dir, "calibration_search_diagnostics.txt")))
cat(sprintf("Saved: %s\n", file.path(out_dir, "calibration_node_summary.csv")))
cat(sprintf("Saved: %s\n", file.path(out_dir, "calibration_global_summary.csv")))
