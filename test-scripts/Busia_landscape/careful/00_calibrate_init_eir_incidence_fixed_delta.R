#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00_calibrate_init_eir_incidence_fixed_delta.R
# ------------------------------------------------------------------------------
# Faster incidence-calibration stage 0 for the Busia careful workflow.
#
# This version keeps the synthetic contact-covariate coefficient fixed at the
# previous/default value from seven_node_covariate_settings()$delta, unless
# BUSIA_FIXED_DELTA is set. It calibrates raw counterfactual no-intervention
# init_EIR to the target mean modeled infection incidence after baseline ITNs
# and ACT clinical treatment are applied at runtime. The target age band is
# children aged 6 months to younger than 10 years:
#
#   n_inc_182_3650 / (n_age_182_3650 / 365)
#
# Targets and reporting:
#   mean incidence = 3.138 infections/person-year  (calibration target)
#   site-level SD  = 1.381 infections/person-year  (reported only)
#   hazard CV      = 0.48                          (reported only)
#   observed range = 1.1-6.9 infections/person-year (reported only)
#
# The simulator exposes modeled infection incidence, not a microscopy-confirmed
# infection-incidence counter; outputs are labelled accordingly.
#
# Outputs, under BUSIA_CAREFUL_OUT_DIR or output/careful:
#   calibrated_init_eir_incidence_fixed_delta.rds
#   incidence_calibration_fixed_delta_log.csv
#   incidence_calibration_fixed_delta_node_summary.csv
#   incidence_calibration_fixed_delta_global_summary.csv
#   incidence_calibration_fixed_delta_no_intervention_node_summary.csv
#   incidence_calibration_fixed_delta_no_intervention_global_summary.csv
#   theta.rds
#
# Set BUSIA_WRITE_STAGE00_COMPAT=true to also write calibrated_init_eir.rds.
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/Busia_landscape/careful/00_calibrate_init_eir_incidence_fixed_delta.R",
    mustWork = TRUE
  )
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

source(file.path(example_dir, "lib", "movement_mu.R"))
source(file.path(example_dir, "lib", "msimGD_truth_generation.R"))
source(file.path(example_dir, "lib", "synthetic_covariate.R"))

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

parse_logical_env <- function(name, default = FALSE) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) {
    return(default)
  }
  tolower(val) %in% c("true", "t", "1", "yes", "y")
}

parse_numeric_env <- function(name, default, positive = FALSE,
                              non_negative = FALSE, allow_inf = FALSE) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) {
    return(default)
  }
  out <- suppressWarnings(as.numeric(val))
  finite_ok <- is.finite(out) || (isTRUE(allow_inf) && is.infinite(out))
  if (length(out) != 1L || !finite_ok ||
      (positive && out <= 0) || (non_negative && out < 0)) {
    stop(sprintf("%s must be a valid numeric scalar.", name), call. = FALSE)
  }
  out
}

parse_integer_env <- function(name, default, minimum = NULL) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) {
    return(as.integer(default))
  }
  out <- suppressWarnings(as.integer(val))
  if (length(out) != 1L || is.na(out) ||
      (!is.null(minimum) && out < as.integer(minimum))) {
    stop(sprintf("%s must be a valid integer.", name), call. = FALSE)
  }
  as.integer(out)
}

parse_numeric_vector_env <- function(name, default) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) {
    return(as.numeric(default))
  }
  stripped <- gsub("[c()]", " ", val)
  parts <- unlist(strsplit(stripped, "[,[:space:]]+"), use.names = FALSE)
  parts <- parts[nzchar(parts)]
  out <- suppressWarnings(as.numeric(parts))
  if (length(out) < 1L || any(!is.finite(out))) {
    stop(sprintf("%s must be a comma- or space-separated numeric vector.",
                 name), call. = FALSE)
  }
  as.numeric(out)
}

override_movement_settings <- function(settings) {
  mu_override <- Sys.getenv("BUSIA_MOVEMENT_MU", unset = NA_character_)
  if (!is.na(mu_override) && nzchar(mu_override)) {
    settings$mu <- parse_numeric_env("BUSIA_MOVEMENT_MU", settings$mu,
                                     positive = TRUE)
  }
  p_move_override <- Sys.getenv("BUSIA_P_MOVE", unset = NA_character_)
  if (!is.na(p_move_override) && nzchar(p_move_override)) {
    settings$p_move <- parse_numeric_env("BUSIA_P_MOVE", settings$p_move,
                                         positive = TRUE)
  }
  if (!is.numeric(settings$p_move) || length(settings$p_move) != 1L ||
      !is.finite(settings$p_move) || settings$p_move <= 0 ||
      settings$p_move >= 1) {
    stop("Movement p_move must be a finite number in (0, 1).", call. = FALSE)
  }
  settings
}

validate_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("%s must be a finite number in [0, 1].", name),
         call. = FALSE)
  }
  x
}

safe_sd <- function(x) if (length(x) < 2L) NA_real_ else stats::sd(x)

safe_mean_na <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 1L || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_min_na <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 1L || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max_na <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 1L || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
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

summarise_incidence_run <- function(data_list, nodes, covariate,
                                    measurement_window_days, age_band_days,
                                    init_EIR, targets) {
  if (length(data_list) != nrow(nodes)) {
    stop(sprintf("Run returned %d node outputs, expected %d.",
                 length(data_list), nrow(nodes)), call. = FALSE)
  }
  inc_col <- sprintf("n_inc_%d_%d", age_band_days[[1]], age_band_days[[2]])
  age_col <- sprintf("n_age_%d_%d", age_band_days[[1]], age_band_days[[2]])
  required <- c("timestep", inc_col, age_col)
  optional <- c("EIR_gamb", "total_M_gamb")

  node_rows <- vector("list", length(data_list))
  unavailable <- character(0)
  for (i in seq_along(data_list)) {
    df <- data_list[[i]]
    missing <- setdiff(required, names(df))
    if (length(missing) > 0L) {
      stop(sprintf("Node %d is missing required columns: %s",
                   i, paste(missing, collapse = ", ")), call. = FALSE)
    }
    max_t <- max(df$timestep)
    start_t <- max_t - measurement_window_days + 1L
    win <- df[df$timestep >= start_t & df$timestep <= max_t, , drop = FALSE]
    if (nrow(win) != measurement_window_days) {
      stop(sprintf("Node %d expected %d measurement days but found %d.",
                   i, measurement_window_days, nrow(win)), call. = FALSE)
    }
    for (col in optional) {
      if (!(col %in% names(df))) {
        unavailable <- c(unavailable, sprintf("node %d: %s unavailable", i, col))
      }
    }

    incident_infections <- sum(win[[inc_col]])
    person_days <- sum(win[[age_col]])
    person_years <- person_days / 365
    incidence <- incident_infections / pmax(person_years, .Machine$double.eps)

    node_rows[[i]] <- data.frame(
      node_id = as.integer(nodes$node[[i]]),
      cluster_id = as.integer(nodes$cluster_id[[i]]),
      human_population = as.integer(nodes$NH_per_node[[i]]),
      population_weight = as.numeric(nodes$NH_per_node[[i]]) / sum(nodes$NH_per_node),
      covariate_z = as.numeric(covariate$z[[i]]),
      underlying_S = as.numeric(covariate$underlying_S[[i]]),
      contact_multiplier = as.numeric(covariate$contact_multiplier[[i]]),
      n_inc_182_3650 = incident_infections,
      n_age_182_3650 = person_days,
      person_years_at_risk = person_years,
      infection_incidence_per_person_year = incidence,
      difference_from_target_mean = incidence - targets$mean,
      within_reference_range =
        incidence >= targets$range[[1]] && incidence <= targets$range[[2]],
      mean_eir = if ("EIR_gamb" %in% names(win)) mean(win$EIR_gamb) else NA_real_,
      mean_adult_female_mosquitoes = if ("total_M_gamb" %in% names(win)) {
        mean(win$total_M_gamb)
      } else {
        NA_real_
      },
      endpoint = inc_col,
      stringsAsFactors = FALSE
    )
  }

  node_summary <- do.call(rbind, node_rows)
  rates <- node_summary$infection_incidence_per_person_year
  unweighted_mean <- mean(rates)
  unweighted_sd <- safe_sd(rates)
  unweighted_cv <- unweighted_sd / pmax(unweighted_mean, .Machine$double.eps)
  pooled_inc <- sum(node_summary$n_inc_182_3650) /
    pmax(sum(node_summary$n_age_182_3650) / 365, .Machine$double.eps)
  mean_rel_error <- abs(unweighted_mean - targets$mean) / targets$mean
  node_eir <- node_summary$mean_eir
  node_m <- node_summary$mean_adult_female_mosquitoes

  global_summary <- data.frame(
    init_EIR = as.numeric(init_EIR),
    fixed_delta = as.numeric(covariate$settings$delta),
    n_nodes = nrow(node_summary),
    measurement_window_days = as.integer(measurement_window_days),
    incidence_age_min_days = as.integer(age_band_days[[1]]),
    incidence_age_max_days = as.integer(age_band_days[[2]]),
    incidence_endpoint = inc_col,
    incidence_endpoint_definition =
      "modeled infection incidence; not microscopy-confirmed incidence",
    unweighted_mean_incidence_per_person_year = unweighted_mean,
    target_mean_incidence_per_person_year = targets$mean,
    mean_difference_from_target = unweighted_mean - targets$mean,
    mean_relative_error = mean_rel_error,
    site_sd_incidence_per_person_year = unweighted_sd,
    target_site_sd_incidence_per_person_year = targets$sd,
    sd_difference_from_target = unweighted_sd - targets$sd,
    between_cluster_hazard_cv = unweighted_cv,
    target_between_cluster_hazard_cv = targets$cv,
    cv_difference_from_target = unweighted_cv - targets$cv,
    pooled_population_weighted_incidence_per_person_year = pooled_inc,
    min_site_incidence_per_person_year = min(rates),
    max_site_incidence_per_person_year = max(rates),
    reference_range_min_incidence_per_person_year = targets$range[[1]],
    reference_range_max_incidence_per_person_year = targets$range[[2]],
    n_sites_inside_reference_range = sum(node_summary$within_reference_range),
    contact_multiplier_min = min(node_summary$contact_multiplier),
    contact_multiplier_mean = mean(node_summary$contact_multiplier),
    contact_multiplier_max = max(node_summary$contact_multiplier),
    contact_multiplier_sd = safe_sd(node_summary$contact_multiplier),
    mean_realized_eir = safe_mean_na(node_eir),
    min_node_mean_realized_eir = safe_min_na(node_eir),
    max_node_mean_realized_eir = safe_max_na(node_eir),
    mean_realized_adult_female_mosquitoes = safe_mean_na(node_m),
    min_node_mean_adult_female_mosquitoes = safe_min_na(node_m),
    max_node_mean_adult_female_mosquitoes = safe_max_na(node_m),
    unavailable_output_fields = paste(unique(unavailable), collapse = "; "),
    stringsAsFactors = FALSE
  )

  list(node_summary = node_summary, global_summary = global_summary)
}

candidate_score <- function(global_summary, targets) {
  abs(global_summary$unweighted_mean_incidence_per_person_year - targets$mean)
}

add_counterfactual_equilibrium_fields <- function(summary, res, init_EIR) {
  eq_total_M <- attr(res, "equilibrium_total_M")
  init_foim <- attr(res, "init_foim")
  n_nodes <- nrow(summary$node_summary)

  if (is.null(eq_total_M) || length(eq_total_M) != n_nodes) {
    eq_total_M <- rep(NA_real_, n_nodes)
  }
  if (is.null(init_foim) || length(init_foim) != n_nodes) {
    init_foim <- rep(NA_real_, n_nodes)
  }

  summary$node_summary$raw_counterfactual_init_EIR <- as.numeric(init_EIR)
  summary$node_summary$counterfactual_no_intervention_equilibrium_total_M <-
    as.numeric(eq_total_M)
  summary$node_summary$counterfactual_no_intervention_init_foim <-
    as.numeric(init_foim)

  summary$global_summary$raw_counterfactual_init_EIR <- as.numeric(init_EIR)
  summary$global_summary$calibrated_eir_definition <-
    paste(
      "raw counterfactual no-intervention init_EIR passed to set_equilibrium();",
      "baseline ITNs and ACT treatment are applied after equilibrium initialization"
    )
  summary$global_summary$mean_counterfactual_no_intervention_equilibrium_total_M <-
    mean(as.numeric(eq_total_M), na.rm = TRUE)
  summary$global_summary$min_counterfactual_no_intervention_equilibrium_total_M <-
    min(as.numeric(eq_total_M), na.rm = TRUE)
  summary$global_summary$max_counterfactual_no_intervention_equilibrium_total_M <-
    max(as.numeric(eq_total_M), na.rm = TRUE)
  summary$global_summary$mean_counterfactual_no_intervention_init_foim <-
    mean(as.numeric(init_foim), na.rm = TRUE)
  summary
}

apply_busia_baseline_interventions <- function(parameters, interventions) {
  if (isTRUE(interventions$itn_enabled)) {
    n_species <- length(parameters$species)
    parameters <- set_bednets(
      parameters = parameters,
      timesteps = interventions$itn_timesteps,
      coverages = rep(interventions$itn_usage, length(interventions$itn_timesteps)),
      retention = interventions$itn_retention_days,
      dn0 = matrix(interventions$itn_dn0,
                   nrow = length(interventions$itn_timesteps),
                   ncol = n_species),
      rn = matrix(interventions$itn_rn,
                  nrow = length(interventions$itn_timesteps),
                  ncol = n_species),
      rnm = matrix(interventions$itn_rnm,
                   nrow = length(interventions$itn_timesteps),
                   ncol = n_species),
      gamman = rep(interventions$itn_gamman_days,
                   length(interventions$itn_timesteps))
    )
    parameters$initial_bednet_coverage <- interventions$itn_usage
    parameters$initial_bednet_time <- interventions$itn_initial_time
    parameters$bednet_skip_timesteps <- interventions$itn_initial_time
  }

  if (isTRUE(interventions$treatment_enabled)) {
    drug_params <- switch(
      interventions$treatment_drug,
      AL = AL_params,
      DHA_PQP = DHA_PQP_params,
      stop(sprintf("Unsupported BUSIA_TREATMENT_DRUG: %s",
                   interventions$treatment_drug), call. = FALSE)
    )
    parameters <- set_drugs(parameters, drugs = list(drug_params))
    parameters <- set_clinical_treatment(
      parameters = parameters,
      drug = 1L,
      timesteps = interventions$treatment_timesteps,
      coverages = rep(interventions$treatment_coverage,
                      length(interventions$treatment_timesteps))
    )
  }

  parameters
}

run_truth_raw_counterfactual_eir <- function(setup, cube, NF, NH, tmax, mu,
                                             p_move, release, theta, init_EIR,
                                             prevalence_rendering_min_age = NULL,
                                             prevalence_rendering_max_age = NULL,
                                             infection_incidence_min_age = NULL,
                                             infection_incidence_max_age = NULL,
                                             clinical_incidence_min_age = NULL,
                                             clinical_incidence_max_age = NULL,
                                             warmup_days = 0L,
                                             parameter_modifier = NULL,
                                             post_equilibrium_modifier = NULL,
                                             seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  prevalence_band <- .msimGD_normalise_optional_age_band(
    prevalence_rendering_min_age,
    prevalence_rendering_max_age,
    "Prevalence rendering"
  )
  infection_band <- .msimGD_normalise_optional_age_band(
    infection_incidence_min_age,
    infection_incidence_max_age,
    "Infection incidence"
  )
  clinical_band <- .msimGD_normalise_optional_age_band(
    clinical_incidence_min_age,
    clinical_incidence_max_age,
    "Clinical incidence"
  )

  prep <- .msimGD_prepare_truth_run(
    setup = setup,
    cube = cube,
    NF = NF,
    NH = NH,
    tmax = tmax,
    mu = mu,
    p_move = p_move,
    release = release,
    theta = theta,
    init_EIR = init_EIR,
    prevalence_band = prevalence_band,
    infection_band = infection_band,
    clinical_band = clinical_band,
    warmup_days = warmup_days,
    parameter_modifier = parameter_modifier,
    baseline_checkpoint_state = NULL,
    baseline_checkpoint_metadata = NULL,
    human_initialization_library = NULL
  )

  if (!is.null(post_equilibrium_modifier)) {
    for (nd in seq_along(prep$parameters)) {
      prep$parameters[[nd]] <- post_equilibrium_modifier(
        prep$parameters[[nd]],
        node_index = as.integer(nd),
        warmup_days = prep$effective_warmup_days
      )
      if (is.null(prep$parameters[[nd]]) || !is.list(prep$parameters[[nd]])) {
        stop("post_equilibrium_modifier must return a parameter list.",
             call. = FALSE)
      }
    }
    prep$baseline_time_dependent_signature <-
      .msimGD_capture_time_dependent_signature(prep$parameters)
    prep$baseline_contact_signature <-
      .msimGD_capture_contact_signature(prep$parameters)
    prep$contact_multiplier_by_node <-
      .msimGD_contact_multiplier_by_node(prep$parameters)
    prep$contact_covariates_by_node <-
      .msimGD_contact_covariates_by_node(prep$parameters)
  }

  result <- run_metapop_simulation(
    timesteps = prep$total_timesteps,
    parameters = prep$parameters,
    mixing_tt = 1,
    export_mixing = list(diag(prep$n_nodes)),
    import_mixing = list(diag(prep$n_nodes)),
    p_captured_tt = 1,
    p_captured = list(matrix(0, prep$n_nodes, prep$n_nodes)),
    p_success = 0
  )

  attr(result, "equilibrium_total_M") <- prep$equilibrium_total_M
  attr(result, "effective_total_M") <- prep$effective_total_M
  attr(result, "init_foim") <- prep$init_foim
  attr(result, "baseline_time_dependent_signature") <-
    prep$baseline_time_dependent_signature
  attr(result, "baseline_contact_signature") <- prep$baseline_contact_signature
  attr(result, "contact_multiplier_by_node") <- prep$contact_multiplier_by_node
  attr(result, "contact_covariates_by_node") <- prep$contact_covariates_by_node
  attr(result, "prevalence_rendering_min_age") <- prevalence_band$min
  attr(result, "prevalence_rendering_max_age") <- prevalence_band$max
  attr(result, "infection_incidence_min_age") <- infection_band$min
  attr(result, "infection_incidence_max_age") <- infection_band$max
  attr(result, "clinical_incidence_min_age") <- clinical_band$min
  attr(result, "clinical_incidence_max_age") <- clinical_band$max
  attr(result, "warmup_days") <- prep$effective_warmup_days
  attr(result, "raw_counterfactual_init_EIR") <- as.numeric(init_EIR)
  attr(result, "post_equilibrium_interventions") <-
    !is.null(post_equilibrium_modifier)

  .msimGD_crop_result(result, warmup_days = prep$effective_warmup_days)
}

write_outputs <- function(best, out_dir, targets, metadata, search_log) {
  if (is.null(best)) {
    return(invisible(FALSE))
  }
  utils::write.csv(
    best$node_summary,
    file.path(out_dir, "incidence_calibration_fixed_delta_node_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    best$global_summary,
    file.path(out_dir, "incidence_calibration_fixed_delta_global_summary.csv"),
    row.names = FALSE
  )
  if (length(search_log) > 0L) {
    utils::write.csv(
      do.call(rbind, search_log),
      file.path(out_dir, "incidence_calibration_fixed_delta_log.csv"),
      row.names = FALSE
    )
  }
  if (!is.null(best$no_intervention_counterfactual)) {
    utils::write.csv(
      best$no_intervention_counterfactual$node_summary,
      file.path(
        out_dir,
        "incidence_calibration_fixed_delta_no_intervention_node_summary.csv"
      ),
      row.names = FALSE
    )
    utils::write.csv(
      best$no_intervention_counterfactual$global_summary,
      file.path(
        out_dir,
        "incidence_calibration_fixed_delta_no_intervention_global_summary.csv"
      ),
      row.names = FALSE
    )
  }
  saveRDS(
    c(
      list(
        init_EIR = best$init_EIR,
        raw_counterfactual_init_EIR = best$init_EIR,
        calibrated_eir_definition =
          best$global_summary$calibrated_eir_definition,
        fixed_delta = best$delta,
        realised_incidence = best$global_summary$unweighted_mean_incidence_per_person_year,
        target_incidence_mean = targets$mean,
        target_incidence_sd = targets$sd,
        target_hazard_cv = targets$cv,
        reference_incidence_range = targets$range,
        covariate_settings = best$covariate$settings,
        covariate_z = best$covariate$z,
        contact_multiplier = best$covariate$contact_multiplier,
        underlying_S = best$covariate$underlying_S,
        node_summary = best$node_summary,
        global_summary = best$global_summary,
        search_log = if (length(search_log) > 0L) do.call(rbind, search_log) else NULL,
        no_intervention_counterfactual =
          best$no_intervention_counterfactual
      ),
      metadata
    ),
    file.path(out_dir, "calibrated_init_eir_incidence_fixed_delta.rds")
  )
  invisible(TRUE)
}

write_compat_artifact <- function(best, out_dir, theta, metadata, targets) {
  if (is.null(best)) {
    return(invisible(FALSE))
  }
  saveRDS(
    c(
      list(
        init_EIR = best$init_EIR,
        raw_counterfactual_init_EIR = best$init_EIR,
        theta = theta,
        calibration_target = "modeled_infection_incidence_182_3650_fixed_delta",
        calibrated_eir_definition =
          best$global_summary$calibrated_eir_definition,
        target_prevalence = NA_real_,
        target_incidence_mean = targets$mean,
        realised_incidence = best$global_summary$unweighted_mean_incidence_per_person_year,
        fixed_delta = best$delta,
        covariate_settings = best$covariate$settings,
        node_summary = best$node_summary,
        global_summary = best$global_summary
      ),
      metadata
    ),
    file.path(out_dir, "calibrated_init_eir.rds")
  )
  invisible(TRUE)
}

# --- Calibration knobs ------------------------------------------------------

targets <- list(
  mean = 3.138,
  sd = 1.381,
  cv = 0.48,
  range = c(1.1, 6.9)
)
INCIDENCE_AGE_MIN_DAYS <- 182L
INCIDENCE_AGE_MAX_DAYS <- 3650L

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

SEARCH_WARMUP_DAYS <- parse_integer_env(
  "BUSIA_SEARCH_WARMUP_DAYS", 3650L, minimum = 1L
)
MEASUREMENT_WINDOW_DAYS <- min(365L, SEARCH_WARMUP_DAYS)
MAX_REFINE_ITER <- parse_integer_env(
  "BUSIA_CALIBRATION_MAX_REFINE", 3L, minimum = 0L
)
RNG_SEED <- parse_integer_env("BUSIA_RNG_SEED", 20260514L)
MEAN_TOLERANCE_REL <- parse_numeric_env(
  "BUSIA_INCIDENCE_MEAN_TOLERANCE_REL", 0.03, positive = TRUE
)
RUNTIME_BUDGET_MINUTES <- parse_numeric_env(
  "BUSIA_CALIBRATION_RUNTIME_MINUTES", Inf,
  positive = TRUE,
  allow_inf = TRUE
)
RUNTIME_BUDGET_SECONDS <- RUNTIME_BUDGET_MINUTES * 60
WRITE_STAGE00_COMPAT <- parse_logical_env("BUSIA_WRITE_STAGE00_COMPAT", FALSE)
RUN_COUNTERFACTUAL_NO_INTERVENTION <- parse_logical_env(
  "BUSIA_RUN_COUNTERFACTUAL_NO_INTERVENTION",
  TRUE
)

EIR_MIN <- parse_numeric_env("BUSIA_EIR_MIN", 0.01, positive = TRUE)
EIR_MAX <- parse_numeric_env("BUSIA_EIR_MAX", 500, positive = TRUE)
if (EIR_MAX <= EIR_MIN) {
  stop("BUSIA_EIR_MAX must be greater than BUSIA_EIR_MIN.", call. = FALSE)
}
EIR_START <- parse_numeric_env("BUSIA_EIR_START", 30, positive = TRUE)
EIR_EXPANSION_FACTOR <- parse_numeric_env(
  "BUSIA_EIR_EXPANSION_FACTOR", 3, positive = TRUE
)
if (EIR_EXPANSION_FACTOR <= 1) {
  stop("BUSIA_EIR_EXPANSION_FACTOR must be greater than 1.", call. = FALSE)
}
MAX_BRACKET_ITER <- parse_integer_env(
  "BUSIA_EIR_MAX_BRACKET_ITER", 6L, minimum = 1L
)

interventions <- list(
  itn_enabled = !parse_logical_env("BUSIA_DISABLE_ITNS", FALSE),
  itn_ownership = validate_probability(
    parse_numeric_env("BUSIA_ITN_OWNERSHIP", 0.882, non_negative = TRUE),
    "BUSIA_ITN_OWNERSHIP"
  ),
  itn_usage = validate_probability(
    parse_numeric_env("BUSIA_ITN_USAGE", 0.684, non_negative = TRUE),
    "BUSIA_ITN_USAGE"
  ),
  itn_distribution_interval_days = parse_integer_env(
    "BUSIA_ITN_DISTRIBUTION_INTERVAL_DAYS", 3L * 365L, minimum = 1L
  ),
  itn_timesteps = seq.int(
    0L,
    SEARCH_WARMUP_DAYS,
    by = parse_integer_env(
      "BUSIA_ITN_DISTRIBUTION_INTERVAL_DAYS", 3L * 365L, minimum = 1L
    )
  ),
  itn_initial_time = 0L,
  itn_retention_days = parse_numeric_env(
    "BUSIA_ITN_RETENTION_DAYS", 5 * 365, positive = TRUE
  ),
  itn_dn0 = validate_probability(
    parse_numeric_env("BUSIA_ITN_DN0", 0.533, non_negative = TRUE),
    "BUSIA_ITN_DN0"
  ),
  itn_rn = validate_probability(
    parse_numeric_env("BUSIA_ITN_RN", 0.56, non_negative = TRUE),
    "BUSIA_ITN_RN"
  ),
  itn_rnm = validate_probability(
    parse_numeric_env("BUSIA_ITN_RNM", 0.24, non_negative = TRUE),
    "BUSIA_ITN_RNM"
  ),
  itn_gamman_days = parse_numeric_env(
    "BUSIA_ITN_GAMMAN_DAYS", 2.64 * 365, positive = TRUE
  ),
  treatment_enabled = !parse_logical_env("BUSIA_DISABLE_TREATMENT", FALSE),
  treatment_coverage = validate_probability(
    parse_numeric_env("BUSIA_TREATMENT_COVERAGE", 0.80, non_negative = TRUE),
    "BUSIA_TREATMENT_COVERAGE"
  ),
  treatment_drug = toupper(Sys.getenv("BUSIA_TREATMENT_DRUG", unset = "AL")),
  treatment_timesteps = 0L
)

cat(sprintf(
  "[%s] careful/00 fixed-delta incidence: target mean %.3f, age %d-%d days\n",
  format(Sys.time(), "%F %T"),
  targets$mean, INCIDENCE_AGE_MIN_DAYS, INCIDENCE_AGE_MAX_DAYS
))
cat(sprintf(
  "Endpoint: n_inc_%d_%d / (n_age_%d_%d / 365); modeled infection incidence, not microscopy-confirmed incidence.\n",
  INCIDENCE_AGE_MIN_DAYS, INCIDENCE_AGE_MAX_DAYS,
  INCIDENCE_AGE_MIN_DAYS, INCIDENCE_AGE_MAX_DAYS
))

# --- Shared Busia setup -----------------------------------------------------

land <- build_busia_landscape()
validate_busia_landscape(land)
NH_per_node <- as.integer(land$nodes$NH_per_node)
validate_busia_demography_for_NH(NH_per_node)

setup <- list(D = land$D, allowed = busia_allowed_matrix(land$n_nodes))
mv_settings <- override_movement_settings(seven_node_movement_settings())
rng <- mu_feasible_range(
  D = setup$D,
  allowed = setup$allowed,
  attractiveness = rep(1, land$n_nodes),
  move_rate = 1
)
if (mv_settings$mu < rng$mu_min || mv_settings$mu > rng$mu_beta0) {
  stop(sprintf(
    paste(
      "Movement mu=%.4f is outside the Busia 58-node feasible range",
      "[%.4f, %.4f]. Edit config/movement.R or set BUSIA_MOVEMENT_MU",
      "to a value inside that range."
    ),
    mv_settings$mu, rng$mu_min, rng$mu_beta0
  ), call. = FALSE)
}

seas <- seven_node_seasonality()
cube <- build_seven_node_drive_cube()
covariate_settings <- seven_node_covariate_settings()
covariate_settings$delta <- parse_numeric_env(
  "BUSIA_FIXED_DELTA",
  covariate_settings$delta,
  non_negative = TRUE
)
cov <- build_seven_node_covariate(setup$D, settings = covariate_settings)

contact_surface <- list(
  type = "contact_surface",
  contact_multiplier = stats::setNames(
    as.numeric(cov$contact_multiplier),
    as.character(seq_len(land$n_nodes))
  )
)

equilibrium_parameter_modifier <- function(parameters, node_index, warmup_days) {
  parameters <- apply_busia_demography(
    parameters,
    total_population = parameters$human_population
  )
  parameters$model_seasonality <- TRUE
  parameters$g0 <- seas$g0
  parameters$g <- seas$g
  parameters$h <- seas$h
  parameters$rainfall_floor <- seas$rainfall_floor
  parameters <- malariasimulationGD::apply_node_contact_surface(
    parameters = parameters,
    contact_surface = contact_surface,
    node_index = as.integer(node_index)
  )
  parameters
}

post_equilibrium_intervention_modifier <- function(parameters, node_index,
                                                  warmup_days) {
  parameters <- apply_busia_baseline_interventions(parameters, interventions)
  parameters
}

metadata <- list(
  theta = theta,
  incidence_age_min_days = INCIDENCE_AGE_MIN_DAYS,
  incidence_age_max_days = INCIDENCE_AGE_MAX_DAYS,
  incidence_endpoint = sprintf("n_inc_%d_%d",
                               INCIDENCE_AGE_MIN_DAYS,
                               INCIDENCE_AGE_MAX_DAYS),
  incidence_age_definition =
    "children aged 6 months to younger than 10 years, rendered as 182-3650 days",
  incidence_endpoint_definition =
    "modeled infection incidence; not microscopy-confirmed incidence",
  search_warmup_days = SEARCH_WARMUP_DAYS,
  measurement_window_days = MEASUREMENT_WINDOW_DAYS,
  rng_seed = RNG_SEED,
  NH_per_node = NH_per_node,
  population_total = sum(NH_per_node),
  population_mean = mean(NH_per_node),
  population_sd = stats::sd(NH_per_node),
  busia_age_structure = read_busia_age_structure(),
  busia_demography = busia_demography_spec(total_population = sum(NH_per_node)),
  fixed_delta = covariate_settings$delta,
  covariate_settings = covariate_settings,
  eir_start = EIR_START,
  eir_expansion_factor = EIR_EXPANSION_FACTOR,
  max_bracket_iter = MAX_BRACKET_ITER,
  max_refine_iter = MAX_REFINE_ITER,
  runtime_budget_minutes = RUNTIME_BUDGET_MINUTES,
  run_counterfactual_no_intervention = RUN_COUNTERFACTUAL_NO_INTERVENTION,
  mean_tolerance_rel = MEAN_TOLERANCE_REL,
  movement_settings = mv_settings,
  movement_feasible_range = rng,
  baseline_interventions = interventions,
  calibrated_eir_definition =
    paste(
      "raw counterfactual no-intervention init_EIR passed to set_equilibrium();",
      "baseline ITNs and ACT treatment are applied after equilibrium initialization"
    )
)

candidate_log <- list()
best <- NULL
budget_limited <- FALSE
t0_script <- Sys.time()

elapsed_seconds <- function() {
  as.numeric(difftime(Sys.time(), t0_script, units = "secs"))
}

can_start_candidate <- function() {
  if (!is.finite(RUNTIME_BUDGET_SECONDS)) {
    return(TRUE)
  }
  elapsed_seconds() < RUNTIME_BUDGET_SECONDS
}

candidate_already_run <- function(init_EIR) {
  if (length(candidate_log) < 1L) {
    return(FALSE)
  }
  log_df <- do.call(rbind, candidate_log)
  any(abs(log_df$init_EIR - init_EIR) / pmax(init_EIR, 1) < 1e-8)
}

current_log <- function() {
  if (length(candidate_log) < 1L) {
    return(NULL)
  }
  do.call(rbind, candidate_log)
}

find_bracket <- function(rows) {
  if (is.null(rows) || nrow(rows) < 2L) {
    return(NULL)
  }
  rows <- rows[order(rows$init_EIR), , drop = FALSE]
  means <- rows$unweighted_mean_incidence_per_person_year
  below_idx <- which(means <= targets$mean)
  above_idx <- which(means >= targets$mean)
  if (length(below_idx) < 1L || length(above_idx) < 1L) {
    return(NULL)
  }
  below <- rows[below_idx[which.max(means[below_idx])], , drop = FALSE]
  above <- rows[above_idx[which.min(means[above_idx])], , drop = FALSE]
  list(below = below, above = above)
}

suggest_next_eir <- function(rows) {
  rows <- rows[order(rows$init_EIR), , drop = FALSE]
  means <- rows$unweighted_mean_incidence_per_person_year
  eirs <- rows$init_EIR

  bracket <- find_bracket(rows)
  if (!is.null(bracket)) {
    below_row <- bracket$below
    above_row <- bracket$above
    if (abs(below_row$init_EIR - above_row$init_EIR) < 1e-12) {
      return(as.numeric(below_row$init_EIR))
    }
    x <- c(below_row$unweighted_mean_incidence_per_person_year,
           above_row$unweighted_mean_incidence_per_person_year)
    y <- log(c(below_row$init_EIR, above_row$init_EIR))
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    if (length(unique(x)) >= 2L) {
      return(as.numeric(exp(stats::approx(x = x, y = y,
                                          xout = targets$mean,
                                          rule = 2)$y)))
    }
  }

  if (targets$mean > max(means)) {
    return(min(max(eirs) * EIR_EXPANSION_FACTOR, EIR_MAX))
  }
  max(min(eirs) / EIR_EXPANSION_FACTOR, EIR_MIN)
}

run_candidate <- function(init_EIR, kind) {
  if (!can_start_candidate()) {
    budget_limited <<- TRUE
    return(NULL)
  }

  cat(sprintf("  [%s] fixed_delta=%5.3f raw_init_EIR=%8.4f ...\n",
              kind, covariate_settings$delta, init_EIR))
  t0 <- Sys.time()
  res <- run_truth_raw_counterfactual_eir(
    setup = setup,
    cube = cube,
    NF = NULL,
    NH = NH_per_node,
    tmax = SEARCH_WARMUP_DAYS,
    mu = mv_settings$mu,
    p_move = mv_settings$p_move,
    release = NULL,
    theta = theta,
    init_EIR = init_EIR,
    prevalence_rendering_min_age = NULL,
    prevalence_rendering_max_age = NULL,
    infection_incidence_min_age = INCIDENCE_AGE_MIN_DAYS,
    infection_incidence_max_age = INCIDENCE_AGE_MAX_DAYS,
    clinical_incidence_min_age = NULL,
    clinical_incidence_max_age = NULL,
    warmup_days = 0L,
    parameter_modifier = equilibrium_parameter_modifier,
    post_equilibrium_modifier = post_equilibrium_intervention_modifier,
    seed = RNG_SEED
  )
  run_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  summary <- summarise_incidence_run(
    data_list = extract_data_list(res),
    nodes = land$nodes,
    covariate = cov,
    measurement_window_days = MEASUREMENT_WINDOW_DAYS,
    age_band_days = c(INCIDENCE_AGE_MIN_DAYS, INCIDENCE_AGE_MAX_DAYS),
    init_EIR = init_EIR,
    targets = targets
  )
  summary <- add_counterfactual_equilibrium_fields(
    summary = summary,
    res = res,
    init_EIR = init_EIR
  )
  global <- summary$global_summary
  global$kind <- kind
  global$post_equilibrium_interventions <- TRUE
  global$itn_enabled <- interventions$itn_enabled
  global$itn_ownership <- interventions$itn_ownership
  global$itn_usage <- interventions$itn_usage
  global$treatment_enabled <- interventions$treatment_enabled
  global$treatment_drug <- interventions$treatment_drug
  global$treatment_coverage <- interventions$treatment_coverage
  global$run_seconds <- run_seconds
  global$elapsed_seconds <- elapsed_seconds()
  summary$global_summary <- global
  candidate_log[[length(candidate_log) + 1L]] <<- global

  score <- candidate_score(global, targets)
  if (is.null(best) || score < candidate_score(best$global_summary, targets)) {
    best <<- list(
      delta = as.numeric(covariate_settings$delta),
      init_EIR = as.numeric(init_EIR),
      covariate = cov,
      node_summary = summary$node_summary,
      global_summary = global
    )
    write_outputs(
      best = best,
      out_dir = out_dir,
      targets = targets,
      metadata = c(metadata, list(budget_limited = budget_limited)),
      search_log = candidate_log
    )
  } else {
    utils::write.csv(
      do.call(rbind, candidate_log),
      file.path(out_dir, "incidence_calibration_fixed_delta_log.csv"),
      row.names = FALSE
    )
  }

  cat(sprintf(
    "       mean=%.3f rel_err=%.3f sd=%.3f cv=%.3f range=[%.3f, %.3f] [%.1f s]\n",
    global$unweighted_mean_incidence_per_person_year,
    global$mean_relative_error,
    global$site_sd_incidence_per_person_year,
    global$between_cluster_hazard_cv,
    global$min_site_incidence_per_person_year,
    global$max_site_incidence_per_person_year,
    run_seconds
  ))
  invisible(summary)
}

run_no_intervention_counterfactual <- function(best) {
  if (is.null(best) || !can_start_candidate()) {
    return(NULL)
  }

  cat(sprintf(
    "  [selected_no_intervention_counterfactual] fixed_delta=%5.3f raw_init_EIR=%8.4f ...\n",
    best$delta,
    best$init_EIR
  ))
  t0 <- Sys.time()
  res <- run_truth_raw_counterfactual_eir(
    setup = setup,
    cube = cube,
    NF = NULL,
    NH = NH_per_node,
    tmax = SEARCH_WARMUP_DAYS,
    mu = mv_settings$mu,
    p_move = mv_settings$p_move,
    release = NULL,
    theta = theta,
    init_EIR = best$init_EIR,
    prevalence_rendering_min_age = NULL,
    prevalence_rendering_max_age = NULL,
    infection_incidence_min_age = INCIDENCE_AGE_MIN_DAYS,
    infection_incidence_max_age = INCIDENCE_AGE_MAX_DAYS,
    clinical_incidence_min_age = NULL,
    clinical_incidence_max_age = NULL,
    warmup_days = 0L,
    parameter_modifier = equilibrium_parameter_modifier,
    post_equilibrium_modifier = NULL,
    seed = RNG_SEED
  )
  run_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  summary <- summarise_incidence_run(
    data_list = extract_data_list(res),
    nodes = land$nodes,
    covariate = cov,
    measurement_window_days = MEASUREMENT_WINDOW_DAYS,
    age_band_days = c(INCIDENCE_AGE_MIN_DAYS, INCIDENCE_AGE_MAX_DAYS),
    init_EIR = best$init_EIR,
    targets = targets
  )
  summary <- add_counterfactual_equilibrium_fields(
    summary = summary,
    res = res,
    init_EIR = best$init_EIR
  )
  summary$global_summary$kind <- "selected_no_intervention_counterfactual"
  summary$global_summary$post_equilibrium_interventions <- FALSE
  summary$global_summary$itn_enabled <- FALSE
  summary$global_summary$itn_ownership <- NA_real_
  summary$global_summary$itn_usage <- NA_real_
  summary$global_summary$treatment_enabled <- FALSE
  summary$global_summary$treatment_drug <- NA_character_
  summary$global_summary$treatment_coverage <- NA_real_
  summary$global_summary$run_seconds <- run_seconds
  summary$global_summary$elapsed_seconds <- elapsed_seconds()
  cat(sprintf(
    "       no-intervention mean=%.3f sd=%.3f cv=%.3f range=[%.3f, %.3f] [%.1f s]\n",
    summary$global_summary$unweighted_mean_incidence_per_person_year,
    summary$global_summary$site_sd_incidence_per_person_year,
    summary$global_summary$between_cluster_hazard_cv,
    summary$global_summary$min_site_incidence_per_person_year,
    summary$global_summary$max_site_incidence_per_person_year,
    run_seconds
  ))
  summary
}

# --- EIR-only search --------------------------------------------------------

cat(sprintf("Fixed covariate delta: %.3f\n", covariate_settings$delta))
cat(sprintf(
  "Calibration EIR definition: raw counterfactual no-intervention init_EIR; interventions are applied after equilibrium initialization.\n"
))
cat(sprintf(
  "Runtime interventions: ITN ownership %.1f%%, ITN active usage %.1f%%, ACT(%s) clinical treatment %.1f%%\n",
  100 * interventions$itn_ownership,
  100 * interventions$itn_usage,
  interventions$treatment_drug,
  100 * interventions$treatment_coverage
))

eir0 <- min(max(EIR_START, EIR_MIN), EIR_MAX)
run_candidate(eir0, "start")

for (it in seq_len(MAX_BRACKET_ITER)) {
  if (budget_limited) {
    break
  }
  log_df <- current_log()
  if (!is.null(find_bracket(log_df))) {
    break
  }
  next_eir <- suggest_next_eir(log_df)
  if (!is.finite(next_eir) || next_eir <= 0 || candidate_already_run(next_eir)) {
    break
  }
  run_candidate(next_eir, sprintf("bracket_%d", it))
}

for (it in seq_len(MAX_REFINE_ITER)) {
  if (budget_limited || length(candidate_log) < 1L) {
    break
  }
  log_df <- current_log()
  best_rel_err <- min(log_df$mean_relative_error)
  if (best_rel_err <= MEAN_TOLERANCE_REL) {
    break
  }
  next_eir <- suggest_next_eir(log_df)
  if (!is.finite(next_eir) || next_eir <= 0 || candidate_already_run(next_eir)) {
    break
  }
  run_candidate(next_eir, sprintf("refine_%d", it))
}

if (is.null(best)) {
  stop("No fixed-delta incidence calibration candidates were evaluated.",
       call. = FALSE)
}

if (RUN_COUNTERFACTUAL_NO_INTERVENTION) {
  counterfactual <- run_no_intervention_counterfactual(best)
  if (!is.null(counterfactual)) {
    best$no_intervention_counterfactual <- counterfactual
  } else {
    budget_limited <- TRUE
    cat("Skipped selected no-intervention counterfactual because the runtime budget was exhausted.\n")
  }
}

metadata$elapsed_seconds <- elapsed_seconds()
metadata$budget_limited <- budget_limited
write_outputs(
  best = best,
  out_dir = out_dir,
  targets = targets,
  metadata = metadata,
  search_log = candidate_log
)

if (WRITE_STAGE00_COMPAT) {
  write_compat_artifact(
    best = best,
    out_dir = out_dir,
    theta = theta,
    metadata = metadata,
    targets = targets
  )
}

cat(sprintf(
  "[%s] fixed-delta incidence calibration %s in %.1f s (%.1f min).\n",
  format(Sys.time(), "%F %T"),
  if (budget_limited) "stopped at runtime budget" else "finished",
  metadata$elapsed_seconds,
  metadata$elapsed_seconds / 60
))
cat(sprintf(
  "Selected fixed_delta = %.3f | raw_init_EIR = %.4f | intervention mean = %.3f | rel_err = %.3f | SD = %.3f | CV = %.3f\n",
  best$delta,
  best$init_EIR,
  best$global_summary$unweighted_mean_incidence_per_person_year,
  best$global_summary$mean_relative_error,
  best$global_summary$site_sd_incidence_per_person_year,
  best$global_summary$between_cluster_hazard_cv
))
if (!is.null(best$no_intervention_counterfactual)) {
  cat(sprintf(
    "Selected raw_init_EIR no-intervention counterfactual mean incidence = %.3f (same EIR, no ITNs/ACT).\n",
    best$no_intervention_counterfactual$global_summary$unweighted_mean_incidence_per_person_year
  ))
}
cat(sprintf("Wrote fixed-delta incidence calibration outputs to: %s\n", out_dir))
if (!WRITE_STAGE00_COMPAT) {
  cat("Did not overwrite calibrated_init_eir.rds. Set BUSIA_WRITE_STAGE00_COMPAT=true to opt in.\n")
}
