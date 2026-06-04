#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00_run_identity_sanity.R
# ------------------------------------------------------------------------------
# Scenario 1: Baseline sanity test for explicit human mobility.
#
# Compares:
#   A) no_human_mobility
#   B) identity_human_mobility, where human_move_probs = diag(n_nodes)
#
# Expected behaviour:
#   - identity mobility should create no trips
#   - visitors_present == 0
#   - residents_away == 0
#   - trips_started == 0
#   - epidemiology and mosquito-genetics outputs should be very similar to
#     no mobility, ideally identical if RNG streams are not affected.
#
# Outputs:
#   output/human_mobility_identity_sanity/
#     timeseries_by_arm.csv
#     carrier_frequency_by_arm.csv
#     release_schedule_by_arm.csv
#     summary_by_arm.csv
#     comparison_summary.csv
#     nodes.csv
#     context.rds
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/seven_node_example/human_mobility/00_run_identity_sanity.R",
    mustWork = FALSE
  )
}

# This assumes the script is placed in:
#   test-scripts/seven_node_example/human_mobility/
# and the existing seven_node_example structure is unchanged.
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
pkg_root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)

if (!requireNamespace("malariasimulationGD", quietly = TRUE)) {
  stop(
    "Install malariasimulationGD first, e.g. Rscript test-scripts/install_local.R",
    call. = FALSE
  )
}
if (!requireNamespace("MGDrivE", quietly = TRUE)) {
  stop(
    "Install MGDrivE first, e.g. Rscript test-scripts/install_local.R",
    call. = FALSE
  )
}

# Source the same library/config files as the existing quick example.
source(file.path(example_dir, "lib", "movement_mu.R"))
source(file.path(example_dir, "lib", "synthetic_covariate.R"))
source(file.path(example_dir, "config", "seven_node_landscape.R"))
source(file.path(example_dir, "config", "movement.R"))
source(file.path(example_dir, "config", "seasonality.R"))
source(file.path(example_dir, "config", "covariate.R"))
source(file.path(example_dir, "config", "homing_drive.R"))
source(file.path(example_dir, "config", "trial_design.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

out_dir <- file.path(example_dir, "output", "human_mobility_identity_sanity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

INIT_EIR_QUICK <- 15
RNG_SEED <- 20260514L

cat(sprintf("[%s] Human mobility identity sanity workflow\n", format(Sys.time(), "%F %T")))
cat(sprintf("  pkg root:    %s\n", pkg_root))
cat(sprintf("  example dir: %s\n", example_dir))
cat(sprintf("  output:      %s\n", out_dir))

# ------------------------------------------------------------------------------
# Landscape, movement, covariates, release design
# ------------------------------------------------------------------------------

set.seed(RNG_SEED)

land <- build_seven_node_landscape()
n_nodes <- land$n_nodes
nodes <- land$nodes
D <- land$D

mv <- build_seven_node_movement(D)
cat(sprintf(
  "  mosquito movement: mu_achieved = %.3f km, p_move = %.4f, origin_rate = %.6f\n",
  mv$mu_achieved, mv$p_move, mv$origin_rate
))

cov <- build_seven_node_covariate(D)
cat(sprintf(
  "  contact multipliers: %s\n",
  paste(sprintf("%.2f", cov$contact_multiplier), collapse = ", ")
))

contact_surface <- list(
  type = "contact_surface",
  contact_multiplier = stats::setNames(
    as.numeric(cov$contact_multiplier),
    as.character(seq_len(n_nodes))
  )
)

td <- seven_node_trial_design()
readout_day <- as.integer(td$release_day + td$horizon_day)
TMAX_RUN <- readout_day

release_nodes <- pick_release_nodes_by_spread(nodes, td$n_release_nodes)
cat(sprintf(
  "  release nodes: %s\n",
  paste(release_nodes, collapse = ", ")
))

seas <- seven_node_seasonality()
cube <- build_seven_node_drive_cube()

# Identity matrices for legacy transmission mixing/capture.
# Explicit mobility must not be combined with non-identity import/export mixing.
mixing_identity <- list(diag(n_nodes))
p_captured_zero <- list(matrix(0, n_nodes, n_nodes))

# ------------------------------------------------------------------------------
# Arm definitions
# ------------------------------------------------------------------------------

arms <- list(
  no_human_mobility = list(
    human_mobility_enabled = FALSE,
    human_move_probs = NULL,
    human_trip_duration_type = NULL,
    human_trip_duration_mean = NULL,
    human_mobility_store_diagnostics = FALSE
  ),
  identity_human_mobility = list(
    human_mobility_enabled = TRUE,
    human_move_probs = diag(n_nodes),
    human_trip_duration_type = "fixed",
    human_trip_duration_mean = 1L,
    human_mobility_store_diagnostics = TRUE
  )
)

# ------------------------------------------------------------------------------
# Build parameters
# ------------------------------------------------------------------------------

build_one_node_params <- function(node_idx, arm_cfg) {
  NH_v <- nodes$NH_per_node[node_idx]
  
  p <- malariasimulationGD::get_parameters(list(
    human_population = NH_v,
    individual_mosquitoes = FALSE,
    native_mosquito_backend = TRUE,
    model_seasonality = TRUE,
    g0 = seas$g0,
    g = seas$g,
    h = seas$h,
    rainfall_floor = seas$rainfall_floor,
    progress_bar = FALSE
  ))
  
  p <- malariasimulationGD::apply_node_contact_surface(
    parameters = p,
    contact_surface = contact_surface,
    node_index = as.integer(node_idx)
  )
  
  p <- malariasimulationGD::set_equilibrium(
    p,
    init_EIR = INIT_EIR_QUICK
  )
  
  p$cube <- cube
  
  # Mosquito movement remains unchanged across arms.
  p$move_probs <- mv$move_probs
  p$move_rates <- mv$move_rates
  
  # Human mobility arm.
  if (isTRUE(arm_cfg$human_mobility_enabled)) {
    p$human_mobility_enabled <- TRUE
    p$human_mobility_mode <- "explicit"
    p$human_move_probs <- arm_cfg$human_move_probs
    p$human_trip_duration_type <- arm_cfg$human_trip_duration_type
    p$human_trip_duration_mean <- arm_cfg$human_trip_duration_mean
    p$human_mobility_store_diagnostics <- arm_cfg$human_mobility_store_diagnostics
  } else {
    p$human_mobility_enabled <- FALSE
  }
  
  # Release schedule.
  if (node_idx %in% release_nodes) {
    p <- malariasimulationGD::set_releases(p, list(
      releasesStart = td$release_day,
      releasesNumber = 1L,
      releaseCount = td$release_size,
      releaseSex = "M",
      releaseGenotype = "HH",
      releasesInterval = 0L
    ))
  }
  
  p
}

build_params_list <- function(arm_cfg) {
  lapply(seq_len(n_nodes), build_one_node_params, arm_cfg = arm_cfg)
}

# ------------------------------------------------------------------------------
# Output extraction helpers
# ------------------------------------------------------------------------------

stack_timeseries <- function(sim, arm) {
  data_list <- sim$data
  
  if (is.null(data_list) || !is.list(data_list) || length(data_list) != n_nodes) {
    stop("run_metapop_simulation did not return expected per-node sim$data list.",
         call. = FALSE)
  }
  
  ts_df <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    d <- data_list[[i]]
    d$node <- i
    d$arm <- arm
    d
  }))
  
  front_cols <- intersect(c("arm", "node", "timestep"), names(ts_df))
  ts_df[, c(front_cols, setdiff(names(ts_df), front_cols)), drop = FALSE]
}

extract_carrier_frequency <- function(sim, arm) {
  data_list <- sim$data
  H_genos_definition <- c("WH", "HH", "HR", "HB")
  
  carrier_rows <- vector("list", n_nodes)
  first_H_per_node <- integer(n_nodes)
  
  for (i in seq_len(n_nodes)) {
    female <- attr(data_list[[i]], "mosquito_genotype_counts_female")
    male <- attr(data_list[[i]], "mosquito_genotype_counts_male")
    
    if (is.null(female) || is.null(male)) {
      stop(sprintf(
        "Node %d is missing genotype attributes. Was cube attached?",
        i
      ), call. = FALSE)
    }
    
    H_genos <- intersect(H_genos_definition, colnames(female))
    
    H_carriers_F <- rowSums(female[, H_genos, drop = FALSE])
    H_carriers_M <- rowSums(male[, H_genos, drop = FALSE])
    
    total_F <- rowSums(female)
    total_M <- rowSums(male)
    
    cf <- (H_carriers_F + H_carriers_M) / pmax(total_F + total_M, 1)
    
    carrier_rows[[i]] <- data.frame(
      arm = arm,
      node = i,
      time = seq_along(cf),
      carrier_freq = as.numeric(cf)
    )
    
    fd <- which((H_carriers_F + H_carriers_M) > 0)
    first_H_per_node[i] <- if (length(fd) > 0L) fd[[1L]] else NA_integer_
  }
  
  list(
    carrier_df = do.call(rbind, carrier_rows),
    first_H_per_node = first_H_per_node
  )
}

extract_release_schedule <- function(sim, arm) {
  data_list <- sim$data
  rows <- list()
  
  for (i in seq_len(n_nodes)) {
    rs <- attr(data_list[[i]], "mosquito_release_schedule")
    if (!is.null(rs) && nrow(rs) > 0L) {
      rs$node <- i
      rs$arm <- arm
      rows[[length(rows) + 1L]] <- rs
    }
  }
  
  if (length(rows) == 0L) {
    data.frame(
      arm = arm,
      node = release_nodes,
      timestep = rep(td$release_day, length(release_nodes)),
      count = rep(td$release_size, length(release_nodes)),
      sex = rep("M", length(release_nodes)),
      genotype = rep("HH", length(release_nodes))
    )
  } else {
    out <- do.call(rbind, rows)
    front_cols <- intersect(c("arm", "node", "timestep"), names(out))
    out[, c(front_cols, setdiff(names(out), front_cols)), drop = FALSE]
  }
}

safe_col_sum <- function(df, col) {
  if (col %in% names(df)) {
    sum(df[[col]], na.rm = TRUE)
  } else {
    0
  }
}

regional_prevalence_at <- function(ts_df, day) {
  idx <- ts_df$timestep == day
  num_col <- "n_detect_lm_730_3650"
  den_col <- "n_age_730_3650"
  
  if (!all(c(num_col, den_col) %in% names(ts_df))) {
    return(NA_real_)
  }
  
  sum(ts_df[[num_col]][idx], na.rm = TRUE) /
    pmax(sum(ts_df[[den_col]][idx], na.rm = TRUE), 1)
}

regional_infections_total <- function(ts_df) {
  safe_col_sum(ts_df, "n_infections")
}

summarise_arm <- function(arm, ts_df, carrier_info, elapsed) {
  carrier_df <- carrier_info$carrier_df
  first_H <- carrier_info$first_H_per_node
  
  first_H_day_overall <- if (any(!is.na(first_H))) {
    min(first_H, na.rm = TRUE)
  } else {
    NA_integer_
  }
  
  peak_carrier_release_nodes <- max(
    carrier_df$carrier_freq[carrier_df$node %in% release_nodes],
    na.rm = TRUE
  )
  
  peak_carrier_nonrelease <- max(
    carrier_df$carrier_freq[!(carrier_df$node %in% release_nodes)],
    na.rm = TRUE
  )
  
  data.frame(
    arm = arm,
    tmax = TMAX_RUN,
    init_eir = INIT_EIR_QUICK,
    n_nodes = n_nodes,
    n_release_nodes = length(release_nodes),
    release_nodes = paste(release_nodes, collapse = ","),
    release_day = td$release_day,
    horizon_day = td$horizon_day,
    readout_day = readout_day,
    release_size = td$release_size,
    first_H_carrier_day = first_H_day_overall,
    peak_carrier_release = peak_carrier_release_nodes,
    peak_carrier_nonrelease = peak_carrier_nonrelease,
    regional_infections_total = regional_infections_total(ts_df),
    regional_pfpr_readout = regional_prevalence_at(ts_df, readout_day),
    total_humans_present = safe_col_sum(ts_df, "humans_present"),
    total_visitors_present = safe_col_sum(ts_df, "visitors_present"),
    total_residents_away = safe_col_sum(ts_df, "residents_away"),
    total_trips_started = safe_col_sum(ts_df, "trips_started"),
    elapsed_seconds = round(elapsed, 1),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Run one arm
# ------------------------------------------------------------------------------

run_arm <- function(arm, arm_cfg) {
  cat(sprintf("\n[%s] Running arm: %s\n", format(Sys.time(), "%F %T"), arm))
  
  # Rebuild parameters fresh for each arm to avoid accidental mutation.
  params_list <- build_params_list(arm_cfg)
  
  # Reset the seed per arm. If identity mobility consumes no extra RNG, this
  # makes equality checks stronger. If it does consume RNG internally, exact
  # equality may not hold, but mobility diagnostics should still be zero.
  set.seed(RNG_SEED)
  
  t0 <- Sys.time()
  sim <- malariasimulationGD::run_metapop_simulation(
    timesteps = TMAX_RUN,
    parameters = params_list,
    mixing_tt = 1L,
    export_mixing = mixing_identity,
    import_mixing = mixing_identity,
    p_captured_tt = 1L,
    p_captured = p_captured_zero,
    p_success = 0,
    return_state = FALSE,
    render_output = TRUE,
    return_summary = TRUE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  cat(sprintf("[%s] Finished arm %s in %.1f seconds\n",
              format(Sys.time(), "%F %T"), arm, elapsed))
  
  ts_df <- stack_timeseries(sim, arm)
  carrier_info <- extract_carrier_frequency(sim, arm)
  release_schedule <- extract_release_schedule(sim, arm)
  summary_df <- summarise_arm(arm, ts_df, carrier_info, elapsed)
  
  list(
    sim = sim,
    timeseries = ts_df,
    carrier_frequency = carrier_info$carrier_df,
    release_schedule = release_schedule,
    summary = summary_df
  )
}

# ------------------------------------------------------------------------------
# Run all arms
# ------------------------------------------------------------------------------

results <- list()

for (arm in names(arms)) {
  results[[arm]] <- run_arm(arm, arms[[arm]])
}

bind_rows_fill <- function(dfs) {
  all_names <- unique(unlist(lapply(dfs, names)))
  
  dfs <- lapply(dfs, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) {
      d[[nm]] <- NA
    }
    d[, all_names, drop = FALSE]
  })
  
  do.call(rbind, dfs)
}


timeseries_all <- bind_rows_fill(lapply(results, `[[`, "timeseries"))
carrier_all <- bind_rows_fill(lapply(results, `[[`, "carrier_frequency"))
release_all <- bind_rows_fill(lapply(results, `[[`, "release_schedule"))
summary_all <- bind_rows_fill(lapply(results, `[[`, "summary"))
# ------------------------------------------------------------------------------
# Comparison checks
# ------------------------------------------------------------------------------

compare_arms <- function(ts_all, arm_a, arm_b) {
  a <- ts_all[ts_all$arm == arm_a, , drop = FALSE]
  b <- ts_all[ts_all$arm == arm_b, , drop = FALSE]
  
  key <- c("node", "timestep")
  
  common_cols <- intersect(names(a), names(b))
  numeric_cols <- common_cols[vapply(a[common_cols], is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, key)
  
  m <- merge(
    a[, c(key, numeric_cols), drop = FALSE],
    b[, c(key, numeric_cols), drop = FALSE],
    by = key,
    suffixes = c("_a", "_b")
  )
  
  rows <- lapply(numeric_cols, function(col) {
    x <- m[[paste0(col, "_a")]]
    y <- m[[paste0(col, "_b")]]
    
    paired <- is.finite(x) & is.finite(y)
    
    if (!any(paired)) {
      return(NULL)
    }
    
    diff <- abs(x[paired] - y[paired])
    
    data.frame(
      variable = col,
      n_paired = sum(paired),
      max_abs_diff = max(diff),
      mean_abs_diff = mean(diff),
      stringsAsFactors = FALSE
    )
  })
  
  rows <- Filter(Negate(is.null), rows)
  
  if (length(rows) == 0L) {
    return(data.frame(
      variable = character(),
      n_paired = integer(),
      max_abs_diff = numeric(),
      mean_abs_diff = numeric()
    ))
  }
  
  do.call(rbind, rows)
}

comparison_summary <- compare_arms(
  timeseries_all,
  "no_human_mobility",
  "identity_human_mobility"
)

# Put most relevant diagnostics first.
priority_vars <- c(
  "n_infections",
  "n_detect_lm_730_3650",
  "n_age_730_3650",
  "humans_present",
  "visitors_present",
  "residents_away",
  "trips_started"
)
comparison_summary$priority <- match(comparison_summary$variable, priority_vars)
comparison_summary$priority[is.na(comparison_summary$priority)] <- 999L
comparison_summary <- comparison_summary[order(
  comparison_summary$priority,
  comparison_summary$variable
), ]
comparison_summary$priority <- NULL

# Expected strict mobility diagnostic checks for identity mobility.
identity_ts <- timeseries_all[timeseries_all$arm == "identity_human_mobility", ]

identity_diag <- data.frame(
  check = c(
    "sum_visitors_present",
    "sum_residents_away",
    "sum_trips_started"
  ),
  value = c(
    safe_col_sum(identity_ts, "visitors_present"),
    safe_col_sum(identity_ts, "residents_away"),
    safe_col_sum(identity_ts, "trips_started")
  )
)

cat("\nIdentity mobility diagnostics:\n")
print(identity_diag)

if (any(identity_diag$value != 0)) {
  warning(
    "Identity mobility produced nonzero visitors/residents_away/trips_started. ",
    "This should be investigated."
  )
}

cat("\nTop comparison differences: no mobility vs identity mobility\n")
print(utils::head(comparison_summary, 20))

# ------------------------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------------------------

utils::write.csv(
  timeseries_all,
  file.path(out_dir, "timeseries_by_arm.csv"),
  row.names = FALSE
)

utils::write.csv(
  carrier_all,
  file.path(out_dir, "carrier_frequency_by_arm.csv"),
  row.names = FALSE
)

utils::write.csv(
  release_all,
  file.path(out_dir, "release_schedule_by_arm.csv"),
  row.names = FALSE
)

utils::write.csv(
  summary_all,
  file.path(out_dir, "summary_by_arm.csv"),
  row.names = FALSE
)

utils::write.csv(
  comparison_summary,
  file.path(out_dir, "comparison_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  identity_diag,
  file.path(out_dir, "identity_mobility_diagnostics.csv"),
  row.names = FALSE
)

utils::write.csv(
  nodes,
  file.path(out_dir, "nodes.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    scenario = "baseline_sanity_no_mobility_vs_identity_mobility",
    arms = arms,
    release_nodes = release_nodes,
    release_day = td$release_day,
    horizon_day = td$horizon_day,
    readout_day = readout_day,
    contact_multiplier = cov$contact_multiplier,
    mosquito_movement = list(
      mu_achieved = mv$mu_achieved,
      mu_min = mv$mu_min,
      mu_beta0 = mv$mu_beta0,
      beta = mv$beta,
      p_move = mv$p_move,
      origin_rate = mv$origin_rate
    ),
    human_mobility_identity = diag(n_nodes),
    identity_mobility_diagnostics = identity_diag
  ),
  file.path(out_dir, "context.rds")
)

cat(sprintf("\n[%s] Done. Outputs written to: %s\n",
            format(Sys.time(), "%F %T"), out_dir))

cat("\nSummary by arm:\n")
print(summary_all)

cat("\nIdentity mobility should have zero visitors, residents away, and trips started.\n")
print(identity_diag)

