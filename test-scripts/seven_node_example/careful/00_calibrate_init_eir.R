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
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/seven_node_example/careful/00_calibrate_init_eir.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(malariasimulationGD))
suppressPackageStartupMessages(library(malariaEquilibrium))

source(file.path(example_dir, "lib", "movement_mu.R"))
source(file.path(example_dir, "lib", "calibrate_eir.R"))
source(file.path(example_dir, "lib", "msimGD_truth_generation.R"))
source(file.path(example_dir, "lib", "synthetic_covariate.R"))
source(file.path(example_dir, "lib", "simulation_calibration.R"))

source(file.path(example_dir, "config", "seven_node_landscape.R"))
source(file.path(example_dir, "config", "movement.R"))
source(file.path(example_dir, "config", "seasonality.R"))
source(file.path(example_dir, "config", "covariate.R"))
source(file.path(example_dir, "config", "homing_drive.R"))

out_dir <- file.path(example_dir, "output", "careful")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Calibration knobs ------------------------------------------------------

TARGET_PFPR    <- 0.30
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
  nu = 1 / (4 / 24)
)
saveRDS(theta, file.path(out_dir, "theta.rds"))

# Short-warmup length used for each calibration probe (years). Production
# uses simulation runs that are similar in scale.
SEARCH_WARMUP_DAYS <- 730L
RNG_SEED <- 20260514L

cat(sprintf("[%s] careful/00: simulation-based calibration to PfPR%d-%d = %.3f\n",
            format(Sys.time(), "%F %T"),
            TARGET_AGE_MIN_YEARS, TARGET_AGE_MAX_YEARS, TARGET_PFPR))

# ============================================================================
# Step 1 — analytical anchor (non-seasonal, fast)
# ============================================================================

land <- build_seven_node_landscape()
NH_per_node <- as.integer(land$nodes$NH_per_node)
NH_mean     <- as.integer(round(mean(NH_per_node)))
seas        <- seven_node_seasonality()

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

allowed_mat <- matrix(TRUE, land$n_nodes, land$n_nodes); diag(allowed_mat) <- FALSE
setup <- list(D = land$D, allowed = allowed_mat)

contact_surface <- list(
  type = "contact_surface",
  contact_multiplier = stats::setNames(
    as.numeric(cov$contact_multiplier),
    as.character(seq_len(land$n_nodes))
  )
)

parameter_modifier <- function(parameters, node_index, warmup_days) {
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
  measurement_window_days = 365L,
  age_band_days      = c(TARGET_AGE_MIN_DAYS, TARGET_AGE_MAX_DAYS),
  starting_init_EIR  = init_EIR_analytical,
  tolerance_rel      = 0.05,
  max_refine_iter    = 3L,
  seed               = RNG_SEED,
  verbose            = TRUE
)
sim_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Step 2 done in %.1f s (%.1f min).\n",
            sim_seconds, sim_seconds / 60))
cat(sprintf("Calibrated init_EIR = %.4f  | realised PfPR = %.4f  (target %.4f)\n",
            cal_sim$init_EIR, cal_sim$realised_pfpr, TARGET_PFPR))

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
    theta                = theta,
    init_EIR_analytical  = init_EIR_analytical,
    analytical_pfpr      = TARGET_PFPR,    # by construction the analytical hit exactly
    search_warmup_days   = SEARCH_WARMUP_DAYS,
    analytical_seconds   = analytical_seconds,
    simulation_seconds   = sim_seconds
  ),
  file.path(out_dir, "calibrated_init_eir.rds")
)
utils::write.csv(cal_sim$search_log,
                 file.path(out_dir, "calibration_log.csv"), row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(out_dir, "calibrated_init_eir.rds")))
cat(sprintf("Saved: %s\n", file.path(out_dir, "calibration_log.csv")))
