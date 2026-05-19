#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 01_warmup_and_checkpoint.R
# ------------------------------------------------------------------------------
# Careful-workflow stage 1: build a baseline checkpoint LIBRARY (6 snapshots,
# 365 days apart, after a 6-year burnin), score each snapshot's stationarity
# (Y/Y rel_change + cycle RMSE on PfPR2-10 and total_M), promote the
# best-scoring snapshot, and confirm seed-robustness with a multi-seed rerun.
#
# Architecture mirrors the production audit-cell stage-3 + stage-4 pipeline:
# init_EIR grid (handled in 00_) -> checkpoint library -> stationarity battery
# -> promotion -> multi-seed validation. It keeps the same defaults and
# diagnostic signals, with example-local tolerance choices for a compact
# single-replicate share example. The implementation lives under
# lib/baseline_library.R rather than vendoring the full helper stack.
#
# Output (under output/careful/):
#   baseline_checkpoint_library.rds    full library object (all snapshots,
#                                       per-snapshot scores, promotion record)
#   baseline_checkpoint.rds            promoted snapshot only (what 02_ uses)
#   stationarity_battery_summary.csv   per-snapshot metrics + pass/fail
#   multi_seed_validation.rds          across-seed seed-robustness check
#   multi_seed_summary.csv             same summary as CSV
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/seven_node_example/careful/01_warmup_and_checkpoint.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(malariasimulationGD))

source(file.path(example_dir, "lib", "movement_mu.R"))
source(file.path(example_dir, "lib", "calibrate_eir.R"))
source(file.path(example_dir, "lib", "msimGD_truth_generation.R"))
source(file.path(example_dir, "lib", "synthetic_covariate.R"))
source(file.path(example_dir, "lib", "baseline_library.R"))
source(file.path(example_dir, "lib", "multi_seed_validation.R"))

source(file.path(example_dir, "config", "seven_node_landscape.R"))
source(file.path(example_dir, "config", "movement.R"))
source(file.path(example_dir, "config", "seasonality.R"))
source(file.path(example_dir, "config", "covariate.R"))
source(file.path(example_dir, "config", "homing_drive.R"))
source(file.path(example_dir, "config", "trial_design.R"))

out_dir <- file.path(example_dir, "output", "careful")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Knobs (production defaults) -------------------------------------------

BURNIN_TIMESTEPS  <- 2190L      # 6 years -- audit-cell convention
N_SNAPSHOTS       <- 6L         # checkpoint_n_snapshots
SNAPSHOT_SPACING  <- 365L       # checkpoint_snapshot_spacing
MIN_PROMOTION     <- 3L         # min snapshots that must PASS for promotion quorum
VALIDATION_YEARS  <- 2L
N_SEEDS           <- 3L         # multi-seed validation
RNG_SEED          <- 20260514L

# --- Pull calibrated init_EIR + theta from stage 0 -------------------------

cal_path <- file.path(out_dir, "calibrated_init_eir.rds")
if (!file.exists(cal_path)) {
  stop("Run 00_calibrate_init_eir.R first; calibrated_init_eir.rds is missing.",
       call. = FALSE)
}
cal <- readRDS(cal_path)
init_EIR_cal <- cal$init_EIR
theta        <- cal$theta

cat(sprintf("[%s] careful/01: baseline library at init_EIR = %.4f (target PfPR = %.3f)\n",
            format(Sys.time(), "%F %T"), init_EIR_cal, cal$target_prevalence))

# --- Build setup, contact surface, parameter_modifier ----------------------

land  <- build_seven_node_landscape()
n_nodes <- land$n_nodes
nodes <- land$nodes
NH_per_node <- as.integer(nodes$NH_per_node)

allowed_mat <- matrix(TRUE, n_nodes, n_nodes); diag(allowed_mat) <- FALSE
setup <- list(D = land$D, allowed = allowed_mat)

mv_settings <- seven_node_movement_settings()
mu_baseline     <- mv_settings$mu
p_move_baseline <- mv_settings$p_move

rng <- mu_feasible_range(D = setup$D, allowed = setup$allowed,
                         attractiveness = rep(1, n_nodes), move_rate = 1)
if (mu_baseline < rng$mu_min || mu_baseline > rng$mu_beta0) {
  stop(sprintf(
    "Preset mu=%.4f is outside the landscape feasible range [%.4f, %.4f].",
    mu_baseline, rng$mu_min, rng$mu_beta0), call. = FALSE)
}

seas <- seven_node_seasonality()
cov  <- build_seven_node_covariate(setup$D)
cube <- build_seven_node_drive_cube()
td   <- seven_node_trial_design()
readout_day <- as.integer(td$release_day + td$horizon_day)

contact_surface <- list(
  type = "contact_surface",
  contact_multiplier = stats::setNames(
    as.numeric(cov$contact_multiplier),
    as.character(seq_len(n_nodes))
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

# --- 1. Build the snapshot library ----------------------------------------

t0_full <- Sys.time()
library_obj <- build_baseline_checkpoint_library(
  burnin_timesteps = BURNIN_TIMESTEPS,
  n_snapshots      = N_SNAPSHOTS,
  snapshot_spacing = SNAPSHOT_SPACING,
  setup = setup, cube = cube, NH = NH_per_node,
  mu = mu_baseline, p_move = p_move_baseline,
  theta = theta, init_EIR = init_EIR_cal,
  parameter_modifier = parameter_modifier,
  seed = RNG_SEED, verbose = TRUE
)
saveRDS(library_obj, file.path(out_dir, "baseline_checkpoint_library.rds"))

# --- 2. Score each snapshot (stationarity battery) ------------------------

scored <- score_stationarity_per_snapshot(
  library = library_obj,
  setup = setup, cube = cube, NH = NH_per_node,
  mu = mu_baseline, p_move = p_move_baseline,
  theta = theta, init_EIR = init_EIR_cal,
  parameter_modifier = parameter_modifier,
  age_band_days = c(cal$age_min_days, cal$age_max_days),
  validation_years = VALIDATION_YEARS,
  verbose = TRUE
)

# --- 3. Promote the best snapshot -----------------------------------------

promoted <- promote_best_snapshot(
  library_obj, scored, min_promotion_snapshots = MIN_PROMOTION
)
cat(sprintf("\n[%s] PROMOTION: chose snapshot %d (t = %d days; quorum_met=%s; %d passers)\n",
            format(Sys.time(), "%F %T"),
            promoted$promoted_index,
            library_obj$snapshots[[promoted$promoted_index]]$metadata$timesteps_from_t0,
            promoted$quorum_met, promoted$n_passers))

# Persist the promoted checkpoint as the canonical artifact 02_ reads.
saveRDS(promoted$promoted_snapshot,
        file.path(out_dir, "baseline_checkpoint.rds"))
utils::write.csv(promoted$summary,
                 file.path(out_dir, "stationarity_battery_summary.csv"),
                 row.names = FALSE)

# --- 4. Multi-seed validation of the promoted snapshot --------------------

ms <- multi_seed_validate_promoted_snapshot(
  promoted_snapshot = promoted$promoted_snapshot,
  setup = setup, cube = cube, NH = NH_per_node,
  mu = mu_baseline, p_move = p_move_baseline,
  theta = theta, parameter_modifier = parameter_modifier,
  n_seeds = N_SEEDS,
  validation_years = VALIDATION_YEARS,
  age_band_days = c(cal$age_min_days, cal$age_max_days),
  verbose = TRUE
)
saveRDS(ms, file.path(out_dir, "multi_seed_validation.rds"))
utils::write.csv(ms$summary, file.path(out_dir, "multi_seed_summary.csv"),
                 row.names = FALSE)

# --- Persist the context our 02_/03_ scripts will use ---------------------

elapsed_full <- as.numeric(difftime(Sys.time(), t0_full, units = "secs"))
saveRDS(list(
  nodes              = nodes,
  contact_multiplier = cov$contact_multiplier,
  movement_settings  = mv_settings,
  release_nodes      = pick_release_nodes_by_spread(nodes, td$n_release_nodes),
  release_day        = td$release_day,
  horizon_day        = td$horizon_day,
  readout_day        = readout_day,
  burnin_timesteps   = BURNIN_TIMESTEPS,
  n_snapshots        = N_SNAPSHOTS,
  snapshot_spacing   = SNAPSHOT_SPACING,
  promoted_index     = promoted$promoted_index,
  promoted_t_days    = library_obj$snapshots[[promoted$promoted_index]]$metadata$timesteps_from_t0,
  promotion_quorum_met = promoted$quorum_met,
  multi_seed_pass    = ms$pass,
  init_EIR           = init_EIR_cal,
  theta              = theta,
  seasonality        = seas,
  seed               = RNG_SEED,
  elapsed_seconds    = elapsed_full
), file.path(out_dir, "context.rds"))

cat(sprintf("\n[%s] careful/01 complete in %.1f min.\n",
            format(Sys.time(), "%F %T"), elapsed_full / 60))
cat(sprintf("  Saved promoted checkpoint:  %s\n",
            file.path(out_dir, "baseline_checkpoint.rds")))
cat(sprintf("  Saved full library:         %s\n",
            file.path(out_dir, "baseline_checkpoint_library.rds")))
cat(sprintf("  Multi-seed pass:            %s\n", ms$pass))
