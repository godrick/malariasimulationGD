#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00_run_quick.R
# ------------------------------------------------------------------------------
# Quick workflow: a single-shot 7-node metapopulation release simulation.
#
# Skips calibration and warmup. Uses a hardcoded init_EIR (set below) and
# starts from set_equilibrium()'s analytic equilibrium. The release fires on
# release_day; the run continues to release_day + horizon_day.
#
# Outputs (under output/quick/):
#   timeseries.csv          per-day, per-node simulation summary
#   carrier_frequency.csv   per-day, per-node adult mosquito drive-carrier frequency
#   release_schedule.csv    the schedule the simulator actually consumed
#   summary.csv             one-row sanity summary
# ------------------------------------------------------------------------------

# Resolve absolute paths regardless of where Rscript is invoked from.
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/seven_node_example/quick/00_run_quick.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
pkg_root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)

if (!requireNamespace("malariasimulationGD", quietly = TRUE)) {
  stop("Install malariasimulationGD first (for example: Rscript test-scripts/install_local.R).",
       call. = FALSE)
}
if (!requireNamespace("MGDrivE", quietly = TRUE)) {
  stop("Install MGDrivE first (for example: Rscript test-scripts/install_local.R).",
       call. = FALSE)
}

# Source lib + config (note order: lib first so config can reference helpers).
source(file.path(example_dir, "lib", "movement_mu.R"))
source(file.path(example_dir, "lib", "synthetic_covariate.R"))
source(file.path(example_dir, "config", "seven_node_landscape.R"))
source(file.path(example_dir, "config", "movement.R"))
source(file.path(example_dir, "config", "seasonality.R"))
source(file.path(example_dir, "config", "covariate.R"))
source(file.path(example_dir, "config", "homing_drive.R"))
source(file.path(example_dir, "config", "trial_design.R"))

out_dir <- file.path(example_dir, "output", "quick")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Configuration -----------------------------------------------------------

INIT_EIR_QUICK <- 15            # annual EIR; tuned so prevalence sits ~0.2 with seasonality
RNG_SEED <- 20260514L

cat(sprintf("[%s] 7-node QUICK workflow\n", format(Sys.time(), "%F %T")))
cat(sprintf("  pkg root:       %s\n", pkg_root))
cat(sprintf("  example dir:    %s\n", example_dir))
cat(sprintf("  output:         %s\n", out_dir))

set.seed(RNG_SEED)

# --- Build geometry and movement --------------------------------------------

land <- build_seven_node_landscape()
n_nodes <- land$n_nodes
nodes <- land$nodes
D <- land$D

mv <- build_seven_node_movement(D)
cat(sprintf("  movement: mu_achieved = %.3f km (band [%.2f, %.2f]), p_move = %.4f (origin rate = %.6f)\n",
            mv$mu_achieved, mv$mu_min, mv$mu_beta0, mv$p_move, mv$origin_rate))

cov <- build_seven_node_covariate(D)
cat(sprintf("  contact multipliers: %s\n",
            paste(sprintf("%.2f", cov$contact_multiplier), collapse = ", ")))
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
cat(sprintf("  release nodes (spread-picked): %s\n",
            paste(release_nodes, collapse = ", ")))

seas <- seven_node_seasonality()

cube <- build_seven_node_drive_cube()

# --- Build per-node parameters ----------------------------------------------

# For each of the 7 nodes:
#   - get_parameters() with the node's NH and seasonality
#   - apply the per-node contact multiplier through the package runtime
#     contact-surface API
#   - set_equilibrium() at the hardcoded init_EIR
#   - attach the drive cube, the movement matrix, and the release schedule
build_one_node_params <- function(node_idx) {
  NH_v <- nodes$NH_per_node[node_idx]

  # Mosquito engine: native count-based tau-leap. We explicitly set
  # individual_mosquitoes = FALSE so it is unambiguous that we are not
  # routing through the legacy event-based individual-mosquito backend.
  p <- malariasimulationGD::get_parameters(list(
    human_population = NH_v,
    individual_mosquitoes = FALSE,
    native_mosquito_backend = TRUE,
    model_seasonality = TRUE,
    g0 = seas$g0, g = seas$g, h = seas$h,
    rainfall_floor = seas$rainfall_floor,
    progress_bar = FALSE
  ))

  p <- malariasimulationGD::apply_node_contact_surface(
    parameters = p,
    contact_surface = contact_surface,
    node_index = as.integer(node_idx)
  )

  # Calibrated equilibrium at the chosen init_EIR
  p <- malariasimulationGD::set_equilibrium(p, init_EIR = INIT_EIR_QUICK)

  # Attach the drive cube
  p$cube <- cube

  # Movement: every node gets the same n x n matrix and per-origin rate.
  p$move_probs <- mv$move_probs
  p$move_rates <- mv$move_rates

  # Release schedule: only the release nodes fire. Match the production
  # audit-cell call signature for set_releases() (explicit releaseGenotype
  # and releasesInterval = 0L so it's a single instant release event).
  if (node_idx %in% release_nodes) {
    p <- malariasimulationGD::set_releases(p, list(
      releasesStart    = td$release_day,
      releasesNumber   = 1L,
      releaseCount     = td$release_size,
      releaseSex       = "M",
      releaseGenotype  = "HH",
      releasesInterval = 0L
    ))
  }

  p
}

params_list <- lapply(seq_len(n_nodes), build_one_node_params)

# --- Run the metapop simulation ---------------------------------------------

# Mixing matrices: identity (no cross-village human biting), since the example
# already has between-village MOSQUITO movement carrying the spread. This is
# the cleanest demonstration of the drive sweeping across villages via the
# mosquito-movement matrix alone.
mixing_identity <- list(diag(n_nodes))
p_captured_zero <- list(matrix(0, n_nodes, n_nodes))

cat(sprintf("[%s] running run_metapop_simulation (timesteps = %d)...\n",
            format(Sys.time(), "%F %T"), TMAX_RUN))

t0 <- Sys.time()
sim <- malariasimulationGD::run_metapop_simulation(
  timesteps      = TMAX_RUN,
  parameters     = params_list,
  mixing_tt      = 1L,
  export_mixing  = mixing_identity,
  import_mixing  = mixing_identity,
  p_captured_tt  = 1L,
  p_captured     = p_captured_zero,
  p_success      = 0,
  return_state   = FALSE,
  render_output  = TRUE,
  return_summary = TRUE
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("[%s] simulation finished in %.1f seconds\n",
            format(Sys.time(), "%F %T"), elapsed))

# --- Extract outputs --------------------------------------------------------

# sim$data is a list of length n_nodes, each element a data.frame with the
# per-day epi + mosquito counts for that node. Genotype-resolved adult counts
# are attached as ATTRIBUTES (not nested objects) on each per-node data.frame.
data_list <- sim$data
if (is.null(data_list) || !is.list(data_list) || length(data_list) != n_nodes) {
  stop("run_metapop_simulation did not return a per-node sim$data list of the expected length.",
       call. = FALSE)
}

# Tag each per-node data.frame with its node index and stack.
ts_df <- do.call(rbind, lapply(seq_along(data_list), function(i) {
  d <- data_list[[i]]
  d$node <- i
  d
}))
front_cols <- intersect(c("node", "timestep"), names(ts_df))
ts_df <- ts_df[, c(front_cols, setdiff(names(ts_df), front_cols))]

# Carrier frequency from the per-node genotype-history attributes. The
# MGDrivE::cubeHomingDrive() genotype set is (W = wildtype, H = drive,
# R = resistance, B = broken):
#   WW, WH, WR, WB, HH, HR, HB, RR, RB, BB
# An H-carrier is any genotype with at least one H allele: WH, HH, HR, HB.
carrier_rows <- list()
release_attr_per_node <- vector("list", n_nodes)
first_H_per_node <- integer(n_nodes)
H_genos_definition <- c("WH", "HH", "HR", "HB")
for (i in seq_len(n_nodes)) {
  female <- attr(data_list[[i]], "mosquito_genotype_counts_female")
  male   <- attr(data_list[[i]], "mosquito_genotype_counts_male")
  if (is.null(female) || is.null(male)) {
    stop(sprintf("Node %d is missing genotype attributes. Was cube attached to its parameter list?", i),
         call. = FALSE)
  }
  H_genos <- intersect(H_genos_definition, colnames(female))
  H_carriers_F <- rowSums(female[, H_genos, drop = FALSE])
  H_carriers_M <- rowSums(male  [, H_genos, drop = FALSE])
  total_F <- rowSums(female)
  total_M <- rowSums(male)
  cf <- (H_carriers_F + H_carriers_M) / pmax(total_F + total_M, 1)
  carrier_rows[[i]] <- data.frame(
    node = i, time = seq_along(cf), carrier_freq = as.numeric(cf)
  )
  fd <- which((H_carriers_F + H_carriers_M) > 0)
  first_H_per_node[i] <- if (length(fd) > 0L) fd[[1L]] else NA_integer_
  release_attr_per_node[[i]] <- attr(data_list[[i]], "mosquito_release_schedule")
}
carrier_df <- do.call(rbind, carrier_rows)

# Release schedule: each release node's parameters has its own schedule
# attribute; combine them into a single CSV.
release_schedule_rows <- list()
for (i in seq_along(release_attr_per_node)) {
  rs <- release_attr_per_node[[i]]
  if (!is.null(rs) && nrow(rs) > 0L) {
    rs$node <- i
    release_schedule_rows[[length(release_schedule_rows) + 1L]] <- rs
  }
}
release_schedule <- if (length(release_schedule_rows) > 0L) {
  do.call(rbind, release_schedule_rows)
} else {
  # Fallback so the CSV is never missing.
  data.frame(node = release_nodes,
             timestep = rep(td$release_day, length(release_nodes)),
             count = rep(td$release_size, length(release_nodes)),
             sex = rep("M", length(release_nodes)),
             genotype = rep("HH", length(release_nodes)))
}
front_cols <- intersect(c("node", "timestep"), names(release_schedule))
release_schedule <- release_schedule[, c(front_cols, setdiff(names(release_schedule), front_cols))]

first_H_day_overall <- if (any(!is.na(first_H_per_node))) {
  min(first_H_per_node, na.rm = TRUE)
} else NA_integer_

peak_carrier_release_nodes <- max(
  carrier_df$carrier_freq[carrier_df$node %in% release_nodes],
  na.rm = TRUE
)
peak_carrier_nonrelease <- max(
  carrier_df$carrier_freq[!(carrier_df$node %in% release_nodes)],
  na.rm = TRUE
)

summary_df <- data.frame(
  workflow              = "quick",
  tmax                  = TMAX_RUN,
  init_eir              = INIT_EIR_QUICK,
  n_nodes               = n_nodes,
  n_release_nodes       = length(release_nodes),
  release_nodes         = paste(release_nodes, collapse = ","),
  release_day           = td$release_day,
  horizon_day           = td$horizon_day,
  readout_day           = readout_day,
  release_size          = td$release_size,
  first_H_carrier_day   = first_H_day_overall,
  peak_carrier_release  = peak_carrier_release_nodes,
  peak_carrier_nonrelease = peak_carrier_nonrelease,
  elapsed_seconds       = round(elapsed, 1),
  stringsAsFactors = FALSE
)

# --- Write CSVs --------------------------------------------------------------

utils::write.csv(ts_df,            file.path(out_dir, "timeseries.csv"),         row.names = FALSE)
utils::write.csv(carrier_df,       file.path(out_dir, "carrier_frequency.csv"),  row.names = FALSE)
utils::write.csv(release_schedule, file.path(out_dir, "release_schedule.csv"),   row.names = FALSE)
utils::write.csv(summary_df,       file.path(out_dir, "summary.csv"),            row.names = FALSE)

# Also save the geometry / movement / covariate inputs so the plot script
# doesn't need to re-source the config.
utils::write.csv(nodes, file.path(out_dir, "nodes.csv"), row.names = FALSE)
saveRDS(list(
  release_nodes = release_nodes,
  release_day = td$release_day,
  horizon_day = td$horizon_day,
  readout_day = readout_day,
  contact_multiplier = cov$contact_multiplier,
  movement = list(mu_achieved = mv$mu_achieved,
                  mu_min = mv$mu_min, mu_beta0 = mv$mu_beta0,
                  beta = mv$beta, p_move = mv$p_move, origin_rate = mv$origin_rate)
), file.path(out_dir, "context.rds"))

cat(sprintf("[%s] Done. Outputs in: %s\n",
            format(Sys.time(), "%F %T"), out_dir))
print(summary_df)
