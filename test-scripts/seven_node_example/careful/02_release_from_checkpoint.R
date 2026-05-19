#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 02_release_from_checkpoint.R
# ------------------------------------------------------------------------------
# Careful-workflow stage 2: load the baseline checkpoint produced by 01_,
# then call the production helper `msimGD_run_truth()` twice:
#
#   (a) Release arm:    release = list(nodes = release_nodes, time, size,
#                                       stage = "M", genotype = "HH")
#   (b) No-release arm: release = NULL
#
# Both arms pass the SAME `baseline_checkpoint` and SAME `seed`, so they
# share the warmup RNG state and any divergence between them is attributable
# to the release.
#
# Output (under output/careful/):
#   release_truth.rds         msimGD_run_truth output (release arm)
#   no_release_truth.rds      msimGD_run_truth output (no_release arm)
#   release_timeseries.csv    per-day, per-node summary (release arm)
#   no_release_timeseries.csv same (no_release arm)
#   carrier_frequency.csv     release-arm drive carrier frequency
#   prevalence_compare.csv    both arms, regional PfPR(age band) per day
#   incidence_compare.csv     both arms, regional clinical incidence in 14-d windows
#   summary.csv               one-row sanity summary
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/seven_node_example/careful/02_release_from_checkpoint.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(malariasimulationGD))

`%||%` <- function(a, b) if (!is.null(a)) a else b

source(file.path(example_dir, "lib", "movement_mu.R"))
source(file.path(example_dir, "lib", "calibrate_eir.R"))
source(file.path(example_dir, "lib", "msimGD_truth_generation.R"))
source(file.path(example_dir, "lib", "synthetic_covariate.R"))

source(file.path(example_dir, "config", "seven_node_landscape.R"))
source(file.path(example_dir, "config", "movement.R"))
source(file.path(example_dir, "config", "seasonality.R"))
source(file.path(example_dir, "config", "covariate.R"))
source(file.path(example_dir, "config", "homing_drive.R"))
source(file.path(example_dir, "config", "trial_design.R"))

out_dir <- file.path(example_dir, "output", "careful")
ck_path <- file.path(out_dir, "baseline_checkpoint.rds")
ctx_path <- file.path(out_dir, "context.rds")
if (!file.exists(ck_path) || !file.exists(ctx_path)) {
  stop("Run 01_warmup_and_checkpoint.R first.", call. = FALSE)
}
checkpoint <- readRDS(ck_path)
context    <- readRDS(ctx_path)

nodes              <- context$nodes
release_nodes      <- context$release_nodes
release_day        <- context$release_day
horizon_day        <- context$horizon_day
readout_day        <- as.integer(context$readout_day %||% (release_day + horizon_day))
init_EIR_cal       <- context$init_EIR
theta              <- context$theta
mv_settings        <- context$movement_settings
contact_multiplier <- context$contact_multiplier
seas               <- context$seasonality
RNG_SEED           <- context$seed

n_nodes <- nrow(nodes)
NH_per_node <- as.integer(nodes$NH_per_node)
land <- list(nodes = nodes,
             D = as.matrix(stats::dist(nodes[, c("x", "y")])),
             n_nodes = n_nodes)

allowed_mat <- matrix(TRUE, n_nodes, n_nodes); diag(allowed_mat) <- FALSE
setup <- list(D = land$D, allowed = allowed_mat)

cube <- build_seven_node_drive_cube()
td <- seven_node_trial_design()
# Run a short margin past the readout day so trailing incidence windows and
# plots have room; summaries still report the configured readout day.
tmax_release <- as.integer(readout_day + 35L)

contact_surface <- list(
  type = "contact_surface",
  contact_multiplier = stats::setNames(
    as.numeric(contact_multiplier),
    as.character(seq_len(n_nodes))
  )
)

parameter_modifier <- function(parameters, node_index, warmup_days) {
  parameters$model_seasonality <- TRUE
  parameters$g0 <- seas$g0
  parameters$g  <- seas$g
  parameters$h  <- seas$h
  parameters$rainfall_floor <- seas$rainfall_floor
  malariasimulationGD::apply_node_contact_surface(
    parameters     = parameters,
    contact_surface = contact_surface,
    node_index     = as.integer(node_index)
  )
}

release_spec <- list(
  nodes    = release_nodes,
  time     = release_day,
  size     = td$release_size,
  stage    = "M",
  genotype = "HH"
)

cat(sprintf("[%s] careful/02: resume from checkpoint, tmax = %d days post-warmup\n",
            format(Sys.time(), "%F %T"), tmax_release))
cat(sprintf("  release_nodes = %s | release_day = %d | release_size = %d\n",
            paste(release_nodes, collapse = ","), release_day, td$release_size))

run_one_arm <- function(release_arg, label) {
  cat(sprintf("[%s] running %s arm...\n",
              format(Sys.time(), "%F %T"), label))
  t0 <- Sys.time()
  res <- msimGD_run_truth(
    setup        = setup,
    cube         = cube,
    NF           = NULL,
    NH           = NH_per_node,
    tmax         = tmax_release,
    mu           = mv_settings$mu,
    p_move       = mv_settings$p_move,
    release      = release_arg,
    theta        = theta,
    init_EIR     = init_EIR_cal,
    prevalence_rendering_min_age = 730L,    # 2 years (PfPR2-10 band)
    prevalence_rendering_max_age = 3650L,   # 10 years
    clinical_incidence_min_age   = 182L,    # 6 months
    clinical_incidence_max_age   = 5475L,   # 15 years (production audit-cell band)
    warmup_days  = 0L,
    parameter_modifier = parameter_modifier,
    baseline_checkpoint = checkpoint,
    seed = RNG_SEED
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%s] %s arm finished in %.1f s\n",
              format(Sys.time(), "%F %T"), label, elapsed))
  list(res = res, elapsed = elapsed)
}

out_release    <- run_one_arm(release_spec, "release")
out_no_release <- run_one_arm(NULL,         "no_release")

saveRDS(out_release$res,    file.path(out_dir, "release_truth.rds"))
saveRDS(out_no_release$res, file.path(out_dir, "no_release_truth.rds"))

# --- Extract per-node timeseries to CSV -------------------------------------

stack_data <- function(res, arm_label) {
  # msimGD_run_truth returns a *list of per-node data.frames* directly
  # (each carrying per-day epi + mosquito counts; genotype counts are on
  # attributes named "mosquito_genotype_counts_female"/"_male").
  data_list <- if (is.list(res) && all(vapply(res, is.data.frame, logical(1)))) {
    res
  } else if (!is.null(res$data)) {
    res$data
  } else if (!is.null(res$sim$data)) {
    res$sim$data
  } else {
    stop("Could not locate per-node data.frames in msimGD_run_truth output.",
         call. = FALSE)
  }
  ts_df <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    d <- data_list[[i]]
    d$node <- i
    d$arm  <- arm_label
    d
  }))
  front <- intersect(c("node", "arm", "timestep"), names(ts_df))
  ts_df[, c(front, setdiff(names(ts_df), front))]
}

# Helper to grab the per-node data list regardless of wrapping
get_data_list <- function(res) {
  if (is.list(res) && all(vapply(res, is.data.frame, logical(1)))) res
  else if (!is.null(res$data)) res$data
  else if (!is.null(res$sim$data)) res$sim$data
  else stop("Could not locate per-node data.frames.", call. = FALSE)
}

ts_release    <- stack_data(out_release$res,    "release")
ts_no_release <- stack_data(out_no_release$res, "no_release")
utils::write.csv(ts_release,
                 file.path(out_dir, "release_timeseries.csv"),    row.names = FALSE)
utils::write.csv(ts_no_release,
                 file.path(out_dir, "no_release_timeseries.csv"), row.names = FALSE)

# --- Carrier frequency (release arm) ----------------------------------------

H_genos_definition <- c("WH", "HH", "HR", "HB")
cf_rows <- list()
release_attr_per_node <- vector("list", n_nodes)
first_H_per_node <- integer(n_nodes)
data_release <- get_data_list(out_release$res)
for (i in seq_len(n_nodes)) {
  female <- attr(data_release[[i]], "mosquito_genotype_counts_female")
  male   <- attr(data_release[[i]], "mosquito_genotype_counts_male")
  H_genos <- intersect(H_genos_definition, colnames(female))
  H_F <- rowSums(female[, H_genos, drop = FALSE])
  H_M <- rowSums(male  [, H_genos, drop = FALSE])
  totF <- rowSums(female); totM <- rowSums(male)
  cf <- (H_F + H_M) / pmax(totF + totM, 1)
  cf_rows[[i]] <- data.frame(node = i, time = seq_along(cf),
                             carrier_freq = as.numeric(cf))
  fd <- which((H_F + H_M) > 0)
  first_H_per_node[i] <- if (length(fd) > 0L) fd[[1L]] else NA_integer_
  release_attr_per_node[[i]] <- attr(data_release[[i]],
                                     "mosquito_release_schedule")
}
carrier_df <- do.call(rbind, cf_rows)
utils::write.csv(carrier_df,
                 file.path(out_dir, "carrier_frequency.csv"), row.names = FALSE)

# --- Regional prevalence comparison (PfPR2-10, microscopy) ------------------

# msimGD's prevalence rendering produces a count detected per node per day in
# columns like `n_detect_lm_730_3650` and `n_age_730_3650`. We aggregate the
# numerator and denominator across the 7 villages, per day, per arm.
regional_prevalence <- function(ts_df, arm_label) {
  agg <- aggregate(ts_df[, c("n_detect_lm_730_3650", "n_age_730_3650")],
                   by = list(time = ts_df$timestep), FUN = sum)
  data.frame(
    time = agg$time,
    arm  = arm_label,
    prevalence = agg$n_detect_lm_730_3650 / pmax(agg$n_age_730_3650, 1),
    stringsAsFactors = FALSE
  )
}
prev_df <- rbind(
  regional_prevalence(ts_release,    "release"),
  regional_prevalence(ts_no_release, "no_release")
)
utils::write.csv(prev_df, file.path(out_dir, "prevalence_compare.csv"),
                 row.names = FALSE)

# --- Regional clinical incidence comparison (14-day window counts) ----------

regional_incidence <- function(ts_df, arm_label, window_days = 14L) {
  daily <- aggregate(ts_df$n_infections,
                     by = list(time = ts_df$timestep), FUN = sum)
  names(daily)[2] <- "count"
  daily$window_end_day <-
    ((daily$time - 1L) %/% window_days) * window_days + window_days
  win <- aggregate(daily$count,
                   by = list(window_end_day = daily$window_end_day), FUN = sum)
  names(win)[2] <- "count"
  win$arm <- arm_label
  win
}
inc_df <- rbind(
  regional_incidence(ts_release,    "release"),
  regional_incidence(ts_no_release, "no_release")
)
utils::write.csv(inc_df, file.path(out_dir, "incidence_compare.csv"),
                 row.names = FALSE)

# --- Summary ---------------------------------------------------------------

summary_df <- data.frame(
  workflow                = "careful",
  warmup_days             = context$burnin_timesteps %||% context$warmup_days,
  promoted_snapshot_index = context$promoted_index %||% NA_integer_,
  promoted_snapshot_t_days = context$promoted_t_days %||% NA_integer_,
  multi_seed_pass         = context$multi_seed_pass %||% NA,
  init_EIR_calibrated     = init_EIR_cal,
  tmax_release            = tmax_release,
  release_day             = release_day,
  horizon_day             = horizon_day,
  readout_day             = readout_day,
  n_nodes                 = n_nodes,
  release_nodes           = paste(release_nodes, collapse = ","),
  release_size            = td$release_size,
  first_H_carrier_day     = min(first_H_per_node, na.rm = TRUE),
  peak_carrier_release    = max(carrier_df$carrier_freq[
                                  carrier_df$node %in% release_nodes],
                                na.rm = TRUE),
  peak_carrier_nonrelease = max(carrier_df$carrier_freq[
                                  !(carrier_df$node %in% release_nodes)],
                                na.rm = TRUE),
  release_arm_seconds     = round(out_release$elapsed,    1),
  no_release_arm_seconds  = round(out_no_release$elapsed, 1),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_df, file.path(out_dir, "summary.csv"), row.names = FALSE)

print(summary_df)
cat(sprintf("Outputs saved under: %s\n", out_dir))
