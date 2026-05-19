#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 01_plot_quick.R
# ------------------------------------------------------------------------------
# Reads the CSVs and context written by 00_run_quick.R and renders the
# four-panel diagnostic figure:
#   A) landscape map (with release villages flagged)
#   B) carrier-frequency heatmap (village x time)
#   C) regional infection incidence (sum across villages, per 14-day window)
#   D) regional PfPR2-10 microscopy prevalence (default msimGD render band)
#
# Outputs (under output/quick/):
#   figure_quick.png
#   figure_quick.pdf
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/seven_node_example/quick/01_plot_quick.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(example_dir, "lib", "plotting.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

out_dir <- file.path(example_dir, "output", "quick")
stopifnot(dir.exists(out_dir))

nodes      <- utils::read.csv(file.path(out_dir, "nodes.csv"),             stringsAsFactors = FALSE)
ts         <- utils::read.csv(file.path(out_dir, "timeseries.csv"),        stringsAsFactors = FALSE)
carrier_df <- utils::read.csv(file.path(out_dir, "carrier_frequency.csv"), stringsAsFactors = FALSE)
context    <- readRDS(file.path(out_dir, "context.rds"))

release_nodes <- as.integer(context$release_nodes)
release_day   <- as.integer(context$release_day)
horizon_day   <- as.integer(context$horizon_day)
readout_day   <- as.integer(context$readout_day %||% (release_day + horizon_day))

n_nodes <- max(nodes$node)
tmax    <- max(ts$timestep)

# ---- Panel A: landscape ----------------------------------------------------
panel_A <- plot_landscape(nodes, release_nodes)

# ---- Panel B: carrier-frequency heatmap -----------------------------------
panel_B <- plot_carrier_heatmap(
  carrier_df,
  release_nodes = release_nodes,
  release_day   = release_day,
  readout_day   = readout_day
)

# ---- Panel C: regional infection incidence (14-day windows) ----------------
# Aggregate per-day n_infections across the 7 villages, then sum into
# non-overlapping 14-day windows so the panel reads cleanly.
window_days <- 14L
ts_regional_daily <- aggregate(
  ts$n_infections,
  by = list(timestep = ts$timestep),
  FUN = sum
)
names(ts_regional_daily)[2] <- "n_infections"
ts_regional_daily$window_end_day <-
  ((ts_regional_daily$timestep - 1L) %/% window_days) * window_days + window_days
inc_windowed <- aggregate(
  ts_regional_daily$n_infections,
  by = list(window_end_day = ts_regional_daily$window_end_day),
  FUN = sum
)
names(inc_windowed)[2] <- "count"
inc_windowed$arm <- "release"
panel_C <- plot_regional_incidence(
  inc_windowed,
  release_day = release_day,
  readout_day = readout_day
)

# ---- Panel D: regional PfPR2-10 microscopy prevalence ---------------------
# msimGD's default render gives n_detect_lm_730_3650 (count microscopy+ in
# 2-10y) and n_age_730_3650 (total in band). Regional prevalence is sum
# numerator / sum denominator across villages, per day.
prev_regional <- aggregate(
  ts[, c("n_detect_lm_730_3650", "n_age_730_3650")],
  by = list(timestep = ts$timestep),
  FUN = sum
)
prev_regional$prevalence <-
  prev_regional$n_detect_lm_730_3650 / pmax(prev_regional$n_age_730_3650, 1)
prev_df <- data.frame(
  time = prev_regional$timestep,
  arm  = "release",
  prevalence = prev_regional$prevalence
)
panel_D <- plot_regional_prevalence(
  prev_df,
  calibration_target = NULL,    # quick path doesn't calibrate to a target
  release_day = release_day,
  readout_day = readout_day
)

# ---- Compose + save --------------------------------------------------------
fig <- compose_four_panels(panel_A, panel_B, panel_C, panel_D)

png_path <- file.path(out_dir, "figure_quick.png")
pdf_path <- file.path(out_dir, "figure_quick.pdf")
ggplot2::ggsave(png_path, fig, width = 12, height = 9, units = "in", dpi = 150)
ggplot2::ggsave(pdf_path, fig, width = 12, height = 9, units = "in")

cat(sprintf("Saved: %s\n", png_path))
cat(sprintf("Saved: %s\n", pdf_path))

# Brief sanity printout
prev_idx <- which.min(abs(prev_regional$timestep - readout_day))
inc_idx <- which.min(abs(inc_windowed$window_end_day - readout_day))
cat(sprintf("\nReadout-day (%d) regional prevalence (2-10y, lm): %.3f\n",
            readout_day, prev_regional$prevalence[prev_idx]))
cat(sprintf("Readout-day (%d) regional 14-day infection count: %d\n",
            readout_day, inc_windowed$count[inc_idx]))
