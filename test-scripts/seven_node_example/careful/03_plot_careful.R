#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 03_plot_careful.R
# ------------------------------------------------------------------------------
# Renders the four-panel diagnostic figure for the careful workflow:
#   A) landscape with release villages flagged
#   B) carrier-frequency heatmap (release arm)
#   C) regional clinical incidence (release vs no-release counterfactual)
#   D) regional PfPR2-10 microscopy prevalence (release vs no-release)
#
# Input artifacts (under output/careful/):
#   context.rds           (release_nodes, release_day, etc.)
#   carrier_frequency.csv (release-arm drive-allele carrier frequency)
#   prevalence_compare.csv
#   incidence_compare.csv
#   release_timeseries.csv
#   nodes.csv (NOT written by current 02_; we rebuild from context$nodes here)
#
# Output:
#   figure_careful.png / .pdf
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/seven_node_example/careful/03_plot_careful.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(example_dir, "lib", "plotting.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

out_dir <- file.path(example_dir, "output", "careful")
stopifnot(dir.exists(out_dir))

context    <- readRDS(file.path(out_dir, "context.rds"))
carrier_df <- utils::read.csv(file.path(out_dir, "carrier_frequency.csv"),
                              stringsAsFactors = FALSE)
prev_df <- utils::read.csv(file.path(out_dir, "prevalence_compare.csv"),
                           stringsAsFactors = FALSE)
inc_df  <- utils::read.csv(file.path(out_dir, "incidence_compare.csv"),
                           stringsAsFactors = FALSE)

nodes         <- context$nodes
release_nodes <- as.integer(context$release_nodes)
release_day   <- as.integer(context$release_day)
horizon_day   <- as.integer(context$horizon_day)
readout_day   <- as.integer(context$readout_day %||% (release_day + horizon_day))

# Pull the calibration target so we can draw it on the prevalence panel.
cal <- readRDS(file.path(out_dir, "calibrated_init_eir.rds"))
target_prev <- cal$target_prevalence

# ---- Panel A: landscape ---------------------------------------------------
panel_A <- plot_landscape(nodes, release_nodes)

# ---- Panel B: carrier-frequency heatmap (release arm) --------------------
panel_B <- plot_carrier_heatmap(
  carrier_df, release_nodes = release_nodes,
  release_day = release_day, readout_day = readout_day
)

# ---- Panel C: regional incidence, both arms ------------------------------
panel_C <- plot_regional_incidence(
  inc_df, release_day = release_day, readout_day = readout_day
)

# ---- Panel D: regional prevalence, both arms -----------------------------
panel_D <- plot_regional_prevalence(
  prev_df, calibration_target = target_prev,
  release_day = release_day, readout_day = readout_day
)

fig <- compose_four_panels(panel_A, panel_B, panel_C, panel_D)

png_path <- file.path(out_dir, "figure_careful.png")
pdf_path <- file.path(out_dir, "figure_careful.pdf")
ggplot2::ggsave(png_path, fig, width = 12, height = 9, units = "in", dpi = 150)
ggplot2::ggsave(pdf_path, fig, width = 12, height = 9, units = "in")
cat(sprintf("Saved: %s\n", png_path))
cat(sprintf("Saved: %s\n", pdf_path))

# Brief sanity
prev_at_readout <- function(arm) {
  idx <- which(prev_df$arm == arm)
  if (length(idx) == 0L) return(NA_real_)
  idx <- idx[which.min(abs(prev_df$time[idx] - readout_day))]
  prev_df$prevalence[idx]
}
cat(sprintf("\nReadout-day (%d) PfPR2-10: release = %.3f | no_release = %.3f\n",
            readout_day, prev_at_readout("release"),
            prev_at_readout("no_release")))
