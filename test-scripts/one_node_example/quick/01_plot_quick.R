#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 01_plot_quick.R
# ------------------------------------------------------------------------------
# Plot the CSVs written by 00_run_quick.R.
#
# Run from the package root:
#   Rscript test-scripts/one_node_example/quick/01_plot_quick.R
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/one_node_example/quick/01_plot_quick.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out_dir <- file.path(example_dir, "output", "quick")
stopifnot(dir.exists(out_dir))

timeseries <- utils::read.csv(file.path(out_dir, "timeseries.csv"), stringsAsFactors = FALSE)
carrier <- utils::read.csv(file.path(out_dir, "carrier_frequency.csv"), stringsAsFactors = FALSE)
context <- readRDS(file.path(out_dir, "context.rds"))

if (!("pfpr_2_10_lm" %in% names(timeseries))) {
  timeseries$pfpr_2_10_lm <-
    timeseries$n_detect_lm_730_3650 / pmax(timeseries$n_age_730_3650, 1)
}

release_day <- as.integer(context$release_day)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  release_line <- geom_vline(
    xintercept = release_day,
    linetype = "22",
    linewidth = 0.45,
    color = "#555555"
  )
  base_theme <- theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot"
    )

  p_carrier <- ggplot(carrier, aes(timestep, adult_drive_carrier_frequency)) +
    geom_line(color = "#8f2d1f", linewidth = 0.8) +
    release_line +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = "Adult Drive Carrier Frequency",
      x = "Day",
      y = "Adult carriers"
    ) +
    base_theme

  p_incidence <- ggplot(timeseries, aes(timestep, n_infections)) +
    geom_line(color = "#1f6f8b", linewidth = 0.7) +
    release_line +
    labs(
      title = "Malaria Incidence",
      x = "Day",
      y = "New infections"
    ) +
    base_theme

  p_prevalence <- ggplot(timeseries, aes(timestep, pfpr_2_10_lm)) +
    geom_line(color = "#3d7d2a", linewidth = 0.8) +
    release_line +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    labs(
      title = "Child Prevalence",
      x = "Day",
      y = "PfPR 2-10, microscopy"
    ) +
    base_theme

  ggplot2::ggsave(
    file.path(out_dir, "adult_drive_carrier_frequency.png"),
    p_carrier,
    width = 8,
    height = 4.5,
    units = "in",
    dpi = 150
  )
  ggplot2::ggsave(
    file.path(out_dir, "incidence_over_time.png"),
    p_incidence,
    width = 8,
    height = 4.5,
    units = "in",
    dpi = 150
  )
  ggplot2::ggsave(
    file.path(out_dir, "child_prevalence_over_time.png"),
    p_prevalence,
    width = 8,
    height = 4.5,
    units = "in",
    dpi = 150
  )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- (p_carrier / p_incidence / p_prevalence) +
      patchwork::plot_annotation(
        title = "One-Node Native Homing-Drive Quick Summary",
        subtitle = sprintf(
          "Release: %s %s mosquitoes on day %d; genotype fitness costs in cube$s",
          context$release_count,
          context$release_genotype,
          release_day
        )
      )
    ggplot2::ggsave(
      file.path(out_dir, "one_node_quick_summary.png"),
      combined,
      width = 9,
      height = 10,
      units = "in",
      dpi = 150
    )
  }
} else {
  message("Package 'ggplot2' is not available; using base R plots.")
  png(file.path(out_dir, "adult_drive_carrier_frequency.png"), width = 1200, height = 700)
  plot(
    carrier$timestep,
    carrier$adult_drive_carrier_frequency,
    type = "l",
    xlab = "Day",
    ylab = "Adult drive carrier frequency",
    main = "Adult Drive Carrier Frequency"
  )
  abline(v = release_day, lty = 2)
  dev.off()

  png(file.path(out_dir, "incidence_over_time.png"), width = 1200, height = 700)
  plot(
    timeseries$timestep,
    timeseries$n_infections,
    type = "l",
    xlab = "Day",
    ylab = "New infections",
    main = "Malaria Incidence"
  )
  abline(v = release_day, lty = 2)
  dev.off()

  png(file.path(out_dir, "child_prevalence_over_time.png"), width = 1200, height = 700)
  plot(
    timeseries$timestep,
    timeseries$pfpr_2_10_lm,
    type = "l",
    xlab = "Day",
    ylab = "PfPR 2-10, microscopy",
    main = "Child Prevalence"
  )
  abline(v = release_day, lty = 2)
  dev.off()
}

cat("Saved figures:\n")
cat(sprintf("  %s\n", file.path(out_dir, "adult_drive_carrier_frequency.png")))
cat(sprintf("  %s\n", file.path(out_dir, "incidence_over_time.png")))
cat(sprintf("  %s\n", file.path(out_dir, "child_prevalence_over_time.png")))
if (file.exists(file.path(out_dir, "one_node_quick_summary.png"))) {
  cat(sprintf("  %s\n", file.path(out_dir, "one_node_quick_summary.png")))
}
