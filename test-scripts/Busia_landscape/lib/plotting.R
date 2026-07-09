# ------------------------------------------------------------------------------
# plotting.R
# ------------------------------------------------------------------------------
# ggplot helpers for the four diagnostic panels:
#   A) landscape map with release villages flagged
#   B) carrier-frequency heatmap (node x time)
#   C) regional clinical incidence (count per 2-week window)
#   D) regional microscopy prevalence in the calibration age band
#
# All four panels return a ggplot object. They are composed by the workflow
# scripts via patchwork.
# ------------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required for plotting.", call. = FALSE)
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package 'scales' is required for plotting.", call. = FALSE)
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required for plotting.", call. = FALSE)
}

theme_seven_node <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#e6ebef", linewidth = 0.3),
      panel.border = ggplot2::element_rect(color = "#c4cdd2", fill = NA, linewidth = 0.4),
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )
}


#' Panel A: landscape map
plot_landscape <- function(nodes_df, release_nodes) {
  nodes_df$role <- ifelse(nodes_df$node %in% release_nodes,
                          "Release", "Non-release")
  nodes_df$role <- factor(nodes_df$role, levels = c("Non-release", "Release"))
  fill_pal  <- c("Non-release" = "#5b6770", "Release" = "#b24d1b")
  color_pal <- c("Non-release" = "#2f3a40", "Release" = "#7a3712")
  ggplot2::ggplot(nodes_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(fill = role, color = role),
                        shape = 21, size = 4.2, stroke = 0.6, alpha = 0.95) +
    ggplot2::geom_text(ggplot2::aes(label = node, color = role),
                       fontface = "bold", size = 3.2,
                       nudge_y = 0.8, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = fill_pal, name = NULL) +
    ggplot2::scale_color_manual(values = color_pal, name = NULL) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "Landscape",
                  x = "East-west position (km)",
                  y = "North-south position (km)") +
    theme_seven_node()
}


#' Panel B: carrier-frequency heatmap (node x time)
plot_carrier_heatmap <- function(carrier_df, release_nodes,
                                 release_day = NULL, readout_day = NULL) {
  carrier_df$node <- as.integer(carrier_df$node)
  carrier_df$time <- as.integer(carrier_df$time)
  n_nodes <- max(carrier_df$node)
  node_labels <- ifelse(seq_len(n_nodes) %in% release_nodes,
                        sprintf("%d *", seq_len(n_nodes)),
                        as.character(seq_len(n_nodes)))
  p <- ggplot2::ggplot(carrier_df, ggplot2::aes(x = time, y = node, fill = carrier_freq)) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::scale_fill_viridis_c(
      name = "Carrier\nfrequency",
      limits = c(0, 1), option = "magma", direction = 1,
      labels = scales::label_percent(accuracy = 1)
    ) +
    ggplot2::scale_y_continuous(breaks = seq_len(n_nodes), labels = node_labels,
                                expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = "Drive-allele carrier frequency",
                  x = "Time (days)", y = "Village")
  if (!is.null(release_day)) {
    p <- p + ggplot2::geom_vline(xintercept = release_day, linetype = "22",
                                 color = "white", linewidth = 0.45, alpha = 0.9)
  }
  if (!is.null(readout_day)) {
    p <- p + ggplot2::geom_vline(xintercept = readout_day,
                                 color = "white", linewidth = 0.3, alpha = 0.45)
  }
  p + theme_seven_node()
}


#' Panel C: regional clinical incidence (counts in a 14-day window)
plot_regional_incidence <- function(incidence_df,
                                    release_day = NULL,
                                    readout_day = NULL) {
  # incidence_df: long, columns (window_end_day, arm, count), where
  # arm in {"release", "no_release"} (no_release column may be absent for quick path)
  arm_levels <- intersect(c("no_release", "release"), unique(incidence_df$arm))
  incidence_df$arm <- factor(incidence_df$arm, levels = arm_levels)
  arm_colors <- c("no_release" = "#1f4f9c", "release" = "#b24d1b")
  p <- ggplot2::ggplot(incidence_df,
                       ggplot2::aes(x = window_end_day, y = count, color = arm)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_color_manual(values = arm_colors[arm_levels],
                                labels = c(no_release = "No release",
                                           release    = "Release")[arm_levels],
                                name = NULL) +
    ggplot2::labs(title = "Regional infection incidence (count / 14 days)",
                  x = "Time (days)", y = "New infections")
  if (!is.null(release_day)) {
    p <- p + ggplot2::geom_vline(xintercept = release_day, linetype = "22",
                                 color = "#2d2d2d", linewidth = 0.45)
  }
  if (!is.null(readout_day)) {
    p <- p + ggplot2::geom_vline(xintercept = readout_day,
                                 color = "#c4c4c4", linewidth = 0.3)
  }
  p + theme_seven_node()
}


#' Panel D: regional microscopy prevalence (proportion)
plot_regional_prevalence <- function(prevalence_df,
                                     calibration_target = NULL,
                                     release_day = NULL,
                                     readout_day = NULL) {
  # prevalence_df: long, columns (time, arm, prevalence)
  arm_levels <- intersect(c("no_release", "release"), unique(prevalence_df$arm))
  prevalence_df$arm <- factor(prevalence_df$arm, levels = arm_levels)
  arm_colors <- c("no_release" = "#1f4f9c", "release" = "#b24d1b")
  p <- ggplot2::ggplot(prevalence_df,
                       ggplot2::aes(x = time, y = prevalence, color = arm)) +
    ggplot2::geom_line(linewidth = 0.9, alpha = 0.95) +
    ggplot2::scale_color_manual(values = arm_colors[arm_levels],
                                labels = c(no_release = "No release",
                                           release    = "Release")[arm_levels],
                                name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(title = "Regional microscopy prevalence (PfPR 2-10y)",
                  x = "Time (days)", y = "Prevalence")
  if (!is.null(calibration_target)) {
    p <- p + ggplot2::geom_hline(yintercept = calibration_target,
                                 linetype = "dotted", color = "#6f6f6f",
                                 linewidth = 0.45)
  }
  if (!is.null(release_day)) {
    p <- p + ggplot2::geom_vline(xintercept = release_day, linetype = "22",
                                 color = "#2d2d2d", linewidth = 0.45)
  }
  if (!is.null(readout_day)) {
    p <- p + ggplot2::geom_vline(xintercept = readout_day,
                                 color = "#c4c4c4", linewidth = 0.3)
  }
  p + theme_seven_node()
}


#' Compose all four panels (2x2 grid) using patchwork.
compose_four_panels <- function(panel_A, panel_B, panel_C, panel_D) {
  ( (panel_A | panel_B) / (panel_C | panel_D) ) +
    patchwork::plot_annotation(tag_levels = "A")
}
