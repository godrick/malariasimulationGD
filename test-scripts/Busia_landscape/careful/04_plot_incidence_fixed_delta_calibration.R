#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 04_plot_incidence_fixed_delta_calibration.R
# ------------------------------------------------------------------------------
# Publication-quality figures for the fixed-delta Busia incidence calibration.
#
# Input directory:
#   BUSIA_CAREFUL_OUT_DIR, default: test-scripts/Busia_landscape/output/careful
#
# Expected inputs from 00_calibrate_init_eir_incidence_fixed_delta.R:
#   incidence_calibration_fixed_delta_log.csv
#   incidence_calibration_fixed_delta_global_summary.csv
#   incidence_calibration_fixed_delta_node_summary.csv
#
# Optional inputs, when BUSIA_RUN_COUNTERFACTUAL_NO_INTERVENTION=true was used:
#   incidence_calibration_fixed_delta_no_intervention_global_summary.csv
#   incidence_calibration_fixed_delta_no_intervention_node_summary.csv
#
# Outputs:
#   BUSIA_PLOT_OUT_DIR, default:
#     <BUSIA_CAREFUL_OUT_DIR>/figures/incidence_fixed_delta
# ------------------------------------------------------------------------------ 

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/Busia_landscape/careful/04_plot_incidence_fixed_delta_calibration.R",
    mustWork = TRUE
  )
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to render calibration figures.",
       call. = FALSE)
}

source(file.path(example_dir, "config", "busia_landscape.R"))

out_dir <- Sys.getenv(
  "BUSIA_CAREFUL_OUT_DIR",
  unset = file.path(example_dir, "output", "careful")
)
plot_dir <- Sys.getenv(
  "BUSIA_PLOT_OUT_DIR",
  unset = file.path(out_dir, "figures", "incidence_fixed_delta")
)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

read_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

required_col <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) < 1L) {
    stop(sprintf(
      "Could not find %s column. Tried: %s",
      label,
      paste(candidates, collapse = ", ")
    ), call. = FALSE)
  }
  hit[[1L]]
}

optional_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) < 1L) NULL else hit[[1L]]
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_min <- function(x) {
  x <- safe_num(x)
  if (length(x) < 1L || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- safe_num(x)
  if (length(x) < 1L || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  x <- safe_num(x)
  if (length(x) < 1L || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

theme_publication <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "#3F454A"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "#24282C"),
      strip.background = ggplot2::element_rect(fill = "#EDF0F2", color = "#AAB2B8"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "#E4E8EB", linewidth = 0.25),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

save_plot <- function(plot, filename, width, height, dpi = 320) {
  png_path <- file.path(plot_dir, paste0(filename, ".png"))
  pdf_path <- file.path(plot_dir, paste0(filename, ".pdf"))
  ggplot2::ggsave(png_path, plot, width = width, height = height,
                  units = "in", dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                  units = "in", bg = "white")
  c(png = png_path, pdf = pdf_path)
}

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

cal_log <- read_required_csv(
  file.path(out_dir, "incidence_calibration_fixed_delta_log.csv"),
  "fixed-delta calibration log"
)
global <- read_required_csv(
  file.path(out_dir, "incidence_calibration_fixed_delta_global_summary.csv"),
  "fixed-delta global summary"
)
nodes <- read_required_csv(
  file.path(out_dir, "incidence_calibration_fixed_delta_node_summary.csv"),
  "fixed-delta node summary"
)
cf_global <- read_optional_csv(
  file.path(out_dir, "incidence_calibration_fixed_delta_no_intervention_global_summary.csv")
)
cf_nodes <- read_optional_csv(
  file.path(out_dir, "incidence_calibration_fixed_delta_no_intervention_node_summary.csv")
)

if (nrow(global) < 1L) {
  stop("Global summary is empty.", call. = FALSE)
}
if (nrow(nodes) != 58L) {
  warning(sprintf("Node summary has %d rows; expected 58.", nrow(nodes)),
          call. = FALSE)
}

eir_col <- required_col(cal_log, c("raw_counterfactual_init_EIR", "init_EIR"),
                        "candidate raw/init EIR")
mean_col <- required_col(cal_log,
                         c("unweighted_mean_incidence_per_person_year"),
                         "candidate mean incidence")
sd_col <- required_col(cal_log, c("site_sd_incidence_per_person_year"),
                       "candidate site-level SD")
cv_col <- required_col(cal_log, c("between_cluster_hazard_cv"),
                       "candidate CV")
rel_col <- optional_col(cal_log, c("mean_relative_error"))
kind_col <- optional_col(cal_log, c("kind"))
run_col <- optional_col(cal_log, c("run_seconds"))
candidate_realized_eir_col <- optional_col(
  cal_log,
  c("mean_realized_eir", "mean_realised_eir", "realized_eir", "realised_eir")
)

node_inc_col <- required_col(nodes, c("infection_incidence_per_person_year"),
                             "node incidence")
node_id_col <- required_col(nodes, c("node_id", "node"), "node id")
cluster_col <- required_col(nodes, c("cluster_id"), "cluster id")
contact_col <- optional_col(nodes, c("contact_multiplier"))
node_eir_col <- optional_col(nodes, c("mean_eir", "mean_realized_eir"))
pop_col <- optional_col(nodes, c("human_population", "NH_per_node"))

target_mean <- safe_num(global$target_mean_incidence_per_person_year[[1L]])
target_sd <- safe_num(global$target_site_sd_incidence_per_person_year[[1L]])
target_cv <- safe_num(global$target_between_cluster_hazard_cv[[1L]])
range_min <- safe_num(global$reference_range_min_incidence_per_person_year[[1L]])
range_max <- safe_num(global$reference_range_max_incidence_per_person_year[[1L]])
selected_eir <- safe_num((global$raw_counterfactual_init_EIR %||% global$init_EIR)[[1L]])
selected_delta <- safe_num((global$fixed_delta %||% global$covariate_delta)[[1L]])
selected_mean <- safe_num(global$unweighted_mean_incidence_per_person_year[[1L]])
selected_sd <- safe_num(global$site_sd_incidence_per_person_year[[1L]])
selected_cv <- safe_num(global$between_cluster_hazard_cv[[1L]])
selected_realized_eir <- if ("mean_realized_eir" %in% names(global)) {
  safe_num(global$mean_realized_eir[[1L]])
} else if (!is.null(node_eir_col)) {
  safe_mean(nodes[[node_eir_col]])
} else {
  NA_real_
}
global$raw_init_EIR_plot <- selected_eir
global$selected_mean_incidence_plot <- selected_mean

cal_log$raw_init_EIR_plot <- safe_num(cal_log[[eir_col]])
cal_log$mean_incidence_plot <- safe_num(cal_log[[mean_col]])
cal_log$sd_incidence_plot <- safe_num(cal_log[[sd_col]])
cal_log$cv_plot <- safe_num(cal_log[[cv_col]])
cal_log$relative_error_plot <- if (!is.null(rel_col)) {
  safe_num(cal_log[[rel_col]])
} else {
  abs(cal_log$mean_incidence_plot - target_mean) / target_mean
}
cal_log$kind_plot <- if (!is.null(kind_col)) as.character(cal_log[[kind_col]]) else "candidate"
cal_log$run_seconds_plot <- if (!is.null(run_col)) safe_num(cal_log[[run_col]]) else NA_real_
cal_log$realized_eir_plot <- if (!is.null(candidate_realized_eir_col)) {
  safe_num(cal_log[[candidate_realized_eir_col]])
} else {
  NA_real_
}
cal_log$candidate_index <- seq_len(nrow(cal_log))

nodes$node_id_plot <- as.integer(nodes[[node_id_col]])
nodes$cluster_id_plot <- as.integer(nodes[[cluster_col]])
nodes$incidence_plot <- safe_num(nodes[[node_inc_col]])
nodes$contact_plot <- if (!is.null(contact_col)) safe_num(nodes[[contact_col]]) else NA_real_
nodes$mean_eir_plot <- if (!is.null(node_eir_col)) safe_num(nodes[[node_eir_col]]) else NA_real_
nodes$population_plot <- if (!is.null(pop_col)) safe_num(nodes[[pop_col]]) else NA_real_
nodes$rank_incidence <- rank(nodes$incidence_plot, ties.method = "first")
nodes$status_plot <- ifelse(
  nodes$incidence_plot < range_min,
  "Below observed range",
  ifelse(nodes$incidence_plot > range_max, "Above observed range", "Within observed range")
)
nodes$status_plot <- factor(
  nodes$status_plot,
  levels = c("Below observed range", "Within observed range", "Above observed range")
)

land <- build_busia_landscape()
nodes <- merge(
  nodes,
  land$nodes[, c("node", "latitude", "longitude", "x", "y")],
  by.x = "node_id_plot",
  by.y = "node",
  all.x = TRUE,
  sort = FALSE
)

summary_table <- data.frame(
  metric = c(
    "raw_counterfactual_init_EIR",
    "fixed_delta",
    "mean_incidence",
    "target_mean_incidence",
    "site_sd_incidence",
    "target_site_sd_incidence",
    "between_cluster_cv",
    "target_between_cluster_cv",
    "mean_realized_eir",
    "min_site_incidence",
    "max_site_incidence",
    "sites_inside_observed_range"
  ),
  value = c(
    selected_eir,
    selected_delta,
    selected_mean,
    target_mean,
    selected_sd,
    target_sd,
    selected_cv,
    target_cv,
    selected_realized_eir,
    safe_min(nodes$incidence_plot),
    safe_max(nodes$incidence_plot),
    sum(nodes$status_plot == "Within observed range", na.rm = TRUE)
  )
)
utils::write.csv(
  summary_table,
  file.path(plot_dir, "fixed_delta_calibration_figure_summary.csv"),
  row.names = FALSE
)

subtitle_text <- sprintf(
  "Fixed delta = %s; calibrated raw EIR = %s; achieved mean = %s infections/person-year",
  format_num(selected_delta, 3),
  format_num(selected_eir, 3),
  format_num(selected_mean, 3)
)
search_x_min <- safe_min(c(cal_log$raw_init_EIR_plot, selected_eir))
search_x_max <- safe_max(c(cal_log$raw_init_EIR_plot, selected_eir))
search_x_min <- ifelse(is.na(search_x_min) || search_x_min <= 0, 0.01, search_x_min)
search_x_max <- ifelse(is.na(search_x_max) || search_x_max <= search_x_min,
                       search_x_min * 10, search_x_max)

p_search <- ggplot2::ggplot(
  cal_log,
  ggplot2::aes(x = raw_init_EIR_plot, y = mean_incidence_plot)
) +
  ggplot2::annotate(
    "rect",
    xmin = search_x_min / 1.2, xmax = search_x_max * 1.2,
    ymin = range_min, ymax = range_max,
    fill = "#D8E7F5", alpha = 0.45
  ) +
  ggplot2::geom_hline(yintercept = target_mean, linewidth = 0.65,
                      linetype = "22", color = "#1F5A8A") +
  ggplot2::geom_line(ggplot2::aes(group = 1L), linewidth = 0.55,
                     color = "#747B82") +
  ggplot2::geom_point(ggplot2::aes(fill = relative_error_plot),
                      shape = 21, size = 3.2, color = "white", stroke = 0.35) +
  ggplot2::geom_point(
    data = global,
    ggplot2::aes(
      x = raw_init_EIR_plot,
      y = selected_mean_incidence_plot
    ),
    inherit.aes = FALSE,
    shape = 23,
    fill = "#C75146",
    color = "black",
    size = 4.1,
    stroke = 0.4
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_fill_gradient(
    low = "#2D7C4F", high = "#C75146",
    name = "Mean relative\nerror"
  ) +
  ggplot2::labs(
    title = "Calibration Search Trajectory",
    subtitle = "Shaded band is the reported observed site range; dashed line is the target mean",
    x = "Raw counterfactual init_EIR (log scale)",
    y = "Mean modeled infection incidence\n(infections/person-year)"
  ) +
  theme_publication()

realized_search_df <- cal_log[is.finite(cal_log$realized_eir_plot) &
                                cal_log$realized_eir_plot > 0, , drop = FALSE]
realized_has_candidates <- nrow(realized_search_df) > 0L
selected_realized_df <- data.frame(
  realized_eir_plot = selected_realized_eir,
  mean_incidence_plot = selected_mean,
  relative_error_plot = abs(selected_mean - target_mean) / target_mean
)
realized_x_min <- safe_min(c(realized_search_df$realized_eir_plot,
                             selected_realized_eir))
realized_x_max <- safe_max(c(realized_search_df$realized_eir_plot,
                             selected_realized_eir))
realized_x_min <- ifelse(is.na(realized_x_min) || realized_x_min <= 0,
                         0.01, realized_x_min)
realized_x_max <- ifelse(is.na(realized_x_max) ||
                           realized_x_max <= realized_x_min,
                         realized_x_min * 10, realized_x_max)

p_search_realized <- ggplot2::ggplot(
  realized_search_df,
  ggplot2::aes(x = realized_eir_plot, y = mean_incidence_plot)
) +
  ggplot2::annotate(
    "rect",
    xmin = realized_x_min / 1.2, xmax = realized_x_max * 1.2,
    ymin = range_min, ymax = range_max,
    fill = "#D8E7F5", alpha = 0.45
  ) +
  ggplot2::geom_hline(yintercept = target_mean, linewidth = 0.65,
                      linetype = "22", color = "#1F5A8A") +
  ggplot2::geom_line(ggplot2::aes(group = 1L), linewidth = 0.55,
                     color = "#747B82") +
  ggplot2::geom_point(ggplot2::aes(fill = relative_error_plot),
                      shape = 21, size = 3.2, color = "white", stroke = 0.35) +
  ggplot2::geom_point(
    data = selected_realized_df,
    ggplot2::aes(x = realized_eir_plot, y = mean_incidence_plot),
    inherit.aes = FALSE,
    shape = 23,
    fill = "#C75146",
    color = "black",
    size = 4.1,
    stroke = 0.4
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_fill_gradient(
    low = "#2D7C4F", high = "#C75146",
    name = "Mean relative\nerror"
  ) +
  ggplot2::labs(
    title = "Calibration Versus Realized Annual EIR",
    subtitle = if (realized_has_candidates) {
      "Candidate-level realized EIR from the calibration log; diamond marks selected result"
    } else {
      "Candidate-level realized EIR was not in the log; diamond uses selected node summary"
    },
    x = "Mean realized annual EIR across sites (log scale)",
    y = "Mean modeled infection incidence\n(infections/person-year)"
  ) +
  theme_publication()

metric_df <- rbind(
  data.frame(
    candidate_index = cal_log$candidate_index,
    kind = cal_log$kind_plot,
    metric = "Mean incidence",
    value = cal_log$mean_incidence_plot,
    target = target_mean
  ),
  data.frame(
    candidate_index = cal_log$candidate_index,
    kind = cal_log$kind_plot,
    metric = "Site SD",
    value = cal_log$sd_incidence_plot,
    target = target_sd
  ),
  data.frame(
    candidate_index = cal_log$candidate_index,
    kind = cal_log$kind_plot,
    metric = "Between-cluster CV",
    value = cal_log$cv_plot,
    target = target_cv
  )
)
metric_df$metric <- factor(metric_df$metric,
                           levels = c("Mean incidence", "Site SD", "Between-cluster CV"))

p_metrics <- ggplot2::ggplot(
  metric_df,
  ggplot2::aes(x = candidate_index, y = value)
) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = target),
                      linetype = "22", color = "#1F5A8A", linewidth = 0.55) +
  ggplot2::geom_line(color = "#5C636A", linewidth = 0.55) +
  ggplot2::geom_point(size = 2.1, color = "#2F6F9F") +
  ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  ggplot2::scale_x_continuous(breaks = cal_log$candidate_index) +
  ggplot2::labs(
    title = "Candidate Diagnostics",
    subtitle = "Dashed line marks each target metric",
    x = "Candidate evaluation order",
    y = NULL
  ) +
  theme_publication()

p_rank <- ggplot2::ggplot(
  nodes[order(nodes$incidence_plot), ],
  ggplot2::aes(x = seq_along(incidence_plot), y = incidence_plot)
) +
  ggplot2::annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = range_min, ymax = range_max,
    fill = "#D8E7F5", alpha = 0.45
  ) +
  ggplot2::geom_hline(yintercept = target_mean, linetype = "22",
                      linewidth = 0.6, color = "#1F5A8A") +
  ggplot2::geom_point(ggplot2::aes(fill = status_plot),
                      shape = 21, size = 2.8, color = "white", stroke = 0.25) +
  ggplot2::scale_fill_manual(
    values = c(
      "Below observed range" = "#7A8A99",
      "Within observed range" = "#2D7C4F",
      "Above observed range" = "#C75146"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Site-Level Incidence Distribution",
    subtitle = subtitle_text,
    x = "Sites ordered by modeled incidence",
    y = "Modeled infection incidence\n(infections/person-year)"
  ) +
  theme_publication()

p_contact <- ggplot2::ggplot(
  nodes,
  ggplot2::aes(x = contact_plot, y = incidence_plot)
) +
  ggplot2::annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = range_min, ymax = range_max,
    fill = "#D8E7F5", alpha = 0.35
  ) +
  ggplot2::geom_hline(yintercept = target_mean, linetype = "22",
                      linewidth = 0.55, color = "#1F5A8A") +
  ggplot2::geom_point(
    ggplot2::aes(size = population_plot, fill = mean_eir_plot),
    shape = 21, color = "white", stroke = 0.25, alpha = 0.92
  ) +
  ggplot2::scale_size_continuous(range = c(1.8, 5), name = "Population") +
  ggplot2::scale_fill_gradient(
    low = "#C6DBEF", high = "#7F1D1D",
    na.value = "#87919A",
    name = "Mean EIR"
  ) +
  ggplot2::labs(
    title = "Incidence Versus Contact Heterogeneity",
    subtitle = "Point size shows simulated node population",
    x = "Node contact multiplier",
    y = "Modeled infection incidence\n(infections/person-year)"
  ) +
  theme_publication()

p_map <- ggplot2::ggplot(
  nodes,
  ggplot2::aes(x = longitude, y = latitude)
) +
  ggplot2::geom_point(
    ggplot2::aes(size = population_plot, fill = incidence_plot),
    shape = 21, color = "white", stroke = 0.3, alpha = 0.95
  ) +
  ggplot2::scale_fill_gradient(
    low = "#D9E8F5", high = "#B13E36",
    name = "Incidence"
  ) +
  ggplot2::scale_size_continuous(range = c(2, 6), name = "Population") +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = "Spatial Pattern Across Busia Sites",
    subtitle = "Color shows modeled infection incidence in children aged 6 months to <10 years",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_publication()

plot_paths <- list(
  search = save_plot(p_search, "01_search_trajectory", 6.8, 4.8),
  search_realized_eir = save_plot(
    p_search_realized,
    "01b_search_trajectory_realized_eir",
    6.8,
    4.8
  ),
  diagnostics = save_plot(p_metrics, "02_candidate_diagnostics", 6.8, 7.2),
  ranked_sites = save_plot(p_rank, "03_site_incidence_distribution", 6.8, 4.8),
  contact = save_plot(p_contact, "04_incidence_contact_eir", 6.8, 4.8),
  map = save_plot(p_map, "05_spatial_incidence_map", 6.4, 5.2)
)

if (!is.null(cf_nodes) && nrow(cf_nodes) > 0L) {
  cf_inc_col <- required_col(cf_nodes, c("infection_incidence_per_person_year"),
                             "counterfactual node incidence")
  cf_node_col <- required_col(cf_nodes, c("node_id", "node"), "counterfactual node id")
  cf_plot <- data.frame(
    node_id_plot = as.integer(cf_nodes[[cf_node_col]]),
    no_intervention = safe_num(cf_nodes[[cf_inc_col]])
  )
  paired <- merge(
    nodes[, c("node_id_plot", "incidence_plot")],
    cf_plot,
    by = "node_id_plot",
    all = FALSE
  )
  paired$node_rank <- rank(paired$no_intervention, ties.method = "first")
  paired_long <- rbind(
    data.frame(
      node_rank = paired$node_rank,
      arm = "With ITNs/ACT",
      incidence = paired$incidence_plot
    ),
    data.frame(
      node_rank = paired$node_rank,
      arm = "No interventions",
      incidence = paired$no_intervention
    )
  )
  paired_long$arm <- factor(paired_long$arm,
                            levels = c("No interventions", "With ITNs/ACT"))

  p_cf <- ggplot2::ggplot(
    paired_long,
    ggplot2::aes(x = node_rank, y = incidence, color = arm)
  ) +
    ggplot2::geom_hline(yintercept = target_mean, linetype = "22",
                        linewidth = 0.55, color = "#1F5A8A") +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::scale_color_manual(
      values = c("No interventions" = "#C75146", "With ITNs/ACT" = "#2F6F9F"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Selected Raw EIR Counterfactual",
      subtitle = "Same raw EIR and fixed delta; only ITNs/ACT differ",
      x = "Sites ordered by no-intervention incidence",
      y = "Modeled infection incidence\n(infections/person-year)"
    ) +
    theme_publication()

  plot_paths$counterfactual <- save_plot(
    p_cf,
    "06_selected_no_intervention_counterfactual",
    6.8,
    4.8
  )
}

if (requireNamespace("cowplot", quietly = TRUE)) {
  top_row <- cowplot::plot_grid(
    p_search,
    p_search_realized,
    labels = c("A", "B"),
    ncol = 2
  )
  middle_row <- cowplot::plot_grid(
    p_rank,
    p_contact,
    labels = c("C", "D"),
    ncol = 2
  )
  bottom_row <- cowplot::plot_grid(
    p_metrics,
    p_map,
    labels = c("E", "F"),
    ncol = 2
  )
  summary_fig <- cowplot::plot_grid(
    top_row,
    middle_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(1, 1, 1.25)
  )
  plot_paths$summary <- save_plot(
    summary_fig,
    "00_fixed_delta_calibration_summary",
    12,
    13.5
  )
} else {
  message("Package 'cowplot' is unavailable; skipped combined multi-panel figure.")
}

cat(sprintf("Read fixed-delta calibration outputs from: %s\n", out_dir))
cat(sprintf("Wrote fixed-delta calibration figures to: %s\n", plot_dir))
cat(sprintf(
  paste(
    "Summary: raw_init_EIR=%s, fixed_delta=%s, mean incidence=%s",
    "target=%s, SD=%s target_SD=%s, CV=%s target_CV=%s, mean realised EIR=%s\n"
  ),
  format_num(selected_eir, 4),
  format_num(selected_delta, 3),
  format_num(selected_mean, 3),
  format_num(target_mean, 3),
  format_num(selected_sd, 3),
  format_num(target_sd, 3),
  format_num(selected_cv, 3),
  format_num(target_cv, 3),
  format_num(selected_realized_eir, 3)
))
