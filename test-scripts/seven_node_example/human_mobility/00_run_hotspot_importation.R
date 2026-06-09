#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00_run_hotspot_importation.R
# ------------------------------------------------------------------------------
# Hotspot importation scenario for explicit human mobility.
#
# Purpose:
#   Test whether residents of low-transmission nodes travelling to a high-EIR
#   hotspot acquire infection risk away from home, with infections recorded back
#   in their home-node human records.
#
# Arms:
#   1. no_human_mobility
#   2. hotspot_importation_mobility
#
# Scenario:
#   - Node HOTSPOT_NODE has high init_EIR.
#   - Other nodes have lower init_EIR.
#   - Non-hotspot residents have probability Q_TO_HOTSPOT of sleeping in the
#     hotspot for one night.
#   - Hotspot residents stay home.
#   - No gene-drive release.
#   - Mosquito movement is disabled/identity to isolate human mobility.
#
# Expected:
#   - hotspot has high visitors_present.
#   - low-transmission origins may show increased infections under mobility.
#   - infection increases appear in travellers' home nodes, not only hotspot.
#   - no mosquito-genotype spread because no release and no mosquito movement.
#
# Outputs:
#   output/human_mobility_hotspot_importation/
#     timeseries_by_arm.csv
#     summary_by_arm.csv
#     node_comparison.csv
#     mobility_daily_by_arm.csv
#     human_move_probs_hotspot.csv
#     figure_hotspot_importation.png
#     figure_hotspot_importation.pdf
#     context.rds
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/seven_node_example/human_mobility/00_run_hotspot_importation.R",
    mustWork = FALSE
  )
}

example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
pkg_root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)

if (!requireNamespace("malariasimulationGD", quietly = TRUE)) {
  stop("Install malariasimulationGD first, e.g. Rscript test-scripts/install_local.R",
       call. = FALSE)
}
if (!requireNamespace("MGDrivE", quietly = TRUE)) {
  stop("Install MGDrivE first, e.g. Rscript test-scripts/install_local.R",
       call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 to render figures.", call. = FALSE)
}

source(file.path(example_dir, "lib", "movement_mu.R"))
source(file.path(example_dir, "lib", "synthetic_covariate.R"))
source(file.path(example_dir, "config", "seven_node_landscape.R"))
source(file.path(example_dir, "config", "movement.R"))
source(file.path(example_dir, "config", "seasonality.R"))
source(file.path(example_dir, "config", "covariate.R"))
source(file.path(example_dir, "config", "homing_drive.R"))
source(file.path(example_dir, "config", "trial_design.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

out_dir <- file.path(example_dir, "output", "human_mobility_hotspot_importation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Scenario configuration
# ------------------------------------------------------------------------------

RNG_SEED <- 20260514L

HOTSPOT_NODE <- 4L
INIT_EIR_LOW <- 5
INIT_EIR_HOTSPOT <- 80

Q_TO_HOTSPOT <- 0.02
TRIP_DURATION_TYPE <- "fixed"
TRIP_DURATION_MEAN <- 1L

# No release. Use a normal horizon.
td <- seven_node_trial_design()
TMAX_RUN <- as.integer(td$release_day + td$horizon_day)
READOUT_DAY <- TMAX_RUN

cat(sprintf("[%s] Hotspot importation mobility scenario\n", format(Sys.time(), "%F %T")))
cat(sprintf("  package root: %s\n", pkg_root))
cat(sprintf("  output:       %s\n", out_dir))
cat(sprintf("  hotspot node: %d\n", HOTSPOT_NODE))
cat(sprintf("  init_EIR low/hotspot: %.1f / %.1f\n", INIT_EIR_LOW, INIT_EIR_HOTSPOT))
cat(sprintf("  q to hotspot: %.3f\n", Q_TO_HOTSPOT))
cat("  gene-drive release: none\n")
cat("  mosquito movement: disabled/identity\n")

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

bind_rows_fill <- function(dfs) {
  all_names <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) d[[nm]] <- NA
    d[, all_names, drop = FALSE]
  })
  do.call(rbind, dfs)
}

safe_col_sum <- function(df, col) {
  if (col %in% names(df)) sum(df[[col]], na.rm = TRUE) else 0
}

mean_daily_sum <- function(ts_df, col) {
  if (!col %in% names(ts_df)) return(NA_real_)
  daily <- aggregate(
    ts_df[[col]],
    by = list(timestep = ts_df$timestep),
    FUN = sum,
    na.rm = TRUE
  )
  mean(daily$x, na.rm = TRUE)
}

mobility_daily_means <- function(ts_df) {
  data.frame(
    mean_daily_humans_present   = mean_daily_sum(ts_df, "humans_present"),
    mean_daily_visitors_present = mean_daily_sum(ts_df, "visitors_present"),
    mean_daily_residents_away   = mean_daily_sum(ts_df, "residents_away"),
    mean_daily_trips_started    = mean_daily_sum(ts_df, "trips_started")
  )
}

make_hotspot_P <- function(n_nodes, hotspot, q_to_hotspot) {
  stopifnot(hotspot >= 1, hotspot <= n_nodes)
  stopifnot(q_to_hotspot >= 0, q_to_hotspot <= 1)
  
  P <- diag(n_nodes)
  
  for (i in seq_len(n_nodes)) {
    P[i, ] <- 0
    if (i == hotspot) {
      P[i, i] <- 1
    } else {
      P[i, i] <- 1 - q_to_hotspot
      P[i, hotspot] <- q_to_hotspot
    }
  }
  
  stopifnot(max(abs(rowSums(P) - 1)) < 1e-12)
  rownames(P) <- colnames(P) <- as.character(seq_len(n_nodes))
  P
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
    suffixes = c("_no", "_mob")
  )
  
  rows <- lapply(numeric_cols, function(col) {
    x <- m[[paste0(col, "_no")]]
    y <- m[[paste0(col, "_mob")]]
    paired <- is.finite(x) & is.finite(y)
    
    if (!any(paired)) return(NULL)
    
    d <- abs(x[paired] - y[paired])
    
    data.frame(
      variable = col,
      n_paired = sum(paired),
      max_abs_diff = max(d),
      mean_abs_diff = mean(d),
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

make_node_comparison <- function(ts_all) {
  metric_cols <- intersect(
    c("n_infections", "visitors_present", "residents_away", "EIR_gamb", "FOIM_gamb"),
    names(ts_all)
  )
  
  sum_cols <- intersect(
    c("n_infections", "visitors_present", "residents_away"),
    metric_cols
  )
  
  node_sum <- aggregate(
    ts_all[, sum_cols, drop = FALSE],
    by = list(arm = ts_all$arm, node = ts_all$node),
    FUN = sum,
    na.rm = TRUE
  )
  
  mean_cols <- intersect(c("EIR_gamb", "FOIM_gamb"), metric_cols)
  if (length(mean_cols) > 0L) {
    node_mean <- aggregate(
      ts_all[, mean_cols, drop = FALSE],
      by = list(arm = ts_all$arm, node = ts_all$node),
      FUN = mean,
      na.rm = TRUE
    )
    node_summary <- merge(node_sum, node_mean, by = c("arm", "node"), all = TRUE)
  } else {
    node_summary <- node_sum
  }
  
  node_no <- node_summary[node_summary$arm == "no_human_mobility", ]
  node_mob <- node_summary[node_summary$arm == "hotspot_importation_mobility", ]
  
  cmp <- merge(node_no, node_mob, by = "node", suffixes = c("_no", "_mob"))
  
  cmp$infection_diff <- cmp$n_infections_mob - cmp$n_infections_no
  cmp$infection_rel_diff <- cmp$infection_diff / pmax(cmp$n_infections_no, 1)
  
  if (all(c("visitors_present_mob", "residents_away_mob") %in% names(cmp))) {
    cmp$net_visitors <- cmp$visitors_present_mob - cmp$residents_away_mob
  }
  
  cmp
}

get_xy_cols <- function(nodes) {
  candidates <- list(
    c("x", "y"),
    c("X", "Y"),
    c("longitude", "latitude"),
    c("lon", "lat")
  )
  
  for (cc in candidates) {
    if (all(cc %in% names(nodes))) {
      return(cc)
    }
  }
  
  numeric_cols <- names(nodes)[vapply(nodes, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("node", "NH", "NH_per_node", "population"))
  
  if (length(numeric_cols) >= 2L) {
    return(numeric_cols[1:2])
  }
  
  stop("Could not infer coordinate columns in `nodes`.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Build landscape and parameters
# ------------------------------------------------------------------------------

set.seed(RNG_SEED)

land <- build_seven_node_landscape()
n_nodes <- land$n_nodes
nodes <- land$nodes
D <- land$D

if (!"node" %in% names(nodes)) {
  nodes$node <- seq_len(n_nodes)
}

cov <- build_seven_node_covariate(D)

contact_surface <- list(
  type = "contact_surface",
  contact_multiplier = stats::setNames(
    as.numeric(cov$contact_multiplier),
    as.character(seq_len(n_nodes))
  )
)

seas <- seven_node_seasonality()
cube <- build_seven_node_drive_cube()

P_hotspot <- make_hotspot_P(
  n_nodes = n_nodes,
  hotspot = HOTSPOT_NODE,
  q_to_hotspot = Q_TO_HOTSPOT
)

cat("\nHuman hotspot movement matrix P:\n")
print(P_hotspot)

cat("\nMovement row summaries:\n")
print(data.frame(
  node = seq_len(n_nodes),
  stay_home = diag(P_hotspot),
  to_hotspot = P_hotspot[, HOTSPOT_NODE],
  row_sum = rowSums(P_hotspot)
))

# Explicit mobility cannot be combined with import/export transmission mixing.
mixing_identity <- list(diag(n_nodes))
p_captured_zero <- list(matrix(0, n_nodes, n_nodes))

# Disable mosquito movement to isolate human movement.
mosquito_move_probs_identity <- diag(n_nodes)
mosquito_move_rates_zero <- rep(0, n_nodes)

arms <- list(
  no_human_mobility = list(
    human_mobility_enabled = FALSE,
    human_move_probs = NULL,
    human_mobility_store_diagnostics = FALSE
  ),
  hotspot_importation_mobility = list(
    human_mobility_enabled = TRUE,
    human_move_probs = P_hotspot,
    human_trip_duration_type = TRIP_DURATION_TYPE,
    human_trip_duration_mean = TRIP_DURATION_MEAN,
    human_mobility_store_diagnostics = TRUE
  )
)

build_one_node_params <- function(node_idx, arm_cfg) {
  NH_v <- nodes$NH_per_node[node_idx]
  
  init_eir_node <- if (node_idx == HOTSPOT_NODE) {
    INIT_EIR_HOTSPOT
  } else {
    INIT_EIR_LOW
  }
  
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
    init_EIR = init_eir_node
  )
  
  # Attach cube so the native genotype backend remains active, but do not release.
  p$cube <- cube
  
  # No mosquito movement.
  p$move_probs <- mosquito_move_probs_identity
  p$move_rates <- mosquito_move_rates_zero
  
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
  
  # Intentionally no gene-drive releases in this scenario.
  p
}

build_params_list <- function(arm_cfg) {
  lapply(seq_len(n_nodes), build_one_node_params, arm_cfg = arm_cfg)
}

# ------------------------------------------------------------------------------
# Extract outputs
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

summarise_arm <- function(arm, ts_df, elapsed) {
  mobility_means <- mobility_daily_means(ts_df)
  
  data.frame(
    arm = arm,
    tmax = TMAX_RUN,
    hotspot_node = HOTSPOT_NODE,
    init_eir_low = INIT_EIR_LOW,
    init_eir_hotspot = INIT_EIR_HOTSPOT,
    q_to_hotspot = if (arm == "hotspot_importation_mobility") Q_TO_HOTSPOT else NA_real_,
    n_nodes = n_nodes,
    readout_day = READOUT_DAY,
    regional_infections_total = regional_infections_total(ts_df),
    regional_pfpr_readout = regional_prevalence_at(ts_df, READOUT_DAY),
    mean_daily_humans_present =
      mobility_means$mean_daily_humans_present,
    mean_daily_visitors_present =
      mobility_means$mean_daily_visitors_present,
    mean_daily_residents_away =
      mobility_means$mean_daily_residents_away,
    mean_daily_trips_started =
      mobility_means$mean_daily_trips_started,
    elapsed_seconds = round(elapsed, 1),
    stringsAsFactors = FALSE
  )
}

run_arm <- function(arm, arm_cfg) {
  cat(sprintf("\n[%s] Running arm: %s\n", format(Sys.time(), "%F %T"), arm))
  
  params_list <- build_params_list(arm_cfg)
  
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
  
  list(
    sim = sim,
    timeseries = ts_df,
    summary = summarise_arm(arm, ts_df, elapsed)
  )
}

# ------------------------------------------------------------------------------
# Run scenario
# ------------------------------------------------------------------------------

results <- list()
for (arm in names(arms)) {
  results[[arm]] <- run_arm(arm, arms[[arm]])
}

timeseries_all <- bind_rows_fill(lapply(results, `[[`, "timeseries"))
summary_all <- bind_rows_fill(lapply(results, `[[`, "summary"))

comparison_summary <- compare_arms(
  timeseries_all,
  "no_human_mobility",
  "hotspot_importation_mobility"
)

node_comparison <- make_node_comparison(timeseries_all)

mobility_cols <- intersect(
  c("humans_present", "visitors_present", "residents_away", "trips_started"),
  names(timeseries_all)
)

mobility_daily_by_arm <- if (length(mobility_cols) > 0L) {
  aggregate(
    timeseries_all[, mobility_cols, drop = FALSE],
    by = list(arm = timeseries_all$arm, timestep = timeseries_all$timestep),
    FUN = sum,
    na.rm = TRUE
  )
} else {
  data.frame()
}

# Regional daily infections.
regional_infections_daily <- aggregate(
  n_infections ~ arm + timestep,
  data = timeseries_all,
  FUN = sum,
  na.rm = TRUE
)

# Regional PfPR.
if (all(c("n_detect_lm_730_3650", "n_age_730_3650") %in% names(timeseries_all))) {
  regional_pfpr_daily <- aggregate(
    timeseries_all[, c("n_detect_lm_730_3650", "n_age_730_3650")],
    by = list(arm = timeseries_all$arm, timestep = timeseries_all$timestep),
    FUN = sum,
    na.rm = TRUE
  )
  regional_pfpr_daily$pfpr_2_10_lm <-
    regional_pfpr_daily$n_detect_lm_730_3650 /
    pmax(regional_pfpr_daily$n_age_730_3650, 1)
} else {
  regional_pfpr_daily <- data.frame()
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))

xy_cols <- get_xy_cols(nodes)
x_col <- xy_cols[1]
y_col <- xy_cols[2]

nodes_plot <- nodes
nodes_plot$is_hotspot <- nodes_plot$node == HOTSPOT_NODE
nodes_plot$label <- paste0("Node ", nodes_plot$node)

# Arrows from non-hotspot origins to hotspot.
arrow_df <- data.frame()
for (i in seq_len(n_nodes)) {
  if (i != HOTSPOT_NODE && P_hotspot[i, HOTSPOT_NODE] > 0) {
    origin <- nodes_plot[nodes_plot$node == i, ]
    dest <- nodes_plot[nodes_plot$node == HOTSPOT_NODE, ]
    arrow_df <- rbind(arrow_df, data.frame(
      x = origin[[x_col]],
      y = origin[[y_col]],
      xend = dest[[x_col]],
      yend = dest[[y_col]],
      prob = P_hotspot[i, HOTSPOT_NODE]
    ))
  }
}

pA <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = arrow_df,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = prob),
    arrow = ggplot2::arrow(length = grid::unit(0.12, "inches")),
    alpha = 0.45
  ) +
  ggplot2::geom_point(
    data = nodes_plot,
    ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], shape = is_hotspot),
    size = 4
  ) +
  ggplot2::geom_text(
    data = nodes_plot,
    ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], label = node),
    vjust = -1,
    size = 3.5
  ) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  ggplot2::labs(
    title = "A. Hotspot importation design",
    subtitle = sprintf("Non-hotspot residents sleep in node %d with probability %.2f", HOTSPOT_NODE, Q_TO_HOTSPOT),
    x = x_col,
    y = y_col,
    linewidth = "P to hotspot",
    shape = "Hotspot"
  ) +
  ggplot2::theme(legend.position = "bottom")

# Node-level infection difference.
plot_node_cmp <- node_comparison
plot_node_cmp$is_hotspot <- plot_node_cmp$node == HOTSPOT_NODE

pB <- ggplot2::ggplot(
  plot_node_cmp,
  ggplot2::aes(x = factor(node), y = infection_diff, fill = is_hotspot)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
    labels = c(`FALSE` = "Other nodes", `TRUE` = "Hotspot")
  ) +
  ggplot2::labs(
    title = "B. Node-level infection change",
    subtitle = "Mobility arm minus no-mobility arm",
    x = "Home node",
    y = "Difference in cumulative infections",
    fill = NULL
  ) +
  ggplot2::theme(legend.position = "bottom")

# Visitors and residents away.
mob_long <- data.frame(
  node = rep(plot_node_cmp$node, 2),
  metric = rep(c("Visitors present", "Residents away"), each = nrow(plot_node_cmp)),
  value = c(plot_node_cmp$visitors_present_mob, plot_node_cmp$residents_away_mob)
)

pC <- ggplot2::ggplot(
  mob_long,
  ggplot2::aes(x = factor(node), y = value, fill = metric)
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::labs(
    title = "C. Mobility accounting by node",
    subtitle = "Cumulative over the simulation",
    x = "Node",
    y = "Person-nights",
    fill = NULL
  ) +
  ggplot2::theme(legend.position = "bottom")

# Regional 14-day infection difference.
window_days <- 14L
regional_wide <- reshape(
  regional_infections_daily,
  idvar = "timestep",
  timevar = "arm",
  direction = "wide"
)

no_col <- "n_infections.no_human_mobility"
mob_col <- "n_infections.hotspot_importation_mobility"

if (all(c(no_col, mob_col) %in% names(regional_wide))) {
  regional_wide$diff <- regional_wide[[mob_col]] - regional_wide[[no_col]]
  regional_wide$window_end_day <-
    ((regional_wide$timestep - 1L) %/% window_days) * window_days + window_days
  
  regional_window <- aggregate(
    diff ~ window_end_day,
    data = regional_wide,
    FUN = sum,
    na.rm = TRUE
  )
} else {
  regional_window <- data.frame(window_end_day = integer(), diff = numeric())
}

pD <- ggplot2::ggplot(
  regional_window,
  ggplot2::aes(x = window_end_day, y = diff)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 1.2) +
  ggplot2::labs(
    title = "D. Regional infection difference over time",
    subtitle = "14-day windows: mobility minus no mobility",
    x = "Window end day",
    y = "Difference in infections"
  )

save_four_panel <- function(filename, width = 12, height = 9) {
  grDevices::png(filename, width = width, height = height, units = "in", res = 150)
  grid::grid.newpage()
  lay <- grid::grid.layout(2, 2)
  grid::pushViewport(grid::viewport(layout = lay))
  
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(pC, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(pD, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  grid::popViewport()
  grDevices::dev.off()
}

save_four_panel_pdf <- function(filename, width = 12, height = 9) {
  grDevices::pdf(filename, width = width, height = height)
  grid::grid.newpage()
  lay <- grid::grid.layout(2, 2)
  grid::pushViewport(grid::viewport(layout = lay))
  
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(pC, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(pD, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  grid::popViewport()
  grDevices::dev.off()
}

png_path <- file.path(out_dir, "figure_hotspot_importation.png")
pdf_path <- file.path(out_dir, "figure_hotspot_importation.pdf")
save_four_panel(png_path)
save_four_panel_pdf(pdf_path)

# ------------------------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------------------------

utils::write.csv(
  timeseries_all,
  file.path(out_dir, "timeseries_by_arm.csv"),
  row.names = FALSE
)

utils::write.csv(
  summary_all,
  file.path(out_dir, "summary_by_arm.csv"),
  row.names = FALSE
)

utils::write.csv(
  node_comparison,
  file.path(out_dir, "node_comparison.csv"),
  row.names = FALSE
)

utils::write.csv(
  comparison_summary,
  file.path(out_dir, "comparison_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  mobility_daily_by_arm,
  file.path(out_dir, "mobility_daily_by_arm.csv"),
  row.names = FALSE
)

utils::write.csv(
  regional_infections_daily,
  file.path(out_dir, "regional_infections_daily.csv"),
  row.names = FALSE
)

utils::write.csv(
  regional_pfpr_daily,
  file.path(out_dir, "regional_pfpr_daily.csv"),
  row.names = FALSE
)

utils::write.csv(
  P_hotspot,
  file.path(out_dir, "human_move_probs_hotspot.csv"),
  row.names = TRUE
)

utils::write.csv(
  nodes,
  file.path(out_dir, "nodes.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    scenario = "hotspot_importation",
    hotspot_node = HOTSPOT_NODE,
    init_eir_low = INIT_EIR_LOW,
    init_eir_hotspot = INIT_EIR_HOTSPOT,
    q_to_hotspot = Q_TO_HOTSPOT,
    P_hotspot = P_hotspot,
    mosquito_movement = "identity/no movement",
    gene_drive_release = "none",
    readout_day = READOUT_DAY,
    figure_png = png_path,
    figure_pdf = pdf_path
  ),
  file.path(out_dir, "context.rds")
)

# ------------------------------------------------------------------------------
# Console summary
# ------------------------------------------------------------------------------

cat(sprintf("\n[%s] Done. Outputs written to: %s\n",
            format(Sys.time(), "%F %T"), out_dir))

cat("\nSummary by arm:\n")
print(summary_all)

cat("\nNode-level comparison:\n")
print(node_comparison[, intersect(
  c(
    "node",
    "n_infections_no",
    "n_infections_mob",
    "infection_diff",
    "infection_rel_diff",
    "visitors_present_mob",
    "residents_away_mob",
    "net_visitors",
    "EIR_gamb_no",
    "EIR_gamb_mob"
  ),
  names(node_comparison)
), drop = FALSE])

cat("\nMovement conservation check:\n")
cat(sprintf(
  "  total visitors_present = %.0f\n",
  sum(node_comparison$visitors_present_mob, na.rm = TRUE)
))
cat(sprintf(
  "  total residents_away   = %.0f\n",
  sum(node_comparison$residents_away_mob, na.rm = TRUE)
))
cat(sprintf(
  "  difference             = %.0f\n",
  sum(node_comparison$visitors_present_mob, na.rm = TRUE) -
    sum(node_comparison$residents_away_mob, na.rm = TRUE)
))

cat(sprintf("\nSaved figure:\n  %s\n  %s\n", png_path, pdf_path))

