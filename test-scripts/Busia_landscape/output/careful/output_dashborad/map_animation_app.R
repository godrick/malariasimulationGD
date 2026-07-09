required_packages <- c("shiny", "ggplot2", "RColorBrewer", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install required packages before running this app: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(shiny)
library(ggplot2)

`%||%` <- function(x, y) if (is.null(x)) y else x

find_data_dir <- function() {
  candidates <- unique(normalizePath(c(
    file.path(getwd(), ".."),
    getwd(),
    file.path(
      getwd(),
      "test-scripts", "Busia_landscape", "output", "careful"
    )
  ), mustWork = FALSE))

  for (candidate in candidates) {
    needed <- c("no_release_timeseries.csv", "release_timeseries.csv", "nodes.csv")
    if (all(file.exists(file.path(candidate, needed)))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    "Could not find no_release_timeseries.csv, release_timeseries.csv, and nodes.csv. ",
    "Run this app from output_dashborad/ or the repository root.",
    call. = FALSE
  )
}

load_optional_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) NULL)
}

read_timeseries <- function(data_dir) {
  no_release <- read.csv(
    file.path(data_dir, "no_release_timeseries.csv"),
    check.names = FALSE
  )
  release <- read.csv(
    file.path(data_dir, "release_timeseries.csv"),
    check.names = FALSE
  )
  out <- rbind(no_release, release)
  out$node <- as.integer(out$node)
  out$arm <- factor(out$arm, levels = c("no_release", "release"))
  out
}

metric_labels <- c(
  n_infections = "New infections",
  n_inc_clinical_182_5475 = "Clinical cases, age 6 months to 15 years",
  p_inc_clinical_182_5475 = "Clinical incidence probability, age 6 months to 15 years",
  infectivity = "Human infectivity",
  infectivity_weighted_I_gamb = "Weighted human infectivity to An. gambiae",
  vector_infectivity_mean_gamb = "Mean vector infectivity, An. gambiae",
  EIR_gamb = "Entomological inoculation rate, An. gambiae",
  n_bitten = "People bitten",
  FOIM_gamb = "Force of infection on mosquitoes, An. gambiae",
  mu_gamb = "Mosquito movement/mortality value",
  S_count = "Susceptible people",
  A_count = "Asymptomatic infections",
  D_count = "Clinical disease",
  U_count = "Subpatent infections",
  Tr_count = "Treated people",
  ica_mean = "Mean acquired clinical immunity",
  icm_mean = "Mean maternal immunity",
  ib_mean = "Mean anti-parasite immunity",
  iva_mean = "Mean variant-specific immunity",
  ivm_mean = "Mean maternal variant immunity",
  id_mean = "Mean detection immunity",
  n_age_0_3650 = "Population age 0 to 10 years",
  n_age_182_5475 = "Population age 6 months to 15 years",
  n_age_730_3650 = "Population age 2 to 10 years",
  n_age_3650_7300 = "Population age 10 to 20 years",
  n_age_7300_10950 = "Population age 20 to 30 years",
  n_age_10950_14600 = "Population age 30 to 40 years",
  n_age_14600_18250 = "Population age 40 to 50 years",
  n_age_18250_21900 = "Population age 50 to 60 years",
  n_age_21900_25550 = "Population age 60 to 70 years",
  n_age_25550_29200 = "Population age 70 to 80 years",
  n_age_29200_73000 = "Population age 80+ years",
  E_gamb_count = "Egg-stage mosquitoes, An. gambiae",
  L_gamb_count = "Larval mosquitoes, An. gambiae",
  P_gamb_count = "Pupal mosquitoes, An. gambiae",
  Sm_gamb_count = "Susceptible adult mosquitoes, An. gambiae",
  Pm_gamb_count = "Exposed adult mosquitoes, An. gambiae",
  Im_gamb_count = "Infectious adult mosquitoes, An. gambiae",
  total_M_gamb = "Total adult mosquitoes, An. gambiae",
  n_detect_lm_730_3650 = "Microscopy-positive children, age 2 to 10 years",
  p_detect_lm_730_3650 = "PfPR by microscopy, age 2 to 10 years",
  n_detect_pcr_730_3650 = "PCR-positive children, age 2 to 10 years",
  natural_deaths = "Natural deaths"
)

point_in_polygon <- function(x, y, poly_x, poly_y) {
  inside <- rep(FALSE, length(x))
  j <- length(poly_x)
  for (i in seq_along(poly_x)) {
    crosses <- ((poly_y[i] > y) != (poly_y[j] > y)) &
      (x < (poly_x[j] - poly_x[i]) * (y - poly_y[i]) /
        (poly_y[j] - poly_y[i]) + poly_x[i])
    inside <- xor(inside, crosses)
    j <- i
  }
  inside
}

make_hull <- function(nodes, margin_fraction = 0.08) {
  hull_index <- grDevices::chull(nodes$x, nodes$y)
  center <- c(mean(nodes$x), mean(nodes$y))
  hull <- nodes[hull_index, c("x", "y")]
  hull$x <- center[[1L]] + (hull$x - center[[1L]]) * (1 + margin_fraction)
  hull$y <- center[[2L]] + (hull$y - center[[2L]]) * (1 + margin_fraction)
  hull <- rbind(hull, hull[1L, ])
  hull
}

make_idw_surface <- function(nodes, value_col = "value", grid_size = 72L,
                             power = 2, margin_fraction = 0.08) {
  hull <- make_hull(nodes, margin_fraction = margin_fraction)
  x_range <- range(hull$x, na.rm = TRUE)
  y_range <- range(hull$y, na.rm = TRUE)
  grid <- expand.grid(
    x = seq(x_range[[1L]], x_range[[2L]], length.out = grid_size),
    y = seq(y_range[[1L]], y_range[[2L]], length.out = grid_size)
  )
  grid$inside <- point_in_polygon(grid$x, grid$y, hull$x, hull$y)

  values <- nodes[[value_col]]
  valid <- is.finite(values)
  if (!any(valid)) {
    grid$value <- NA_real_
    return(grid[grid$inside, c("x", "y", "value")])
  }

  xs <- nodes$x[valid]
  ys <- nodes$y[valid]
  values <- values[valid]
  estimate <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    distances <- sqrt((grid$x[[i]] - xs)^2 + (grid$y[[i]] - ys)^2)
    exact <- which(distances < 1e-9)
    if (length(exact) > 0L) {
      estimate[[i]] <- values[[exact[[1L]]]]
    } else {
      weights <- 1 / (distances^power)
      estimate[[i]] <- sum(weights * values) / sum(weights)
    }
  }

  grid$value <- estimate
  grid[grid$inside, c("x", "y", "value")]
}

data_dir <- find_data_dir()
timeseries <- read_timeseries(data_dir)
nodes <- read.csv(file.path(data_dir, "nodes.csv"), check.names = FALSE)
context <- load_optional_rds(file.path(data_dir, "context.rds"))
release_nodes <- as.integer(context$release_nodes %||% integer())
if (length(release_nodes) == 0L) {
  release_nodes <- integer()
}
nodes$release_site <- nodes$node %in% release_nodes

numeric_columns <- names(timeseries)[vapply(timeseries, is.numeric, logical(1))]
metric_columns <- setdiff(numeric_columns, c("node", "timestep"))
default_metric <- intersect(c("p_detect_lm_730_3650", "EIR_gamb"), metric_columns)
if (length(default_metric) == 0L) {
  default_metric <- metric_columns[[1L]]
} else {
  default_metric <- default_metric[[1L]]
}

metric_choices <- stats::setNames(
  metric_columns,
  ifelse(
    metric_columns %in% names(metric_labels),
    paste0(metric_labels[metric_columns], " (", metric_columns, ")"),
    metric_columns
  )
)

palette_info <- RColorBrewer::brewer.pal.info
palette_choices <- rownames(palette_info[palette_info$category %in% c("seq", "div"), ])

ui <- fluidPage(
  tags$head(tags$style(HTML("
    :root {
      --bg: #f7f8fb;
      --panel: #ffffff;
      --ink: #22262f;
      --muted: #697282;
      --line: #dce2ec;
      --accent: #356f72;
      --accent-soft: #e5f0ef;
    }
    body {
      background: var(--bg);
      color: var(--ink);
      font-family: Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    }
    .app-shell {
      max-width: 1500px;
      margin: 18px auto 34px auto;
      padding: 0 18px;
    }
    .topbar {
      display: flex;
      justify-content: space-between;
      align-items: flex-end;
      gap: 16px;
      margin-bottom: 16px;
    }
    .title {
      font-size: 26px;
      font-weight: 780;
      letter-spacing: 0;
      margin: 0;
    }
    .subtitle {
      color: var(--muted);
      margin-top: 4px;
    }
    .control-panel, .map-panel, .summary-card {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      box-shadow: 0 10px 28px rgba(34, 38, 47, 0.06);
    }
    .control-panel {
      padding: 18px;
    }
    .map-panel {
      padding: 16px 18px 10px 18px;
    }
    .summary-grid {
      display: grid;
      grid-template-columns: repeat(4, minmax(150px, 1fr));
      gap: 12px;
      margin-bottom: 16px;
    }
    .summary-card {
      padding: 13px 15px;
    }
    .summary-label {
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: .06em;
      margin-bottom: 4px;
    }
    .summary-value {
      font-size: 21px;
      font-weight: 750;
    }
    .section-title {
      font-size: 17px;
      font-weight: 750;
      margin: 0 0 12px 0;
    }
    label {
      color: #343a46;
      font-weight: 650;
      font-size: 13px;
    }
    .selectize-input, .form-control {
      border-radius: 7px;
      border-color: #cfd7e5;
      box-shadow: none;
    }
    .selectize-input.focus, .form-control:focus {
      border-color: var(--accent);
      box-shadow: 0 0 0 3px rgba(53, 111, 114, .12);
    }
    .btn, .btn-default {
      border-radius: 7px;
    }
    .help-text {
      color: var(--muted);
      font-size: 12px;
      line-height: 1.45;
      margin-top: -3px;
      margin-bottom: 14px;
    }
    .irs-bar, .irs-single {
      background: var(--accent);
      border-color: var(--accent);
    }
    .irs-from, .irs-to {
      background: var(--accent);
    }
    @media (max-width: 900px) {
      .topbar {
        align-items: flex-start;
        flex-direction: column;
      }
      .summary-grid {
        grid-template-columns: repeat(2, minmax(140px, 1fr));
      }
    }
    @media (max-width: 600px) {
      .summary-grid {
        grid-template-columns: 1fr;
      }
    }
  "))),
  div(
    class = "app-shell",
    div(
      class = "topbar",
      div(
        h1(class = "title", "Busia Metric Map Animation"),
        div(class = "subtitle", "Use the play button on the day slider to watch node-level values move over time.")
      ),
      downloadButton("download_frame", "Download current frame")
    ),
    div(
      class = "summary-grid",
      div(class = "summary-card", div(class = "summary-label", "Scenario"), div(class = "summary-value", textOutput("scenario_text", inline = TRUE))),
      div(class = "summary-card", div(class = "summary-label", "Day"), div(class = "summary-value", textOutput("day_text", inline = TRUE))),
      div(class = "summary-card", div(class = "summary-label", "Node mean"), div(class = "summary-value", textOutput("mean_text", inline = TRUE))),
      div(class = "summary-card", div(class = "summary-label", "Release-site mean"), div(class = "summary-value", textOutput("release_mean_text", inline = TRUE)))
    ),
    fluidRow(
      column(
        width = 3,
        div(
          class = "control-panel",
          div(class = "section-title", "Map controls"),
          selectInput(
            "arm",
            "Scenario",
            choices = c("No release" = "no_release", "Release" = "release"),
            selected = "release"
          ),
          selectizeInput(
            "metric",
            "Metric",
            choices = metric_choices,
            selected = default_metric,
            multiple = FALSE
          ),
          sliderInput(
            "timestep",
            "Day",
            min = min(timeseries$timestep, na.rm = TRUE),
            max = max(timeseries$timestep, na.rm = TRUE),
            value = min(timeseries$timestep, na.rm = TRUE),
            step = 1,
            sep = ",",
            animate = animationOptions(interval = 112, loop = TRUE)
          ),
          selectInput(
            "palette",
            "ColorBrewer palette",
            choices = palette_choices,
            selected = "YlOrRd"
          ),
          checkboxInput("reverse_palette", "Reverse palette", value = FALSE),
          radioButtons(
            "scale_mode",
            "Color range",
            choices = c("Fixed over selected scenario" = "scenario", "Current day only" = "day"),
            selected = "scenario"
          ),
          checkboxInput("show_surface", "Show interpolated heatmap", value = TRUE),
          sliderInput("surface_alpha", "Heatmap opacity", min = 0.25, max = 1, value = 0.82, step = 0.05),
          sliderInput("grid_size", "Heatmap detail", min = 35, max = 100, value = 72, step = 5),
          div(
            class = "help-text",
            "The heatmap is an inverse-distance interpolation over node coordinates, clipped to the node hull."
          )
        )
      ),
      column(
        width = 9,
        div(
          class = "map-panel",
          plotOutput("map_plot", height = "720px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$arm, {
    available_days <- sort(unique(timeseries$timestep[as.character(timeseries$arm) == input$arm]))
    updateSliderInput(
      session,
      "timestep",
      min = min(available_days),
      max = max(available_days),
      value = min(max(input$timestep, min(available_days)), max(available_days))
    )
  })

  metric_label <- reactive({
    metric <- input$metric %||% default_metric
    if (metric %in% names(metric_labels)) {
      return(unname(metric_labels[[metric]]))
    }
    metric
  })

  scenario_metric_data <- reactive({
    req(input$arm, input$metric)
    timeseries[
      as.character(timeseries$arm) == input$arm,
      c("node", "arm", "timestep", input$metric),
      drop = FALSE
    ]
  })

  current_nodes <- reactive({
    df <- scenario_metric_data()
    day <- input$timestep %||% min(df$timestep, na.rm = TRUE)
    day <- min(max(day, min(df$timestep, na.rm = TRUE)), max(df$timestep, na.rm = TRUE))
    day_df <- df[df$timestep == day, , drop = FALSE]
    names(day_df)[names(day_df) == input$metric] <- "value"
    out <- merge(nodes, day_df[, c("node", "value"), drop = FALSE], by = "node", all.x = TRUE)
    out$release_site <- out$node %in% release_nodes
    out
  })

  color_limits <- reactive({
    req(input$metric)
    if (identical(input$scale_mode, "day")) {
      values <- current_nodes()$value
    } else {
      values <- scenario_metric_data()[[input$metric]]
    }
    values <- values[is.finite(values)]
    if (length(values) == 0L) {
      return(c(0, 1))
    }
    limits <- range(values, na.rm = TRUE)
    if (isTRUE(all.equal(limits[[1L]], limits[[2L]]))) {
      limits <- limits + c(-0.5, 0.5)
    }
    limits
  })

  surface_data <- reactive({
    req(input$show_surface)
    make_idw_surface(
      current_nodes(),
      value_col = "value",
      grid_size = as.integer(input$grid_size %||% 72L)
    )
  })

  summary_stats <- reactive({
    dat <- current_nodes()
    values <- dat$value
    release_values <- dat$value[dat$release_site]
    list(
      mean = mean(values, na.rm = TRUE),
      release_mean = if (length(release_values) > 0L) mean(release_values, na.rm = TRUE) else NA_real_
    )
  })

  output$scenario_text <- renderText({
    if (identical(input$arm, "release")) "Release" else "No release"
  })

  output$day_text <- renderText(scales::comma(input$timestep %||% 0))

  output$mean_text <- renderText({
    value <- summary_stats()$mean
    if (!is.finite(value)) "NA" else scales::comma(value, accuracy = 0.001)
  })

  output$release_mean_text <- renderText({
    value <- summary_stats()$release_mean
    if (!is.finite(value)) "NA" else scales::comma(value, accuracy = 0.001)
  })

  build_map <- reactive({
    dat <- current_nodes()
    hull <- make_hull(nodes)
    limits <- color_limits()
    palette_name <- input$palette %||% "YlOrRd"
    palette_values <- RColorBrewer::brewer.pal(
      RColorBrewer::brewer.pal.info[palette_name, "maxcolors"],
      palette_name
    )
    if (isTRUE(input$reverse_palette)) {
      palette_values <- rev(palette_values)
    }

    p <- ggplot() +
      geom_polygon(
        data = hull,
        aes(x = x, y = y),
        fill = "#eef2f6",
        color = "#cbd3df",
        linewidth = 0.45
      )

    if (isTRUE(input$show_surface)) {
      p <- p +
        geom_raster(
          data = surface_data(),
          aes(x = x, y = y, fill = value),
          alpha = input$surface_alpha %||% 0.82,
          interpolate = TRUE
        )
    }

    p +
      geom_point(
        data = dat[!dat$release_site, , drop = FALSE],
        aes(x = x, y = y, fill = value),
        shape = 21,
        size = 3.7,
        stroke = 0.55,
        color = "white",
        alpha = 0.98
      ) +
      geom_point(
        data = dat[dat$release_site, , drop = FALSE],
        aes(x = x, y = y, fill = value),
        shape = 24,
        size = 5.3,
        stroke = 0.85,
        color = "#2a2f3a",
        alpha = 0.98
      ) +
      geom_point(
        data = dat[dat$release_site, , drop = FALSE],
        aes(x = x, y = y),
        shape = 24,
        size = 6.2,
        stroke = 0.45,
        color = "#2a2f3a",
        fill = NA,
        alpha = 0.35
      ) +
      scale_fill_gradientn(
        colours = palette_values,
        limits = limits,
        oob = scales::squish,
        labels = scales::comma,
        name = metric_label()
      ) +
      coord_equal(expand = FALSE) +
      labs(
        title = paste0(metric_label(), " over Busia nodes"),
        subtitle = sprintf(
          "%s scenario | day %s | release sites shown as triangles",
          if (identical(input$arm, "release")) "Release" else "No-release",
          scales::comma(input$timestep %||% 0)
        ),
        x = "Projected east-west distance from landscape center (km)",
        y = "Projected north-south distance from landscape center (km)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 20, color = "#22262f"),
        plot.subtitle = element_text(color = "#697282", margin = margin(b = 12)),
        axis.title = element_text(color = "#4b5563"),
        axis.text = element_text(color = "#697282"),
        panel.grid.major = element_line(color = "#e5e9f0", linewidth = 0.35),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
  })

  output$map_plot <- renderPlot({
    build_map()
  }, res = 125)

  output$download_frame <- downloadHandler(
    filename = function() {
      sprintf(
        "busia_map_%s_%s_day_%s.png",
        input$arm %||% "scenario",
        input$metric %||% "metric",
        input$timestep %||% "day"
      )
    },
    content = function(file) {
      ggplot2::ggsave(
        filename = file,
        plot = build_map(),
        width = 11,
        height = 8.2,
        dpi = 180
      )
    }
  )
}

shinyApp(ui, server)
