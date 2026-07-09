required_packages <- c("shiny", "ggplot2", "RColorBrewer", "DT", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install required packages before running this dashboard: ",
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
    if (all(file.exists(file.path(
      candidate,
      c("no_release_timeseries.csv", "release_timeseries.csv")
    )))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    "Could not find no_release_timeseries.csv and release_timeseries.csv. ",
    "Run the app from output_dashborad/ or the repository root.",
    call. = FALSE
  )
}

data_dir <- find_data_dir()

load_optional_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) NULL)
}

format_value <- function(value) {
  if (is.null(value)) {
    return("NA")
  }
  if (is.logical(value) && length(value) == 1L) {
    return(ifelse(isTRUE(value), "TRUE", "FALSE"))
  }
  if (is.numeric(value)) {
    value <- as.numeric(value)
    if (length(value) == 1L) {
      return(scales::comma(value, accuracy = ifelse(abs(value) >= 10, 1, 0.001)))
    }
    if (length(value) <= 8L) {
      return(paste(scales::comma(value, accuracy = 0.001), collapse = ", "))
    }
    return(sprintf(
      "%s values; range %s to %s",
      scales::comma(length(value)),
      scales::comma(min(value, na.rm = TRUE), accuracy = 0.001),
      scales::comma(max(value, na.rm = TRUE), accuracy = 0.001)
    ))
  }
  if (is.character(value)) {
    if (length(value) <= 6L) {
      return(paste(value, collapse = ", "))
    }
    return(sprintf("%s values", scales::comma(length(value))))
  }
  if (is.data.frame(value)) {
    return(sprintf("%s rows x %s columns", scales::comma(nrow(value)), ncol(value)))
  }
  if (is.list(value)) {
    return(sprintf("list with %s entries", length(value)))
  }
  paste(value, collapse = ", ")
}

read_timeseries <- function() {
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

timeseries <- read_timeseries()
context <- load_optional_rds(file.path(data_dir, "context.rds"))
theta <- load_optional_rds(file.path(data_dir, "theta.rds")) %||% context$theta
calibration <- load_optional_rds(file.path(data_dir, "calibrated_init_eir.rds"))
nodes <- if (file.exists(file.path(data_dir, "nodes.csv"))) {
  read.csv(file.path(data_dir, "nodes.csv"), check.names = FALSE)
} else {
  context$nodes %||% data.frame(node = sort(unique(timeseries$node)))
}

numeric_columns <- names(timeseries)[vapply(timeseries, is.numeric, logical(1))]
measure_columns <- setdiff(numeric_columns, c("node", "timestep"))

default_measures <- intersect(
  c("p_detect_lm_730_3650", "EIR_gamb", "n_inc_clinical_182_5475"),
  measure_columns
)
if (length(default_measures) == 0L) {
  default_measures <- head(measure_columns, 3L)
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

metric_choices <- stats::setNames(
  measure_columns,
  ifelse(
    measure_columns %in% names(metric_labels),
    paste0(metric_labels[measure_columns], " (", measure_columns, ")"),
    measure_columns
  )
)

count_like <- function(variable) {
  grepl("(^n_|_count$|count$|deaths$|total_|^S_count$|^A_count$|^D_count$|^U_count$|^Tr_count$)", variable)
}

to_long <- function(df, measures) {
  pieces <- lapply(measures, function(measure) {
    data.frame(
      node = df$node,
      arm = df$arm,
      timestep = df$timestep,
      variable = measure,
      value = df[[measure]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

aggregate_series <- function(df, measures, method) {
  proportion_specs <- list(
    p_detect_lm_730_3650 = c("n_detect_lm_730_3650", "n_age_730_3650"),
    p_inc_clinical_182_5475 = c("n_inc_clinical_182_5475", "n_age_182_5475")
  )

  aggregate_one <- function(rows, variable) {
    values <- rows[[variable]]
    if (method == "Auto") {
      spec <- proportion_specs[[variable]]
      if (!is.null(spec) && all(spec %in% names(rows))) {
        denominator <- sum(rows[[spec[[2L]]]], na.rm = TRUE)
        if (denominator == 0) {
          return(NA_real_)
        }
        return(sum(rows[[spec[[1L]]]], na.rm = TRUE) / denominator)
      }
      if (count_like(variable)) {
        return(sum(values, na.rm = TRUE))
      }
      return(mean(values, na.rm = TRUE))
    }
    if (method == "Sum") {
      return(sum(values, na.rm = TRUE))
    }
    if (method == "Median") {
      return(stats::median(values, na.rm = TRUE))
    }
    mean(values, na.rm = TRUE)
  }

  df$key <- paste(df$arm, df$timestep, sep = "\r")
  split_rows <- split(df, df$key, drop = TRUE)

  pieces <- lapply(measures, function(measure) {
    values <- vapply(split_rows, aggregate_one, numeric(1), variable = measure)
    keys <- do.call(rbind, strsplit(names(values), "\r", fixed = TRUE))

    data.frame(
      arm = factor(keys[, 1], levels = c("no_release", "release")),
      timestep = as.integer(keys[, 2]),
      variable = measure,
      value = as.numeric(values),
      scope = "All nodes",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, pieces)
}

make_parameter_rows <- function() {
  rows <- list()
  add <- function(category, parameter, value, description) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Category = category,
      Parameter = parameter,
      Value = format_value(value),
      Description = description,
      stringsAsFactors = FALSE
    )
  }

  add("Study design", "Arms", "release, no_release", "Simulation arms available in the time-series CSV files.")
  add("Study design", "Number of nodes", context$n_nodes %||% length(unique(timeseries$node)), "Landscape nodes included in the metapopulation simulation.")
  add("Study design", "Release nodes", context$release_nodes, "Nodes receiving the gene-drive release.")
  add("Study design", "Release day", context$release_day, "Post-warmup day when release mosquitoes are introduced.")
  add("Study design", "Readout day", context$readout_day, "Primary readout day saved by the careful workflow.")
  add("Study design", "Horizon day", context$horizon_day, "Follow-up window after release.")
  add("Study design", "Random seed", context$seed, "Seed used for the paired release and no-release simulations.")

  add("Human attributes", "Total human population", sum(nodes$NH_per_node %||% context$NH_per_node, na.rm = TRUE), "Total people represented across all nodes.")
  add("Human attributes", "Node population", nodes$NH_per_node %||% context$NH_per_node, "Human population assigned to each node.")
  add("Human attributes", "Age structure", context$busia_age_structure %||% calibration$busia_age_structure, "Busia age groups used to set demography.")
  add("Human attributes", "Contact multiplier", context$contact_multiplier, "Node-level relative human blood-meal contact multiplier.")
  add("Human attributes", "Target PfPR", calibration$target_prevalence %||% context$target_prevalence, "Calibration target for parasite prevalence.")
  add("Human attributes", "Realised PfPR", calibration$realised_pfpr, "Realised calibrated prevalence from the careful calibration stage.")

  add("Mosquito biology", "qE", theta$qE, "Daily egg development probability.")
  add("Mosquito biology", "nE", theta$nE, "Number of egg-stage compartments.")
  add("Mosquito biology", "qL", theta$qL, "Daily larval development probability.")
  add("Mosquito biology", "nL", theta$nL, "Number of larval-stage compartments.")
  add("Mosquito biology", "qP", theta$qP, "Daily pupal development probability.")
  add("Mosquito biology", "nP", theta$nP, "Number of pupal-stage compartments.")
  add("Mosquito biology", "muE", theta$muE, "Egg mortality rate.")
  add("Mosquito biology", "muL", theta$muL, "Larval mortality rate.")
  add("Mosquito biology", "muP", theta$muP, "Pupal mortality rate.")
  add("Mosquito biology", "muF", theta$muF, "Adult female mortality rate.")
  add("Mosquito biology", "muM", theta$muM, "Adult male mortality rate.")
  add("Mosquito biology", "beta", theta$beta, "Density-dependent aquatic carrying-capacity parameter.")
  add("Mosquito biology", "nEIP", theta$nEIP, "Extrinsic incubation period compartments.")
  add("Mosquito biology", "nu", theta$nu, "Daily oviposition/blood-feeding cycle parameter.")

  add("Landscape", "Latitude range", nodes$latitude, "North-south coordinate span of the Busia cluster landscape.")
  add("Landscape", "Longitude range", nodes$longitude, "East-west coordinate span of the Busia cluster landscape.")
  add("Landscape", "Movement mu", context$movement_settings$mu, "Mosquito movement scale supplied to the careful workflow.")
  add("Landscape", "Movement probability", context$movement_settings$p_move, "Daily mosquito movement probability supplied to the careful workflow.")
  add("Landscape", "Seasonality g0", context$seasonality$g0, "Mean seasonal forcing level.")
  add("Landscape", "Seasonality g", context$seasonality$g, "Cosine coefficients for seasonal rainfall forcing.")
  add("Landscape", "Seasonality h", context$seasonality$h, "Sine coefficients for seasonal rainfall forcing.")
  add("Landscape", "Rainfall floor", context$seasonality$rainfall_floor, "Minimum rainfall multiplier used by the seasonal forcing function.")

  add("Calibration and warmup", "Initial EIR", context$init_EIR %||% calibration$init_EIR, "Calibrated initial EIR used to initialize endemic transmission.")
  add("Calibration and warmup", "Burn-in timesteps", context$burnin_timesteps, "Warmup duration used before the release/no-release comparison.")
  add("Calibration and warmup", "Snapshots", context$n_snapshots, "Number of baseline checkpoint candidates retained.")
  add("Calibration and warmup", "Promoted snapshot index", context$promoted_index, "Checkpoint selected for the careful release comparison.")
  add("Calibration and warmup", "Multi-seed pass", context$multi_seed_pass, "Whether the promoted checkpoint passed the multi-seed validation check.")

  do.call(rbind, rows)
}

parameter_rows <- make_parameter_rows()

brewer_palettes <- rownames(RColorBrewer::brewer.pal.info[
  RColorBrewer::brewer.pal.info$category %in% c("qual", "seq", "div"),
])

ui <- navbarPage(
  title = div(class = "brand", "Busia Careful Output"),
  id = "page",
  header = tags$head(tags$style(HTML("
    :root {
      --bg: #f6f7fb;
      --panel: #ffffff;
      --ink: #20242c;
      --muted: #6b7280;
      --line: #d9dee8;
      --accent: #2f6f73;
      --accent-soft: #e4f1ef;
    }
    body {
      background: var(--bg);
      color: var(--ink);
      font-family: Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    }
    .navbar {
      background: var(--panel);
      border: 0;
      border-bottom: 1px solid var(--line);
      box-shadow: 0 8px 30px rgba(32, 36, 44, 0.05);
    }
    .navbar-default .navbar-brand, .navbar-default .navbar-nav > li > a {
      color: var(--ink);
    }
    .navbar-default .navbar-nav > .active > a,
    .navbar-default .navbar-nav > .active > a:focus,
    .navbar-default .navbar-nav > .active > a:hover {
      background: var(--accent-soft);
      color: var(--accent);
    }
    .brand {
      font-weight: 700;
      letter-spacing: 0;
    }
    .page-shell {
      max-width: 1500px;
      margin: 20px auto 36px auto;
      padding: 0 18px;
    }
    .control-panel, .plot-panel, .summary-card {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      box-shadow: 0 10px 28px rgba(32, 36, 44, 0.06);
    }
    .control-panel {
      padding: 18px;
    }
    .plot-panel {
      padding: 16px 18px 10px 18px;
    }
    .summary-grid {
      display: grid;
      grid-template-columns: repeat(4, minmax(160px, 1fr));
      gap: 12px;
      margin-bottom: 16px;
    }
    .summary-card {
      padding: 14px 16px;
    }
    .summary-label {
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: .06em;
      margin-bottom: 4px;
    }
    .summary-value {
      font-size: 22px;
      font-weight: 700;
    }
    label {
      color: #333946;
      font-weight: 650;
      font-size: 13px;
    }
    .selectize-input, .form-control {
      border-radius: 7px;
      border-color: #cfd6e3;
      box-shadow: none;
    }
    .selectize-input.focus, .form-control:focus {
      border-color: var(--accent);
      box-shadow: 0 0 0 3px rgba(47, 111, 115, .12);
    }
    .btn-default, .btn {
      border-radius: 7px;
    }
    .help-text {
      color: var(--muted);
      font-size: 12px;
      line-height: 1.45;
      margin-top: -4px;
      margin-bottom: 14px;
    }
    .section-title {
      font-size: 18px;
      font-weight: 750;
      margin: 0 0 10px 0;
    }
    .section-subtitle {
      color: var(--muted);
      margin-bottom: 18px;
    }
    .dataTables_wrapper {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 12px;
      box-shadow: 0 10px 28px rgba(32, 36, 44, 0.05);
    }
    @media (max-width: 900px) {
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
  tabPanel(
    "Timeseries",
    div(
      class = "page-shell",
      div(
        class = "summary-grid",
        div(class = "summary-card", div(class = "summary-label", "Nodes"), div(class = "summary-value", textOutput("node_count", inline = TRUE))),
        div(class = "summary-card", div(class = "summary-label", "Timesteps"), div(class = "summary-value", textOutput("time_count", inline = TRUE))),
        div(class = "summary-card", div(class = "summary-label", "Metrics"), div(class = "summary-value", textOutput("metric_count", inline = TRUE))),
        div(class = "summary-card", div(class = "summary-label", "Rows"), div(class = "summary-value", textOutput("row_count", inline = TRUE)))
      ),
      fluidRow(
        column(
          width = 3,
          div(
            class = "control-panel",
            div(class = "section-title", "Plot controls"),
            checkboxGroupInput(
              "arms",
              "Simulation arms",
              choices = c("No release" = "no_release", "Release" = "release"),
              selected = c("no_release", "release")
            ),
            radioButtons(
              "scope",
              "Node scope",
              choices = c("Aggregate all nodes" = "aggregate", "Individual node data" = "node"),
              selected = "aggregate"
            ),
            conditionalPanel(
              "input.scope == 'node'",
              selectizeInput(
                "nodes",
                "Node(s)",
                choices = sort(unique(timeseries$node)),
                selected = head(sort(unique(timeseries$node)), 1L),
                multiple = TRUE,
                options = list(plugins = list("remove_button"), maxItems = 12)
              )
            ),
            conditionalPanel(
              "input.scope == 'aggregate'",
              selectInput(
                "aggregate_method",
                "Aggregation",
                choices = c("Auto", "Mean", "Sum", "Median"),
                selected = "Auto"
              ),
              div(class = "help-text", "Auto sums count-like metrics and averages rates or mean-valued metrics.")
            ),
            selectizeInput(
              "measures",
              "Columns to plot",
              choices = metric_choices,
              selected = default_measures,
              multiple = TRUE,
              options = list(plugins = list("remove_button"))
            ),
            fluidRow(
              column(
                6,
                numericInput(
                  "time_min",
                  "Start day",
                  value = min(timeseries$timestep, na.rm = TRUE),
                  min = min(timeseries$timestep, na.rm = TRUE),
                  max = max(timeseries$timestep, na.rm = TRUE),
                  step = 1
                )
              ),
              column(
                6,
                numericInput(
                  "time_max",
                  "End day",
                  value = max(timeseries$timestep, na.rm = TRUE),
                  min = min(timeseries$timestep, na.rm = TRUE),
                  max = max(timeseries$timestep, na.rm = TRUE),
                  step = 1
                )
              )
            ),
            selectInput(
              "palette",
              "ColorBrewer palette",
              choices = brewer_palettes,
              selected = "Dark2"
            ),
            checkboxInput("facet", "Separate panels by column", value = length(default_measures) > 1L),
            checkboxInput("points", "Show points", value = FALSE),
            sliderInput("line_width", "Line width", min = 0.3, max = 2, value = 0.9, step = 0.1)
          )
        ),
        column(
          width = 9,
          div(
            class = "plot-panel",
            div(class = "section-title", "Timeseries plot"),
            plotOutput("timeseries_plot", height = "650px"),
            downloadButton("download_plot", "Download plot")
          )
        )
      )
    )
  ),
  tabPanel(
    "Parameters",
    div(
      class = "page-shell",
      div(class = "section-title", "Run parameters"),
      div(class = "section-subtitle", "Values are read from the careful workflow outputs beside the CSV files."),
      DT::DTOutput("parameter_table")
    )
  ),
  tabPanel(
    "Data Preview",
    div(
      class = "page-shell",
      div(class = "section-title", "Filtered data"),
      div(class = "section-subtitle", "Preview follows the same arm, scope, node, time, and column selections as the plot."),
      DT::DTOutput("preview_table")
    )
  )
)

server <- function(input, output, session) {
  output$node_count <- renderText(scales::comma(length(unique(timeseries$node))))
  output$time_count <- renderText(scales::comma(length(unique(timeseries$timestep))))
  output$metric_count <- renderText(scales::comma(length(measure_columns)))
  output$row_count <- renderText(scales::comma(nrow(timeseries)))

  observeEvent(input$scope, {
    if (identical(input$scope, "node") && length(input$nodes) == 0L) {
      updateSelectizeInput(
        session,
        "nodes",
        selected = head(sort(unique(timeseries$node)), 1L)
      )
    }
  })

  selected_data <- reactive({
    req(input$arms, input$measures)
    start_day <- min(input$time_min, input$time_max, na.rm = TRUE)
    end_day <- max(input$time_min, input$time_max, na.rm = TRUE)
    measures <- intersect(input$measures, measure_columns)
    req(length(measures) > 0L)

    df <- timeseries[
      as.character(timeseries$arm) %in% input$arms &
        timeseries$timestep >= start_day &
        timeseries$timestep <= end_day,
      ,
      drop = FALSE
    ]

    if (identical(input$scope, "node")) {
      req(input$nodes)
      df <- df[df$node %in% as.integer(input$nodes), , drop = FALSE]
      long <- to_long(df, measures)
      long$scope <- paste0("Node ", long$node)
      return(long)
    }

    aggregate_series(df, measures, input$aggregate_method %||% "Auto")
  })

  plot_data <- reactive({
    df <- selected_data()
    req(nrow(df) > 0L)
    df$variable_label <- ifelse(
      df$variable %in% names(metric_labels),
      metric_labels[df$variable],
      df$variable
    )
    df$series <- if (length(unique(df$variable)) == 1L) {
      paste(df$scope, df$arm, sep = " | ")
    } else if (identical(input$scope, "node")) {
      paste(df$variable_label, df$scope, df$arm, sep = " | ")
    } else {
      paste(df$variable_label, df$arm, sep = " | ")
    }
    df
  })

  build_plot <- reactive({
    df <- plot_data()
    palette <- input$palette %||% "Dark2"
    palette_info <- RColorBrewer::brewer.pal.info[palette, ]
    max_colors <- palette_info$maxcolors
    series_count <- length(unique(df$series))
    colors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(min(max_colors, max(3L, min(series_count, max_colors))), palette)
    )(max(3L, series_count))

    p <- ggplot(
      df,
      aes(
        x = timestep,
        y = value,
        color = series,
        group = interaction(series, variable, scope, arm)
      )
    ) +
      geom_line(linewidth = input$line_width %||% 0.9, alpha = 0.95) +
      scale_color_manual(values = colors) +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        x = "Timestep (days)",
        y = NULL,
        color = NULL,
        title = if (identical(input$scope, "aggregate")) {
          "Aggregate time-series across all nodes"
        } else {
          "Individual node time-series"
        },
        subtitle = paste(unique(as.character(df$arm)), collapse = " vs ")
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 18, color = "#20242c"),
        plot.subtitle = element_text(color = "#6b7280", margin = margin(b = 12)),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "#edf0f5"),
        panel.grid.major.y = element_line(color = "#e4e8f0"),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "bottom",
        legend.key.width = unit(18, "pt"),
        legend.text = element_text(size = 10),
        strip.text = element_text(face = "bold", color = "#20242c"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )

    if (isTRUE(input$points)) {
      p <- p + geom_point(size = 1.1, alpha = 0.65)
    }

    if (isTRUE(input$facet) && length(unique(df$variable_label)) > 1L) {
      p <- p + facet_wrap(~ variable_label, scales = "free_y", ncol = 1)
    }

    p
  })

  output$timeseries_plot <- renderPlot({
    build_plot()
  }, res = 120)

  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("busia_timeseries_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggplot2::ggsave(
        filename = file,
        plot = build_plot(),
        width = 12,
        height = if (isTRUE(input$facet)) 8 else 6.5,
        dpi = 180
      )
    }
  )

  output$parameter_table <- DT::renderDT({
    DT::datatable(
      parameter_rows,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 20,
        autoWidth = TRUE,
        order = list(list(0, "asc"))
      )
    )
  })

  output$preview_table <- DT::renderDT({
    df <- plot_data()
    df <- df[order(df$arm, df$variable, df$scope, df$timestep), , drop = FALSE]
    df$metric <- df$variable_label
    df <- df[, c("arm", "timestep", "scope", "metric", "variable", "value"), drop = FALSE]
    DT::datatable(
      df,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })
}

shinyApp(ui, server)
