#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 01_compare_legacy_native.R
# ------------------------------------------------------------------------------
# Compare one-node no-release baseline outputs from the legacy mosquito model
# and the native deterministic mosquito backend.
#
# Run from the package root:
#   Rscript test-scripts/one_node_example/baseline/01_compare_legacy_native.R
#
# Optional environment overrides:
#   MSIMGD_COMPARE_SEEDS=1,2,3
#   MSIMGD_COMPARE_TIMESTEPS=1095
#   MSIMGD_COMPARE_INIT_EIR=10
#   MSIMGD_COMPARE_HUMAN_POPULATION=1000
#   MSIMGD_COMPARE_NATIVE_NE=2
#   MSIMGD_COMPARE_NATIVE_NL=3
#   MSIMGD_COMPARE_NATIVE_NP=2
#   MSIMGD_COMPARE_NATIVE_NEIP=50
#   MSIMGD_COMPARE_NATIVE_NU=1
#   MSIMGD_COMPARE_NATIVE_TAU_STEP=0.1
# ------------------------------------------------------------------------------

start_time <- Sys.time()

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/one_node_example/baseline/01_compare_legacy_native.R",
    mustWork = TRUE
  )
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
pkg_root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)
out_dir <- file.path(example_dir, "output", "baseline", "legacy_native_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_local_package <- function(pkg_root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(pkg_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (requireNamespace("malariasimulationGD", quietly = TRUE)) {
    library(malariasimulationGD)
    return(invisible(TRUE))
  }
  stop(
    "Install malariasimulationGD or pkgload first. From the package root, try: Rscript test-scripts/install_local.R",
    call. = FALSE
  )
}

parse_seed_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  seeds <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (anyNA(seeds)) {
    stop(sprintf("%s must be a comma-separated list of integers.", name), call. = FALSE)
  }
  seeds
}

parse_numeric_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  parsed <- as.numeric(value)
  if (length(parsed) != 1L || !is.finite(parsed)) {
    stop(sprintf("%s must be a finite number.", name), call. = FALSE)
  }
  parsed
}

base_overrides <- function(human_population) {
  list(
    human_population = human_population,
    individual_mosquitoes = FALSE,
    model_seasonality = FALSE,
    human_mobility_enabled = FALSE,
    human_move_probs = NULL,
    human_move_rates = NULL,
    move_probs = matrix(1, 1, 1),
    move_rates = 0,
    bednets = FALSE,
    spraying = FALSE,
    progress_bar = FALSE,
    cube = NULL
  )
}

add_rendering_outputs <- function(parameters) {
  parameters$incidence_rendering_min_ages <- 0
  parameters$incidence_rendering_max_ages <- 100 * 365 - 1
  parameters
}

build_parameters <- function(backend, human_population, init_EIR,
                             native_stage_config = list()) {
  backend <- match.arg(backend, c("legacy", "native"))
  overrides <- base_overrides(human_population)
  overrides$native_mosquito_backend <- identical(backend, "native")

  parameters <- malariasimulationGD::get_parameters(overrides)
  parameters <- add_rendering_outputs(parameters)
  if (identical(backend, "native")) {
    parameters <- apply_native_stage_overrides(parameters, native_stage_config)
  }
  malariasimulationGD::set_equilibrium(
    parameters,
    init_EIR = init_EIR,
    native_total_M = identical(backend, "native")
  )
}

apply_native_stage_overrides <- function(parameters, config) {
  if (!is.null(config$nE)) {
    parameters$native_mosquito_nE <- as.integer(config$nE)
  }
  if (!is.null(config$nL)) {
    parameters$native_mosquito_nL <- as.integer(config$nL)
  }
  if (!is.null(config$nP)) {
    parameters$native_mosquito_nP <- as.integer(config$nP)
  }
  if (!is.null(config$nEIP)) {
    parameters$native_mosquito_nEIP <- as.integer(config$nEIP)
  }
  if (!is.null(config$nu)) {
    parameters$native_mosquito_nu <- as.numeric(config$nu)
  }
  if (!is.null(config$tau_step)) {
    parameters$mosquito_tau_step <- as.numeric(config$tau_step)
  }
  parameters
}

native_stage_config_from_env <- function() {
  value_or_null <- function(name) {
    value <- Sys.getenv(name, unset = "")
    if (!nzchar(value)) {
      return(NULL)
    }
    parsed <- as.numeric(value)
    if (length(parsed) != 1L || !is.finite(parsed)) {
      stop(sprintf("%s must be a finite number.", name), call. = FALSE)
    }
    parsed
  }
  list(
    nE = value_or_null("MSIMGD_COMPARE_NATIVE_NE"),
    nL = value_or_null("MSIMGD_COMPARE_NATIVE_NL"),
    nP = value_or_null("MSIMGD_COMPARE_NATIVE_NP"),
    nEIP = value_or_null("MSIMGD_COMPARE_NATIVE_NEIP"),
    nu = value_or_null("MSIMGD_COMPARE_NATIVE_NU"),
    tau_step = value_or_null("MSIMGD_COMPARE_NATIVE_TAU_STEP")
  )
}

safe_divide <- function(numerator, denominator) {
  numerator / pmax(denominator, 1)
}

adult_female_total <- function(data) {
  data$Sm_gamb_count + data$Pm_gamb_count + data$Im_gamb_count
}

prepare_timeseries <- function(data, backend, seed, parameters) {
  data$backend <- backend
  data$seed <- seed
  data$target_init_EIR <- parameters$init_EIR
  data$parameter_total_M <- parameters$total_M
  data$adult_female_total <- adult_female_total(data)
  data$pfpr_2_10_lm <- safe_divide(data$n_detect_lm_730_3650, data$n_age_730_3650)
  data$annualized_daily_eir_per_person <- 365 * data$EIR_gamb / parameters$human_population
  data
}

summarise_run <- function(data, backend, seed, parameters) {
  final_year <- utils::tail(data, min(365L, nrow(data)))
  human_total <- data$S_count + data$A_count + data$D_count + data$U_count + data$Tr_count
  key_columns <- c(
    "E_gamb_count", "L_gamb_count", "P_gamb_count",
    "Sm_gamb_count", "Pm_gamb_count", "Im_gamb_count",
    "total_M_gamb", "EIR_gamb", "FOIM_gamb", "mu_gamb",
    "n_age_730_3650", "n_detect_lm_730_3650", "n_inc_0_36499"
  )
  key_values <- as.matrix(data[intersect(key_columns, names(data))])

  data.frame(
    seed = seed,
    backend = backend,
    timesteps = nrow(data),
    human_population = parameters$human_population,
    target_init_EIR = parameters$init_EIR,
    parameter_total_M = parameters$total_M,
    mean_E_final_year = mean(final_year$E_gamb_count),
    mean_L_final_year = mean(final_year$L_gamb_count),
    mean_P_final_year = mean(final_year$P_gamb_count),
    mean_Sm_final_year = mean(final_year$Sm_gamb_count),
    mean_Pm_final_year = mean(final_year$Pm_gamb_count),
    mean_Im_final_year = mean(final_year$Im_gamb_count),
    mean_adult_female_final_year = mean(final_year$adult_female_total),
    cv_adult_female_final_year = stats::sd(final_year$adult_female_total) /
      mean(final_year$adult_female_total),
    annual_eir_per_person_final_year = sum(final_year$EIR_gamb) /
      parameters$human_population,
    mean_pfpr_2_10_lm_final_year = mean(final_year$pfpr_2_10_lm),
    annual_infection_incidence_final_year = sum(final_year$n_inc_0_36499) /
      parameters$human_population,
    final_pfpr_2_10_lm = utils::tail(data$pfpr_2_10_lm, 1),
    final_E = utils::tail(data$E_gamb_count, 1),
    final_L = utils::tail(data$L_gamb_count, 1),
    final_P = utils::tail(data$P_gamb_count, 1),
    final_Sm = utils::tail(data$Sm_gamb_count, 1),
    final_Pm = utils::tail(data$Pm_gamb_count, 1),
    final_Im = utils::tail(data$Im_gamb_count, 1),
    human_total_min = min(human_total),
    human_total_max = max(human_total),
    key_outputs_finite_nonnegative = all(is.finite(key_values)) && all(key_values >= 0),
    stringsAsFactors = FALSE
  )
}

run_one <- function(seed, backend, parameters, timesteps) {
  set.seed(seed)
  sim <- malariasimulationGD::run_simulation(
    timesteps = timesteps,
    parameters = parameters
  )
  sim <- prepare_timeseries(sim, backend, seed, parameters)
  list(
    timeseries = sim,
    summary = summarise_run(sim, backend, seed, parameters)
  )
}

compare_summaries <- function(summary_by_backend) {
  metric_names <- setdiff(
    names(summary_by_backend),
    c("seed", "backend", "timesteps", "human_population", "target_init_EIR",
      "key_outputs_finite_nonnegative")
  )

  rows <- vector("list", length(unique(summary_by_backend$seed)) * length(metric_names))
  row_i <- 0L
  for (seed in sort(unique(summary_by_backend$seed))) {
    legacy <- summary_by_backend[summary_by_backend$seed == seed &
                                   summary_by_backend$backend == "legacy", , drop = FALSE]
    native <- summary_by_backend[summary_by_backend$seed == seed &
                                   summary_by_backend$backend == "native", , drop = FALSE]
    if (nrow(legacy) != 1L || nrow(native) != 1L) {
      next
    }
    for (metric in metric_names) {
      legacy_value <- legacy[[metric]]
      native_value <- native[[metric]]
      if (!is.numeric(legacy_value) || !is.numeric(native_value)) {
        next
      }
      row_i <- row_i + 1L
      absolute_difference <- native_value - legacy_value
      relative_difference <- if (is.finite(legacy_value) && abs(legacy_value) > 0) {
        absolute_difference / legacy_value
      } else {
        NA_real_
      }
      rows[[row_i]] <- data.frame(
        seed = seed,
        metric = metric,
        legacy = legacy_value,
        native = native_value,
        native_minus_legacy = absolute_difference,
        relative_difference = relative_difference,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows[seq_len(row_i)])
}

write_optional_plots <- function(timeseries, out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' is not available; skipping optional PNG plots.")
    return(invisible(FALSE))
  }

  library(ggplot2)
  timeseries$backend <- factor(timeseries$backend, levels = c("legacy", "native"))
  timeseries$seed <- factor(timeseries$seed)
  base_theme <- theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title.position = "plot")

  write_plot <- function(filename, plot) {
    ggplot2::ggsave(
      file.path(out_dir, filename),
      plot,
      width = 8,
      height = 4.5,
      units = "in",
      dpi = 150
    )
  }

  write_plot(
    "eir_compare.png",
    ggplot(timeseries, aes(timestep, annualized_daily_eir_per_person, color = backend)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~ seed) +
      labs(title = "Annualized Daily EIR", x = "Day", y = "Annualized EIR/person") +
      base_theme
  )
  write_plot(
    "pfpr_2_10_compare.png",
    ggplot(timeseries, aes(timestep, pfpr_2_10_lm, color = backend)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~ seed) +
      labs(title = "PfPR 2-10 by Microscopy", x = "Day", y = "PfPR 2-10") +
      base_theme
  )
  write_plot(
    "adult_female_total_compare.png",
    ggplot(timeseries, aes(timestep, adult_female_total, color = backend)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~ seed, scales = "free_y") +
      labs(title = "Adult Female Mosquitoes", x = "Day", y = "Adult females") +
      base_theme
  )

  compartments <- rbind(
    data.frame(seed = timeseries$seed, backend = timeseries$backend,
               timestep = timeseries$timestep, compartment = "E",
               count = timeseries$E_gamb_count),
    data.frame(seed = timeseries$seed, backend = timeseries$backend,
               timestep = timeseries$timestep, compartment = "L",
               count = timeseries$L_gamb_count),
    data.frame(seed = timeseries$seed, backend = timeseries$backend,
               timestep = timeseries$timestep, compartment = "P",
               count = timeseries$P_gamb_count),
    data.frame(seed = timeseries$seed, backend = timeseries$backend,
               timestep = timeseries$timestep, compartment = "Sm",
               count = timeseries$Sm_gamb_count),
    data.frame(seed = timeseries$seed, backend = timeseries$backend,
               timestep = timeseries$timestep, compartment = "Pm",
               count = timeseries$Pm_gamb_count),
    data.frame(seed = timeseries$seed, backend = timeseries$backend,
               timestep = timeseries$timestep, compartment = "Im",
               count = timeseries$Im_gamb_count)
  )
  write_plot(
    "entomological_compartments_compare.png",
    ggplot(compartments, aes(timestep, count, color = backend)) +
      geom_line(linewidth = 0.55) +
      facet_grid(compartment ~ seed, scales = "free_y") +
      labs(title = "Entomological Compartments", x = "Day", y = "Count") +
      base_theme
  )

  invisible(TRUE)
}

load_local_package(pkg_root)

seeds <- parse_seed_env("MSIMGD_COMPARE_SEEDS", 1:3)
timesteps <- as.integer(parse_numeric_env("MSIMGD_COMPARE_TIMESTEPS", 3 * 365))
init_EIR <- parse_numeric_env("MSIMGD_COMPARE_INIT_EIR", 10)
human_population <- as.integer(parse_numeric_env("MSIMGD_COMPARE_HUMAN_POPULATION", 1000))

if (timesteps < 365L) {
  stop("MSIMGD_COMPARE_TIMESTEPS must be at least 365.", call. = FALSE)
}
if (human_population <= 0L) {
  stop("MSIMGD_COMPARE_HUMAN_POPULATION must be positive.", call. = FALSE)
}

cat(sprintf("[%s] legacy vs native no-release comparison\n", format(Sys.time(), "%F %T")))
cat(sprintf("  package root: %s\n", pkg_root))
cat(sprintf("  output dir:   %s\n", out_dir))
cat(sprintf("  seeds:        %s\n", paste(seeds, collapse = ", ")))
cat(sprintf("  timesteps:    %d\n", timesteps))
cat(sprintf("  human pop:    %d\n", human_population))
cat(sprintf("  init EIR:     %.2f annual infectious bites/person/year\n", init_EIR))

native_stage_config <- native_stage_config_from_env()
parameters <- list(
  legacy = build_parameters("legacy", human_population, init_EIR),
  native = build_parameters(
    "native",
    human_population,
    init_EIR,
    native_stage_config = native_stage_config
  )
)

runs <- list()
run_i <- 0L
for (seed in seeds) {
  for (backend in c("legacy", "native")) {
    run_i <- run_i + 1L
    runs[[run_i]] <- run_one(
      seed = seed,
      backend = backend,
      parameters = parameters[[backend]],
      timesteps = timesteps
    )
  }
}

timeseries <- do.call(rbind, lapply(runs, `[[`, "timeseries"))
summary_by_backend <- do.call(rbind, lapply(runs, `[[`, "summary"))
comparison <- compare_summaries(summary_by_backend)

context <- list(
  scenario = "one_node_no_release_legacy_vs_native",
  package_root = pkg_root,
  seeds = seeds,
  timesteps = timesteps,
  init_EIR = init_EIR,
  human_population = human_population,
  legacy_total_M = parameters$legacy$total_M,
  native_total_M = parameters$native$total_M,
  native_mosquito_nE = parameters$native$native_mosquito_nE,
  native_mosquito_nL = parameters$native$native_mosquito_nL,
  native_mosquito_nP = parameters$native$native_mosquito_nP,
  native_mosquito_nEIP = parameters$native$native_mosquito_nEIP,
  native_mosquito_nu = parameters$native$native_mosquito_nu,
  mosquito_tau_step = parameters$native$mosquito_tau_step,
  legacy_native_backend = isTRUE(parameters$legacy$native_mosquito_backend),
  native_native_backend = isTRUE(parameters$native$native_mosquito_backend),
  no_releases = TRUE,
  no_landscape = TRUE,
  no_mosquito_movement = TRUE,
  no_human_movement = TRUE,
  no_seasonality = TRUE,
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
)

utils::write.csv(timeseries, file.path(out_dir, "timeseries_by_backend.csv"), row.names = FALSE)
utils::write.csv(summary_by_backend, file.path(out_dir, "summary_by_backend.csv"), row.names = FALSE)
utils::write.csv(comparison, file.path(out_dir, "native_minus_legacy_comparison.csv"), row.names = FALSE)
saveRDS(context, file.path(out_dir, "context.rds"))
write_optional_plots(timeseries, out_dir)

cat("Saved outputs:\n")
cat(sprintf("  %s\n", file.path(out_dir, "timeseries_by_backend.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "summary_by_backend.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "native_minus_legacy_comparison.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "context.rds")))

cat("Final-year summaries by backend:\n")
print(
  summary_by_backend[
    ,
    c(
      "seed", "backend", "parameter_total_M",
      "annual_eir_per_person_final_year",
      "mean_pfpr_2_10_lm_final_year",
      "mean_adult_female_final_year",
      "mean_E_final_year", "mean_L_final_year", "mean_P_final_year",
      "mean_Sm_final_year", "mean_Pm_final_year", "mean_Im_final_year"
    )
  ],
  row.names = FALSE
)

if (!all(summary_by_backend$key_outputs_finite_nonnegative)) {
  failed <- summary_by_backend[!summary_by_backend$key_outputs_finite_nonnegative, ]
  print(failed[, c("seed", "backend")], row.names = FALSE)
  stop("Some backend runs produced negative or non-finite key outputs.", call. = FALSE)
}

if (any(summary_by_backend$human_total_min != human_population) ||
    any(summary_by_backend$human_total_max != human_population)) {
  stop("At least one backend run failed human population accounting.", call. = FALSE)
}

cat("Comparison checks passed: key outputs are finite/non-negative and human totals are conserved.\n")
