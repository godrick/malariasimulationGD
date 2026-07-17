#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 02_compare_eir_pfpr_mapping.R
# ------------------------------------------------------------------------------
# Compare the configured EIR -> realised PfPR2-10 mapping for the legacy and
# native mosquito backends in the one-node, no-release baseline. The native
# model is given a singleton WW cube, so only wild-type mosquitoes can exist.
# The legacy model is the aggregate wild-type baseline and does not track
# mosquito genotypes.
#
# For every configured EIR, backend, and seed, the script reports:
#   * configured annual EIR (the value passed to set_equilibrium());
#   * realised annual EIR per person over the final measurement window; and
#   * realised microscopy PfPR2-10 averaged over the same window.
#
# Run from the package root:
#   Rscript test-scripts/one_node_example/baseline/02_compare_eir_pfpr_mapping.R
#
# Optional environment overrides:
#   MSIMGD_MAPPING_EIRS=0.1,0.5,1,2,5,10,20,40
#   MSIMGD_MAPPING_SEEDS=1,2,3
#   MSIMGD_MAPPING_TIMESTEPS=1095
#   MSIMGD_MAPPING_MEASUREMENT_DAYS=365
#   MSIMGD_MAPPING_HUMAN_POPULATION=1000
#   MSIMGD_MAPPING_NATIVE_NE=2
#   MSIMGD_MAPPING_NATIVE_NL=3
#   MSIMGD_MAPPING_NATIVE_NP=2
#   MSIMGD_MAPPING_NATIVE_NEIP=50
#   MSIMGD_MAPPING_NATIVE_NU=1
#   MSIMGD_MAPPING_NATIVE_TAU_STEP=0.1
# ------------------------------------------------------------------------------

start_time <- Sys.time()

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/one_node_example/baseline/02_compare_eir_pfpr_mapping.R",
    mustWork = TRUE
  )
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
pkg_root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)
out_dir <- file.path(example_dir, "output", "baseline", "eir_pfpr_mapping")
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
    paste(
      "Install malariasimulationGD or pkgload first.",
      "From the package root, try: Rscript test-scripts/install_local.R"
    ),
    call. = FALSE
  )
}

parse_numeric_vector_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  parsed <- as.numeric(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (length(parsed) == 0L || any(!is.finite(parsed))) {
    stop(sprintf("%s must be a comma-separated list of finite numbers.", name),
         call. = FALSE)
  }
  parsed
}

parse_integer_vector_env <- function(name, default) {
  parsed <- parse_numeric_vector_env(name, default)
  if (any(parsed != as.integer(parsed))) {
    stop(sprintf("%s must be a comma-separated list of integers.", name),
         call. = FALSE)
  }
  as.integer(parsed)
}

parse_numeric_env <- function(name, default) {
  parsed <- parse_numeric_vector_env(name, default)
  if (length(parsed) != 1L) {
    stop(sprintf("%s must contain one number.", name), call. = FALSE)
  }
  parsed
}

make_wildtype_cube <- function() {
  genotype <- "WW"
  list(
    ih = array(
      1,
      dim = c(1L, 1L, 1L),
      dimnames = list(genotype, genotype, genotype)
    ),
    tau = array(
      1,
      dim = c(1L, 1L, 1L),
      dimnames = list(genotype, genotype, genotype)
    ),
    eta = matrix(
      1,
      nrow = 1L,
      ncol = 1L,
      dimnames = list(genotype, genotype)
    ),
    b = stats::setNames(1, genotype),
    c = stats::setNames(1, genotype),
    phi = stats::setNames(0.5, genotype),
    omega = stats::setNames(1, genotype),
    xiF = stats::setNames(1, genotype),
    xiM = stats::setNames(1, genotype),
    s = stats::setNames(1, genotype),
    genotypesID = genotype,
    wildType = genotype
  )
}

base_overrides <- function(backend, human_population) {
  backend <- match.arg(backend, c("legacy", "native"))
  list(
    human_population = human_population,
    native_mosquito_backend = identical(backend, "native"),
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
    prevalence_rendering_min_ages = 2 * 365,
    prevalence_rendering_max_ages = 10 * 365,
    # The legacy compartmental baseline is implicitly wild-type only. The
    # singleton cube makes that restriction explicit and verifiable natively.
    cube = if (identical(backend, "native")) make_wildtype_cube() else NULL
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
    parse_numeric_env(name, NA_real_)
  }
  list(
    nE = value_or_null("MSIMGD_MAPPING_NATIVE_NE"),
    nL = value_or_null("MSIMGD_MAPPING_NATIVE_NL"),
    nP = value_or_null("MSIMGD_MAPPING_NATIVE_NP"),
    nEIP = value_or_null("MSIMGD_MAPPING_NATIVE_NEIP"),
    nu = value_or_null("MSIMGD_MAPPING_NATIVE_NU"),
    tau_step = value_or_null("MSIMGD_MAPPING_NATIVE_TAU_STEP")
  )
}

build_parameters <- function(backend, human_population, init_eir,
                             native_stage_config) {
  parameters <- malariasimulationGD::get_parameters(
    base_overrides(backend, human_population)
  )
  if (identical(backend, "native")) {
    parameters <- apply_native_stage_overrides(parameters, native_stage_config)
  }
  malariasimulationGD::set_equilibrium(
    parameters,
    init_EIR = init_eir,
    native_total_M = identical(backend, "native")
  )
}

summarise_run <- function(run, backend, seed, init_eir, parameters,
                          measurement_days) {
  data <- run$data
  measurement <- utils::tail(data, measurement_days)
  pfpr_daily <- measurement$n_detect_lm_730_3650 /
    pmax(measurement$n_age_730_3650, 1)

  native_wildtype_verified <- NA
  if (identical(backend, "native")) {
    female <- run$mosquito_genotypes$female
    male <- run$mosquito_genotypes$male
    native_wildtype_verified <-
      identical(colnames(female), "WW") &&
      identical(colnames(male), "WW") &&
      all(is.finite(female[, "WW"])) &&
      all(is.finite(male[, "WW"])) &&
      all(female[, "WW"] >= 0) &&
      all(male[, "WW"] >= 0)
  }

  data.frame(
    seed = seed,
    backend = backend,
    init_eir = init_eir,
    realised_eir = sum(measurement$EIR_gamb) / parameters$human_population,
    realised_pfpr_2_10 = mean(pfpr_daily),
    pfpr_2_10_sd_daily = stats::sd(pfpr_daily),
    pfpr_2_10_min_daily = min(pfpr_daily),
    pfpr_2_10_max_daily = max(pfpr_daily),
    parameter_total_M = parameters$total_M,
    measurement_days = nrow(measurement),
    native_wildtype_verified = native_wildtype_verified,
    stringsAsFactors = FALSE
  )
}

run_one <- function(backend, seed, init_eir, human_population, timesteps,
                    measurement_days, native_stage_config) {
  parameters <- build_parameters(
    backend = backend,
    human_population = human_population,
    init_eir = init_eir,
    native_stage_config = native_stage_config
  )
  set.seed(seed)
  run <- malariasimulationGD::run_resumable_simulation(
    timesteps = timesteps,
    parameters = parameters
  )
  summarise_run(
    run = run,
    backend = backend,
    seed = seed,
    init_eir = init_eir,
    parameters = parameters,
    measurement_days = measurement_days
  )
}

summarise_mapping <- function(mapping_by_seed) {
  groups <- split(
    mapping_by_seed,
    interaction(mapping_by_seed$backend, mapping_by_seed$init_eir, drop = TRUE)
  )
  rows <- lapply(groups, function(x) {
    data.frame(
      backend = x$backend[[1L]],
      init_eir = x$init_eir[[1L]],
      n_seeds = nrow(x),
      mean_realised_eir = mean(x$realised_eir),
      sd_realised_eir = stats::sd(x$realised_eir),
      mean_realised_pfpr_2_10 = mean(x$realised_pfpr_2_10),
      sd_realised_pfpr_2_10 = stats::sd(x$realised_pfpr_2_10),
      min_realised_pfpr_2_10 = min(x$realised_pfpr_2_10),
      max_realised_pfpr_2_10 = max(x$realised_pfpr_2_10),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(result$backend, result$init_eir), ]
}

write_mapping_plot <- function(mapping_by_seed, mapping_summary, out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' is not available; skipping the PNG plot.")
    return(invisible(FALSE))
  }

  mapping_by_seed$backend <- factor(
    mapping_by_seed$backend,
    levels = c("legacy", "native"),
    labels = c("Legacy", "Native (WW only)")
  )
  mapping_summary$backend <- factor(
    mapping_summary$backend,
    levels = c("legacy", "native"),
    labels = c("Legacy", "Native (WW only)")
  )

  plot <- ggplot2::ggplot(
    mapping_by_seed,
    ggplot2::aes(x = init_eir, y = realised_pfpr_2_10, colour = backend)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = interaction(backend, seed)),
      linewidth = 0.45,
      alpha = 0.25
    ) +
    ggplot2::geom_point(alpha = 0.4, size = 1.6) +
    ggplot2::geom_ribbon(
      data = mapping_summary,
      ggplot2::aes(
        x = init_eir,
        y = mean_realised_pfpr_2_10,
        ymin = min_realised_pfpr_2_10,
        ymax = max_realised_pfpr_2_10,
        fill = backend,
        group = backend
      ),
      inherit.aes = FALSE,
      colour = NA,
      alpha = 0.12
    ) +
    ggplot2::geom_line(
      data = mapping_summary,
      ggplot2::aes(
        y = mean_realised_pfpr_2_10,
        group = backend
      ),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = mapping_summary,
      ggplot2::aes(y = mean_realised_pfpr_2_10),
      size = 2.2
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values = c("#D55E00", "#0072B2")) +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#0072B2")) +
    ggplot2::labs(
      title = "EIR to realised PfPR2-10 mapping",
      subtitle = "One node, no releases, no seasonality; ribbons show the seed range",
      x = "Configured annual EIR (infectious bites per person-year; log scale)",
      y = "Realised microscopy PfPR2-10 (final-window mean)",
      colour = "Mosquito model",
      fill = "Mosquito model"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      legend.position = "top"
    )

  ggplot2::ggsave(
    filename = file.path(out_dir, "legacy_native_eir_pfpr_mapping.png"),
    plot = plot,
    width = 8,
    height = 5,
    units = "in",
    dpi = 180
  )
  invisible(TRUE)
}

load_local_package(pkg_root)

init_eirs <- sort(unique(parse_numeric_vector_env(
  "MSIMGD_MAPPING_EIRS",
  c(0.1, 0.5, 1, 2, 5, 10, 20, 40)
)))
seeds <- unique(parse_integer_vector_env("MSIMGD_MAPPING_SEEDS", 1:3))
timesteps <- as.integer(parse_numeric_env("MSIMGD_MAPPING_TIMESTEPS", 3 * 365))
measurement_days <- as.integer(parse_numeric_env(
  "MSIMGD_MAPPING_MEASUREMENT_DAYS",
  365
))
human_population <- as.integer(parse_numeric_env(
  "MSIMGD_MAPPING_HUMAN_POPULATION",
  1000
))

if (any(init_eirs <= 0)) {
  stop("MSIMGD_MAPPING_EIRS values must all be positive.", call. = FALSE)
}
if (length(seeds) == 0L || anyNA(seeds)) {
  stop("MSIMGD_MAPPING_SEEDS must contain at least one integer.", call. = FALSE)
}
if (timesteps < 1L) {
  stop("MSIMGD_MAPPING_TIMESTEPS must be positive.", call. = FALSE)
}
if (measurement_days < 1L || measurement_days > timesteps) {
  stop("MSIMGD_MAPPING_MEASUREMENT_DAYS must be between 1 and timesteps.",
       call. = FALSE)
}
if (human_population < 1L) {
  stop("MSIMGD_MAPPING_HUMAN_POPULATION must be positive.", call. = FALSE)
}

native_stage_config <- native_stage_config_from_env()

cat(sprintf("[%s] legacy vs native EIR-PfPR2-10 mapping\n",
            format(Sys.time(), "%F %T")))
cat(sprintf("  package root:       %s\n", pkg_root))
cat(sprintf("  output dir:         %s\n", out_dir))
cat(sprintf("  configured EIRs:    %s\n", paste(init_eirs, collapse = ", ")))
cat(sprintf("  seeds:              %s\n", paste(seeds, collapse = ", ")))
cat(sprintf("  timesteps:          %d\n", timesteps))
cat(sprintf("  measurement window: final %d days\n", measurement_days))
cat(sprintf("  human population:   %d\n", human_population))

rows <- vector(
  "list",
  length(init_eirs) * length(seeds) * 2L
)
row_i <- 0L
for (init_eir in init_eirs) {
  for (seed in seeds) {
    for (backend in c("legacy", "native")) {
      row_i <- row_i + 1L
      cat(sprintf(
        "  run %d/%d: backend=%s seed=%d init_EIR=%g\n",
        row_i,
        length(rows),
        backend,
        seed,
        init_eir
      ))
      rows[[row_i]] <- run_one(
        backend = backend,
        seed = seed,
        init_eir = init_eir,
        human_population = human_population,
        timesteps = timesteps,
        measurement_days = measurement_days,
        native_stage_config = native_stage_config
      )
    }
  }
}

mapping_by_seed <- do.call(rbind, rows)
mapping_summary <- summarise_mapping(mapping_by_seed)

native_rows <- mapping_by_seed$backend == "native"
if (!all(mapping_by_seed$native_wildtype_verified[native_rows])) {
  stop("At least one native run failed the singleton-WW verification.", call. = FALSE)
}
if (any(!is.finite(mapping_by_seed$realised_eir)) ||
    any(mapping_by_seed$realised_eir < 0) ||
    any(!is.finite(mapping_by_seed$realised_pfpr_2_10)) ||
    any(mapping_by_seed$realised_pfpr_2_10 < 0) ||
    any(mapping_by_seed$realised_pfpr_2_10 > 1)) {
  stop("At least one run produced an invalid realised EIR or PfPR2-10.",
       call. = FALSE)
}

context <- list(
  scenario = "one_node_wildtype_legacy_vs_native_eir_pfpr_mapping",
  package_root = pkg_root,
  init_eirs = init_eirs,
  seeds = seeds,
  timesteps = timesteps,
  measurement_days = measurement_days,
  human_population = human_population,
  native_genotype = "WW",
  native_stage_config = native_stage_config,
  no_releases = TRUE,
  no_landscape = TRUE,
  no_mosquito_movement = TRUE,
  no_human_movement = TRUE,
  no_seasonality = TRUE,
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
)

by_seed_path <- file.path(out_dir, "eir_pfpr_mapping_by_seed.csv")
summary_path <- file.path(out_dir, "eir_pfpr_mapping_summary.csv")
context_path <- file.path(out_dir, "context.rds")
plot_path <- file.path(out_dir, "legacy_native_eir_pfpr_mapping.png")

utils::write.csv(mapping_by_seed, by_seed_path, row.names = FALSE)
utils::write.csv(mapping_summary, summary_path, row.names = FALSE)
saveRDS(context, context_path)
plot_written <- write_mapping_plot(mapping_by_seed, mapping_summary, out_dir)

cat("Saved outputs:\n")
cat(sprintf("  %s\n", by_seed_path))
cat(sprintf("  %s\n", summary_path))
cat(sprintf("  %s\n", context_path))
if (isTRUE(plot_written)) {
  cat(sprintf("  %s\n", plot_path))
}

cat("Mean mapping by backend and configured EIR:\n")
print(
  mapping_summary[
    ,
    c(
      "backend", "init_eir", "mean_realised_eir",
      "mean_realised_pfpr_2_10", "sd_realised_pfpr_2_10"
    )
  ],
  row.names = FALSE
)
cat("Native wild-type-only verification passed for every native run.\n")
