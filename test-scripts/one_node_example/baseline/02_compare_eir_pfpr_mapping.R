#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 02_compare_eir_pfpr_mapping.R
# ------------------------------------------------------------------------------
# Compare the configured EIR -> realised PfPR2-10 mapping for the legacy and
# native mosquito backends in the one-node, no-release baseline. The native
# model is given a singleton WW cube, so only wild-type mosquitoes can exist.
# The legacy model is the aggregate wild-type baseline and does not track
# mosquito genotypes. Both backends use 45% effective AL clinical-treatment
# coverage and the regional age structure in config/age_structure_summary.csv.
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
#   MSIMGD_MAPPING_NATIVE_NEIP=5
#   MSIMGD_MAPPING_NATIVE_NU=1
#   MSIMGD_MAPPING_NATIVE_TAU_STEP=0.1
# ------------------------------------------------------------------------------

start_time <- Sys.time()
falciparum_b0 <- 0.59
clinical_treatment_coverage <- 0.45
clinical_treatment_start_day <- 1L

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
age_structure_csv <- file.path(example_dir, "config", "age_structure_summary.csv")
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

load_age_structure <- function(path, terminal_age_years = 200) {
  path <- normalizePath(path, mustWork = TRUE)
  age_data <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("age_group", "regional_proportion", "daily_deathrate")
  missing <- setdiff(required, names(age_data))
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "Age-structure CSV is missing required columns: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (nrow(age_data) < 2L || anyDuplicated(age_data$age_group)) {
    stop("Age-structure CSV must contain unique, ordered age groups.", call. = FALSE)
  }

  labels <- trimws(age_data$age_group)
  open_ended <- grepl("\\+$", labels)
  if (sum(open_ended) != 1L || !open_ended[[length(open_ended)]]) {
    stop("The final age group must be the only open-ended group (for example, 80+).",
         call. = FALSE)
  }

  age_min_years <- suppressWarnings(as.numeric(sub("[-+].*$", "", labels)))
  age_max_years <- numeric(length(labels))
  closed <- !open_ended
  age_max_years[closed] <- suppressWarnings(
    as.numeric(sub("^.*-", "", labels[closed])) + 1
  )
  age_max_years[open_ended] <- terminal_age_years

  if (any(!is.finite(age_min_years)) || any(!is.finite(age_max_years)) ||
      age_min_years[[1L]] != 0 ||
      any(age_min_years[-1L] != age_max_years[-length(age_max_years)]) ||
      any(age_max_years <= age_min_years)) {
    stop("Age groups must be contiguous, ordered intervals beginning at age zero.",
         call. = FALSE)
  }

  proportions <- as.numeric(age_data$regional_proportion)
  deathrates <- as.numeric(age_data$daily_deathrate)
  if (any(!is.finite(proportions)) || any(proportions <= 0) ||
      !isTRUE(all.equal(sum(proportions), 1, tolerance = 1e-6))) {
    stop("regional_proportion must contain positive finite values summing to one.",
         call. = FALSE)
  }
  if (any(!is.finite(deathrates)) || any(deathrates <= 0) ||
      any(deathrates >= 1)) {
    stop("daily_deathrate must contain finite probabilities in (0, 1).",
         call. = FALSE)
  }

  list(
    path = path,
    labels = labels,
    regional_proportions = proportions,
    age_min_days = as.integer(age_min_years * 365),
    age_max_days = as.integer(age_max_years * 365),
    daily_deathrates = deathrates,
    source = age_data
  )
}

make_wildtype_cube <- function(b0) {
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
    # cube$b is an absolute mosquito-to-human transmission probability. Match
    # it to the falciparum human-model reference so WW has multiplier b / b0 = 1.
    b = stats::setNames(b0, genotype),
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

base_overrides <- function(backend, human_population, b0) {
  backend <- match.arg(backend, c("legacy", "native"))
  list(
    human_population = human_population,
    b0 = b0,
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
    cube = if (identical(backend, "native")) make_wildtype_cube(b0) else NULL
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
  value_or_default <- function(name, default) {
    value <- Sys.getenv(name, unset = "")
    if (!nzchar(value)) {
      return(default)
    }
    parse_numeric_env(name, default)
  }
  value_or_null <- function(name) {
    value <- Sys.getenv(name, unset = "")
    if (!nzchar(value)) {
      return(NULL)
    }
    parse_numeric_env(name, NA_real_)
  }
  list(
    nE = value_or_default("MSIMGD_MAPPING_NATIVE_NE", 2L),
    nL = value_or_default("MSIMGD_MAPPING_NATIVE_NL", 3L),
    nP = value_or_default("MSIMGD_MAPPING_NATIVE_NP", 2L),
    nEIP = value_or_default("MSIMGD_MAPPING_NATIVE_NEIP", 5L),
    nu = value_or_null("MSIMGD_MAPPING_NATIVE_NU"),
    tau_step = value_or_null("MSIMGD_MAPPING_NATIVE_TAU_STEP")
  )
}

build_parameters <- function(backend, human_population, init_eir, b0,
                             native_stage_config, age_structure,
                             treatment_coverage, treatment_start_day) {
  parameters <- malariasimulationGD::get_parameters(
    base_overrides(backend, human_population, b0)
  )
  if (identical(backend, "native")) {
    parameters <- apply_native_stage_overrides(parameters, native_stage_config)
  }

  parameters <- malariasimulationGD::set_demography(
    parameters,
    agegroups = age_structure$age_max_days,
    timesteps = 0,
    deathrates = matrix(age_structure$daily_deathrates, nrow = 1L)
  )
  parameters <- malariasimulationGD::set_drugs(
    parameters,
    list(malariasimulationGD::AL_params)
  )
  parameters <- malariasimulationGD::set_clinical_treatment(
    parameters,
    drug = 1L,
    timesteps = treatment_start_day,
    coverages = treatment_coverage
  )

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

run_one <- function(backend, seed, init_eir, human_population, b0, timesteps,
                    measurement_days, native_stage_config, age_structure,
                    treatment_coverage, treatment_start_day) {
  parameters <- build_parameters(
    backend = backend,
    human_population = human_population,
    init_eir = init_eir,
    b0 = b0,
    native_stage_config = native_stage_config,
    age_structure = age_structure,
    treatment_coverage = treatment_coverage,
    treatment_start_day = treatment_start_day
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
      min_realised_eir = min(x$realised_eir),
      max_realised_eir = max(x$realised_eir),
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

write_mapping_plots <- function(mapping_by_seed, mapping_summary, out_dir,
                                treatment_coverage) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' is not available; skipping the PNG plots.")
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

  pfpr_plot <- ggplot2::ggplot(
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
      subtitle = sprintf(
        "One node, AL ft=%.2f, no releases or seasonality; ribbons show the seed range",
        treatment_coverage
      ),
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
    plot = pfpr_plot,
    width = 8,
    height = 5,
    units = "in",
    dpi = 180
  )

  eir_plot <- ggplot2::ggplot(
    mapping_by_seed,
    ggplot2::aes(x = init_eir, y = realised_eir, colour = backend)
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      colour = "grey35",
      linetype = "dashed",
      linewidth = 0.65
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
        y = mean_realised_eir,
        ymin = min_realised_eir,
        ymax = max_realised_eir,
        fill = backend,
        group = backend
      ),
      inherit.aes = FALSE,
      colour = NA,
      alpha = 0.12
    ) +
    ggplot2::geom_line(
      data = mapping_summary,
      ggplot2::aes(y = mean_realised_eir, group = backend),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = mapping_summary,
      ggplot2::aes(y = mean_realised_eir),
      size = 2.2
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = c("#D55E00", "#0072B2")) +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#0072B2")) +
    ggplot2::labs(
      title = "Configured versus realised EIR",
      subtitle = paste(
        sprintf("One node, AL ft=%.2f, no releases or seasonality;", treatment_coverage),
        "dashed line is 1:1 and ribbons show the seed range"
      ),
      x = "Configured annual EIR (infectious bites per person-year; log scale)",
      y = "Realised annual EIR (final-window total per person; log scale)",
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
    filename = file.path(out_dir, "legacy_native_configured_vs_realised_eir.png"),
    plot = eir_plot,
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
age_structure <- load_age_structure(age_structure_csv)
aquatic_stage_counts <- unlist(native_stage_config[c("nE", "nL", "nP")])
if (any(aquatic_stage_counts != as.integer(aquatic_stage_counts)) ||
    any(aquatic_stage_counts < 1L)) {
  stop("Native nE, nL, and nP must be positive integers.", call. = FALSE)
}
if (native_stage_config$nEIP != as.integer(native_stage_config$nEIP) ||
    native_stage_config$nEIP < 0L) {
  stop("Native nEIP must be a non-negative integer.", call. = FALSE)
}

cat(sprintf("[%s] legacy vs native EIR-PfPR2-10 mapping\n",
            format(Sys.time(), "%F %T")))
cat(sprintf("  package root:       %s\n", pkg_root))
cat(sprintf("  output dir:         %s\n", out_dir))
cat(sprintf("  configured EIRs:    %s\n", paste(init_eirs, collapse = ", ")))
cat(sprintf("  seeds:              %s\n", paste(seeds, collapse = ", ")))
cat(sprintf("  timesteps:          %d\n", timesteps))
cat(sprintf("  measurement window: final %d days\n", measurement_days))
cat(sprintf("  human population:   %d\n", human_population))
cat(sprintf("  falciparum b0/WW b: %.2f\n", falciparum_b0))
cat(sprintf(
  "  clinical treatment: AL with ft=%.2f from day %d\n",
  clinical_treatment_coverage,
  clinical_treatment_start_day
))
cat(sprintf("  age demography:     %s\n", age_structure$path))
cat(sprintf(
  "  native stages:      nE=%d, nL=%d, nP=%d, nEIP=%d\n",
  native_stage_config$nE,
  native_stage_config$nL,
  native_stage_config$nP,
  native_stage_config$nEIP
))

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
        b0 = falciparum_b0,
        timesteps = timesteps,
        measurement_days = measurement_days,
        native_stage_config = native_stage_config,
        age_structure = age_structure,
        treatment_coverage = clinical_treatment_coverage,
        treatment_start_day = clinical_treatment_start_day
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
  falciparum_b0 = falciparum_b0,
  native_wildtype_b = falciparum_b0,
  native_genotype = "WW",
  native_stage_config = native_stage_config,
  clinical_treatment = list(
    drug = "AL",
    effective_coverage = clinical_treatment_coverage,
    start_day = clinical_treatment_start_day
  ),
  age_structure = list(
    csv = age_structure$path,
    labels = age_structure$labels,
    regional_proportions = age_structure$regional_proportions,
    age_min_days = age_structure$age_min_days,
    age_max_days = age_structure$age_max_days,
    daily_deathrates = age_structure$daily_deathrates,
    expected_counts = human_population * age_structure$regional_proportions
  ),
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
eir_plot_path <- file.path(out_dir, "legacy_native_configured_vs_realised_eir.png")

utils::write.csv(mapping_by_seed, by_seed_path, row.names = FALSE)
utils::write.csv(mapping_summary, summary_path, row.names = FALSE)
saveRDS(context, context_path)
plots_written <- write_mapping_plots(
  mapping_by_seed,
  mapping_summary,
  out_dir,
  clinical_treatment_coverage
)

cat("Saved outputs:\n")
cat(sprintf("  %s\n", by_seed_path))
cat(sprintf("  %s\n", summary_path))
cat(sprintf("  %s\n", context_path))
if (isTRUE(plots_written)) {
  cat(sprintf("  %s\n", plot_path))
  cat(sprintf("  %s\n", eir_plot_path))
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
