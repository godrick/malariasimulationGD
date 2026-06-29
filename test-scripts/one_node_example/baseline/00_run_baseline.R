#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00_run_baseline.R
# ------------------------------------------------------------------------------
# High-transmission, one-node, no-release baseline diagnostic for the native
# mosquito backend with genotype tracking enabled.
#
# Run from the package root:
#   Rscript test-scripts/one_node_example/baseline/00_run_baseline.R
#
# Optional environment overrides:
#   MSIMGD_BASELINE_SEEDS=1,2,3,4,5
#   MSIMGD_BASELINE_TIMESTEPS=1095
#   MSIMGD_BASELINE_INIT_EIR=10
#   MSIMGD_BASELINE_HUMAN_POPULATION=1000
# ------------------------------------------------------------------------------

start_time <- Sys.time()

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/one_node_example/baseline/00_run_baseline.R",
    mustWork = TRUE
  )
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
pkg_root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)
out_dir <- file.path(example_dir, "output", "baseline")
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

build_no_release_baseline_cube <- function() {
  genotypes <- c("WW", "HH")
  G <- length(genotypes)
  ih <- array(
    0,
    dim = c(G, G, G),
    dimnames = list(genotypes, genotypes, genotypes)
  )
  ih["WW", "WW", "WW"] <- 1
  ih["WW", "HH", "HH"] <- 1
  ih["HH", "WW", "HH"] <- 1
  ih["HH", "HH", "HH"] <- 1

  list(
    ih = ih,
    tau = array(1, dim = c(G, G, G), dimnames = list(genotypes, genotypes, genotypes)),
    eta = matrix(1, nrow = G, ncol = G, dimnames = list(genotypes, genotypes)),
    b = stats::setNames(rep(1, G), genotypes),
    c = stats::setNames(rep(1, G), genotypes),
    phi = stats::setNames(rep(0.5, G), genotypes),
    omega = stats::setNames(rep(1, G), genotypes),
    xiF = stats::setNames(rep(1, G), genotypes),
    xiM = stats::setNames(rep(1, G), genotypes),
    s = stats::setNames(rep(1, G), genotypes),
    genotypesID = genotypes,
    wildType = "WW"
  )
}

build_no_release_baseline_parameters <- function(human_population, init_EIR) {
  parameters <- malariasimulationGD::get_parameters(list(
    human_population = human_population,
    native_mosquito_backend = TRUE,
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
    cube = build_no_release_baseline_cube()
  ))

  parameters$incidence_rendering_min_ages <- 0
  parameters$incidence_rendering_max_ages <- 100 * 365 - 1
  parameters$clinical_incidence_rendering_min_ages <- c(0, 5 * 365)
  parameters$clinical_incidence_rendering_max_ages <- c(5 * 365 - 1, 100 * 365 - 1)

  malariasimulationGD::set_equilibrium(
    parameters,
    init_EIR = init_EIR,
    native_total_M = TRUE
  )
}

check_row <- function(seed, check, passed, value = NA_real_, lower = NA_real_,
                      upper = NA_real_, details = "") {
  data.frame(
    seed = seed,
    check = check,
    passed = isTRUE(passed),
    value = as.numeric(value),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    details = details,
    stringsAsFactors = FALSE
  )
}

all_finite_nonnegative <- function(data, columns) {
  present <- intersect(columns, names(data))
  if (length(present) == 0L) {
    return(FALSE)
  }
  values <- as.matrix(data[present])
  all(is.finite(values)) && all(values >= 0)
}

run_one_seed <- function(seed, parameters, timesteps) {
  set.seed(seed)
  run <- malariasimulationGD::run_resumable_simulation(
    timesteps = timesteps,
    parameters = parameters
  )
  data <- run$data
  female <- run$mosquito_genotypes$female
  male <- run$mosquito_genotypes$male
  aquatic_E <- attr(data, "mosquito_aquatic_genotype_E")
  aquatic_L <- attr(data, "mosquito_aquatic_genotype_L")
  aquatic_P <- attr(data, "mosquito_aquatic_genotype_P")

  final_year <- utils::tail(data, min(365L, nrow(data)))
  adult_female_total <- rowSums(female)
  adult_male_total <- rowSums(male)
  adult_total <- adult_female_total + adult_male_total
  non_wt_adults <- female[, "HH"] + male[, "HH"]
  data$seed <- seed
  data$pfpr_2_10_lm <- data$n_detect_lm_730_3650 / pmax(data$n_age_730_3650, 1)
  data$annualized_daily_eir <- 365 * data$EIR_gamb / parameters$human_population
  data$adult_female_total <- adult_female_total
  data$adult_male_total <- adult_male_total
  data$adult_total <- adult_total
  data$non_wt_adults <- non_wt_adults
  data$non_wt_frequency <- non_wt_adults / pmax(adult_total, 1)

  final_idx <- tail(seq_len(nrow(data)), nrow(final_year))
  annual_eir <- sum(data$EIR_gamb[final_idx]) / parameters$human_population
  pfpr_2_10_lm <- mean(data$pfpr_2_10_lm[final_idx])
  annual_infection_incidence <- sum(data$n_inc_0_36499[final_idx]) /
    parameters$human_population
  under5_clinical_incidence <- sum(data$n_inc_clinical_0_1824[final_idx]) /
    mean(data$n_age_0_1824[final_idx])
  adult_female_cv <- stats::sd(data$total_M_gamb[final_idx]) /
    mean(data$total_M_gamb[final_idx])
  human_total <- data$S_count + data$A_count + data$D_count + data$U_count + data$Tr_count
  release_columns <- grep("^n_released_", names(data), value = TRUE)
  aquatic_non_wt <- sum(aquatic_E[, "HH"]) + sum(aquatic_L[, "HH"]) + sum(aquatic_P[, "HH"])

  summary <- data.frame(
    seed = seed,
    timesteps = timesteps,
    init_EIR = parameters$init_EIR,
    human_population = parameters$human_population,
    total_M = parameters$total_M,
    annual_eir_final_year = annual_eir,
    pfpr_2_10_lm_final_year = pfpr_2_10_lm,
    annual_infection_incidence_final_year = annual_infection_incidence,
    under5_clinical_incidence_final_year = under5_clinical_incidence,
    adult_female_cv_final_year = adult_female_cv,
    max_non_wt_frequency = max(data$non_wt_frequency),
    cumulative_infections = sum(data$n_infections),
    stringsAsFactors = FALSE
  )

  checks <- rbind(
    check_row(
      seed,
      "no_release_schedule",
      is.null(attr(data, "mosquito_release_schedule")),
      value = as.numeric(!is.null(attr(data, "mosquito_release_schedule"))),
      details = "No mosquito release schedule should be attached."
    ),
    check_row(
      seed,
      "no_release_columns",
      length(release_columns) == 0L,
      value = length(release_columns),
      details = "No n_released_* columns should be present."
    ),
    check_row(
      seed,
      "adult_non_wildtype_zero",
      sum(non_wt_adults) == 0,
      value = sum(non_wt_adults),
      details = "HH adults should never appear without a release."
    ),
    check_row(
      seed,
      "aquatic_non_wildtype_zero",
      aquatic_non_wt == 0,
      value = aquatic_non_wt,
      details = "HH aquatic states should never appear without a release."
    ),
    check_row(
      seed,
      "female_genotypes_match_rendered_adults",
      isTRUE(all.equal(adult_female_total, data$total_M_gamb, tolerance = 1e-8)),
      value = max(abs(adult_female_total - data$total_M_gamb)),
      upper = 1e-8,
      details = "Wildtype female genotype totals should match rendered adult females."
    ),
    check_row(
      seed,
      "male_totals_match_native_baseline_ratio",
      isTRUE(all.equal(adult_male_total, adult_female_total, tolerance = 1e-8)),
      value = max(abs(adult_male_total - adult_female_total)),
      upper = 1e-8,
      details = "No-release native baseline uses a 1:1 adult male:female abundance."
    ),
    check_row(
      seed,
      "human_compartment_total_constant",
      all(human_total == parameters$human_population),
      value = max(abs(human_total - parameters$human_population)),
      upper = 0,
      details = "S + A + D + U + Tr should equal the configured population every day."
    ),
    check_row(
      seed,
      "key_outputs_finite_nonnegative",
      all_finite_nonnegative(
        data,
        c(
          "E_gamb_count", "L_gamb_count", "P_gamb_count",
          "Sm_gamb_count", "Pm_gamb_count", "Im_gamb_count", "total_M_gamb",
          "EIR_gamb", "FOIM_gamb", "mu_gamb", "infectivity",
          "n_infections", "n_bitten",
          "S_count", "A_count", "D_count", "U_count", "Tr_count",
          "n_age_730_3650", "n_detect_lm_730_3650", "n_detect_pcr_730_3650",
          "n_inc_0_36499", "n_inc_clinical_0_1824"
        )
      ),
      details = "Key mosquito, human, prevalence, and incidence columns should be finite and non-negative."
    ),
    check_row(seed, "annual_eir_high_transmission", annual_eir > 5 && annual_eir < 20,
              value = annual_eir, lower = 5, upper = 20),
    check_row(seed, "pfpr_2_10_lm_endemic", pfpr_2_10_lm > 0.25 && pfpr_2_10_lm < 0.85,
              value = pfpr_2_10_lm, lower = 0.25, upper = 0.85),
    check_row(seed, "annual_infection_incidence_plausible",
              annual_infection_incidence > 1 && annual_infection_incidence < 10,
              value = annual_infection_incidence, lower = 1, upper = 10),
    check_row(seed, "under5_clinical_incidence_plausible",
              under5_clinical_incidence > 0.2 && under5_clinical_incidence < 4,
              value = under5_clinical_incidence, lower = 0.2, upper = 4),
    check_row(seed, "adult_female_baseline_stable", adult_female_cv < 1e-4,
              value = adult_female_cv, lower = 0, upper = 1e-4)
  )

  list(timeseries = data, summary = summary, checks = checks)
}

write_optional_plots <- function(timeseries, out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' is not available; skipping optional PNG plots.")
    return(invisible(FALSE))
  }

  library(ggplot2)
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
    "annualized_eir_over_time.png",
    ggplot(timeseries, aes(timestep, annualized_daily_eir, color = seed)) +
      geom_line(linewidth = 0.6) +
      labs(title = "Annualized Daily EIR", x = "Day", y = "Annualized EIR/person") +
      base_theme
  )
  write_plot(
    "pfpr_2_10_lm_over_time.png",
    ggplot(timeseries, aes(timestep, pfpr_2_10_lm, color = seed)) +
      geom_line(linewidth = 0.6) +
      labs(title = "PfPR 2-10 by Microscopy", x = "Day", y = "PfPR 2-10") +
      base_theme
  )
  write_plot(
    "incidence_over_time.png",
    ggplot(timeseries, aes(timestep, n_infections, color = seed)) +
      geom_line(linewidth = 0.6) +
      labs(title = "Daily New Infections", x = "Day", y = "New infections") +
      base_theme
  )

  adult_counts <- rbind(
    data.frame(
      seed = timeseries$seed,
      timestep = timeseries$timestep,
      sex = "Female",
      adults = timeseries$adult_female_total
    ),
    data.frame(
      seed = timeseries$seed,
      timestep = timeseries$timestep,
      sex = "Male",
      adults = timeseries$adult_male_total
    )
  )
  write_plot(
    "adult_totals_over_time.png",
    ggplot(adult_counts, aes(timestep, adults, color = sex, linetype = seed)) +
      geom_line(linewidth = 0.6) +
      labs(title = "Adult Mosquito Totals", x = "Day", y = "Adults") +
      base_theme
  )
  write_plot(
    "non_wildtype_frequency_over_time.png",
    ggplot(timeseries, aes(timestep, non_wt_frequency, color = seed)) +
      geom_line(linewidth = 0.6) +
      labs(title = "Non-Wildtype Adult Frequency", x = "Day", y = "HH adult frequency") +
      base_theme
  )

  invisible(TRUE)
}

load_local_package(pkg_root)

seeds <- parse_seed_env("MSIMGD_BASELINE_SEEDS", 1:5)
timesteps <- as.integer(parse_numeric_env("MSIMGD_BASELINE_TIMESTEPS", 3 * 365))
init_EIR <- parse_numeric_env("MSIMGD_BASELINE_INIT_EIR", 10)
human_population <- as.integer(parse_numeric_env("MSIMGD_BASELINE_HUMAN_POPULATION", 1000))

if (timesteps < 365L) {
  stop("MSIMGD_BASELINE_TIMESTEPS must be at least 365.", call. = FALSE)
}
if (human_population <= 0L) {
  stop("MSIMGD_BASELINE_HUMAN_POPULATION must be positive.", call. = FALSE)
}

cat(sprintf("[%s] high-transmission no-release baseline diagnostic\n", format(Sys.time(), "%F %T")))
cat(sprintf("  package root: %s\n", pkg_root))
cat(sprintf("  output dir:   %s\n", out_dir))
cat(sprintf("  seeds:        %s\n", paste(seeds, collapse = ", ")))
cat(sprintf("  timesteps:    %d\n", timesteps))
cat(sprintf("  human pop:    %d\n", human_population))
cat(sprintf("  init EIR:     %.2f annual infectious bites/person/year\n", init_EIR))

parameters <- build_no_release_baseline_parameters(
  human_population = human_population,
  init_EIR = init_EIR
)

runs <- lapply(seeds, run_one_seed, parameters = parameters, timesteps = timesteps)
timeseries <- do.call(rbind, lapply(runs, `[[`, "timeseries"))
summary_by_seed <- do.call(rbind, lapply(runs, `[[`, "summary"))
baseline_checks <- do.call(rbind, lapply(runs, `[[`, "checks"))

context <- list(
  scenario = "one_node_native_no_release_high_transmission_baseline",
  package_root = pkg_root,
  seeds = seeds,
  timesteps = timesteps,
  init_EIR = init_EIR,
  human_population = human_population,
  total_M = parameters$total_M,
  native_total_M_equilibrium = isTRUE(parameters$native_total_M_equilibrium),
  genotypes = parameters$cube$genotypesID,
  wild_type = parameters$cube$wildType,
  no_releases = TRUE,
  no_landscape = TRUE,
  no_mosquito_movement = TRUE,
  no_human_movement = TRUE,
  no_seasonality = TRUE,
  thresholds = list(
    annual_eir = c(lower = 5, upper = 20),
    pfpr_2_10_lm = c(lower = 0.25, upper = 0.85),
    annual_infection_incidence = c(lower = 1, upper = 10),
    under5_clinical_incidence = c(lower = 0.2, upper = 4),
    adult_female_cv = c(lower = 0, upper = 1e-4)
  ),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
)

utils::write.csv(timeseries, file.path(out_dir, "timeseries.csv"), row.names = FALSE)
utils::write.csv(summary_by_seed, file.path(out_dir, "summary_by_seed.csv"), row.names = FALSE)
utils::write.csv(baseline_checks, file.path(out_dir, "baseline_checks.csv"), row.names = FALSE)
saveRDS(context, file.path(out_dir, "context.rds"))
write_optional_plots(timeseries, out_dir)

cat("Saved outputs:\n")
cat(sprintf("  %s\n", file.path(out_dir, "timeseries.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "summary_by_seed.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "baseline_checks.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "context.rds")))

cat("Final-year metric ranges across seeds:\n")
cat(sprintf(
  "  annual EIR/person:             %.3f - %.3f\n",
  min(summary_by_seed$annual_eir_final_year),
  max(summary_by_seed$annual_eir_final_year)
))
cat(sprintf(
  "  PfPR 2-10 microscopy:          %.3f - %.3f\n",
  min(summary_by_seed$pfpr_2_10_lm_final_year),
  max(summary_by_seed$pfpr_2_10_lm_final_year)
))
cat(sprintf(
  "  infection incidence/person-yr: %.3f - %.3f\n",
  min(summary_by_seed$annual_infection_incidence_final_year),
  max(summary_by_seed$annual_infection_incidence_final_year)
))
cat(sprintf(
  "  under-5 clinical incidence/yr: %.3f - %.3f\n",
  min(summary_by_seed$under5_clinical_incidence_final_year),
  max(summary_by_seed$under5_clinical_incidence_final_year)
))

if (!all(baseline_checks$passed)) {
  failed <- baseline_checks[!baseline_checks$passed, c("seed", "check", "value", "lower", "upper")]
  print(failed, row.names = FALSE)
  quit(status = 1)
}

cat("All baseline checks passed.\n")
