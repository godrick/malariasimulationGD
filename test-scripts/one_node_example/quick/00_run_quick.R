#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 00_run_quick.R
# ------------------------------------------------------------------------------
# Minimal one-node native backend homing-drive example.
#
# Run from the package root:
#   Rscript test-scripts/one_node_example/quick/00_run_quick.R
#
# Outputs:
#   test-scripts/one_node_example/output/quick/timeseries.csv
#   test-scripts/one_node_example/output/quick/carrier_frequency.csv
#   test-scripts/one_node_example/output/quick/release_schedule.csv
#   test-scripts/one_node_example/output/quick/summary.csv
#   test-scripts/one_node_example/output/quick/context.rds
# ------------------------------------------------------------------------------

start_time <- Sys.time()


args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("test-scripts/one_node_example/quick/00_run_quick.R",
                mustWork = TRUE)
}
example_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
pkg_root <- normalizePath(file.path(example_dir, "..", ".."), mustWork = TRUE)
out_dir <- file.path(example_dir, "output", "quick")
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

load_local_package(pkg_root)
source(file.path(example_dir, "config", "homing_drive.R"))

RNG_SEED <- 20260610L
human_population <- 1000L
init_EIR <- 10
timesteps <- 365*120
release_day <- 60*365L
release_count <- 100L
release_genotype <- "HH"
release_sex <- "M"
allele_s <- c(W = 1, H = 0.8, B = 0.2, R = 1)

set.seed(RNG_SEED)

cat(sprintf("[%s] one-node native quick example\n", format(Sys.time(), "%F %T")))
cat(sprintf("  package root: %s\n", pkg_root))
cat(sprintf("  output dir:   %s\n", out_dir))
cat(sprintf("  human pop:    %d\n", human_population))
cat(sprintf("  init EIR:     %.2f annual infectious bites/person/year\n", init_EIR))

cube <- build_one_node_drive_cube(allele_s = allele_s)
cat("  genotype fitness cube$s:\n")
print(round(cube$s, 3))

parameters <- malariasimulationGD::get_parameters(list(
  human_population = human_population,
  native_mosquito_backend = TRUE,
  individual_mosquitoes = TRUE,
  model_seasonality = FALSE,
  human_mobility_enabled = FALSE,
  human_move_probs = NULL,
  human_move_rates = NULL,
  move_probs = matrix(1, 1, 1),
  move_rates = 0,
  bednets = FALSE,
  spraying = FALSE,
  progress_bar = FALSE,
  cube = cube,
  debug_genotypes = TRUE
))

if ("native_total_M" %in% names(formals(malariasimulationGD::set_equilibrium))) {
  parameters <- malariasimulationGD::set_equilibrium(
    parameters,
    init_EIR = init_EIR,
    native_total_M = TRUE
  )
} else {
  # Older branches do not have native_total_M yet. This path keeps the example
  # runnable, but users should enable native_total_M = TRUE when available.
  parameters <- malariasimulationGD::set_equilibrium(parameters, init_EIR = init_EIR)
}

parameters <- malariasimulationGD::set_releases(parameters, list(
  releasesStart = release_day,
  releasesNumber = 1L,
  releaseCount = release_count,
  releaseSex = release_sex,
  releaseGenotype = release_genotype,
  releasesInterval = 0L
))

cat(sprintf("  solved total_M: %.2f adult females\n", parameters$total_M))
cat(sprintf("  release:       %d %s %s mosquitoes on day %d\n",
            release_count, release_genotype, release_sex, release_day))

t0 <- Sys.time()
sim <- malariasimulationGD::run_simulation(
  timesteps = timesteps,
  parameters = parameters
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("[%s] simulation finished in %.1f seconds\n",
            format(Sys.time(), "%F %T"), elapsed))

if (!("pfpr_2_10_lm" %in% names(sim))) {
  sim$pfpr_2_10_lm <- sim$n_detect_lm_730_3650 / pmax(sim$n_age_730_3650, 1)
}

female <- attr(sim, "mosquito_genotype_counts_female")
male <- attr(sim, "mosquito_genotype_counts_male")
if (is.null(female) || is.null(male)) {
  stop("Simulation output is missing native mosquito genotype count attributes.",
       call. = FALSE)
}

drive_genotypes <- grepl("H", colnames(female), fixed = TRUE)
female_drive <- rowSums(female[, drive_genotypes, drop = FALSE])
male_drive <- rowSums(male[, drive_genotypes, drop = FALSE])
female_total <- rowSums(female)
male_total <- rowSums(male)
adult_drive_carriers <- female_drive + male_drive
adult_total <- female_total + male_total

carrier_frequency <- data.frame(
  timestep = sim$timestep,
  adult_drive_carrier_frequency = adult_drive_carriers / pmax(adult_total, 1),
  adult_drive_carriers = adult_drive_carriers,
  adult_total = adult_total,
  female_drive_carrier_frequency = female_drive / pmax(female_total, 1),
  male_drive_carrier_frequency = male_drive / pmax(male_total, 1),
  female_drive_carriers = female_drive,
  female_total = female_total,
  male_drive_carriers = male_drive,
  male_total = male_total,
  stringsAsFactors = FALSE
)

release_schedule <- attr(sim, "mosquito_release_schedule")
if (is.null(release_schedule)) {
  release_schedule <- data.frame(
    timestep = release_day,
    species = parameters$species[[1]],
    sex = release_sex,
    genotype = release_genotype,
    count = release_count,
    stringsAsFactors = FALSE
  )
}

summary <- data.frame(
  human_population = human_population,
  init_EIR = init_EIR,
  timesteps = timesteps,
  release_day = release_day,
  release_count = release_count,
  release_genotype = release_genotype,
  total_M = parameters$total_M,
  native_total_M_equilibrium = isTRUE(parameters$native_total_M_equilibrium),
  peak_adult_drive_carrier_frequency = max(carrier_frequency$adult_drive_carrier_frequency, na.rm = TRUE),
  final_adult_drive_carrier_frequency = utils::tail(carrier_frequency$adult_drive_carrier_frequency, 1),
  cumulative_infections = sum(sim$n_infections, na.rm = TRUE),
  final_pfpr_2_10_lm = utils::tail(sim$pfpr_2_10_lm, 1),
  elapsed_seconds = round(elapsed, 2),
  stringsAsFactors = FALSE
)

context <- list(
  scenario = "one_node_native_homing_drive",
  package_root = pkg_root,
  human_population = human_population,
  init_EIR = init_EIR,
  timesteps = timesteps,
  release_day = release_day,
  release_count = release_count,
  release_genotype = release_genotype,
  release_sex = release_sex,
  allele_s = allele_s,
  cube_s = cube$s,
  drive_genotypes = colnames(female)[drive_genotypes],
  native_total_M_equilibrium = isTRUE(parameters$native_total_M_equilibrium),
  no_landscape = TRUE,
  no_mosquito_movement = TRUE,
  no_human_movement = TRUE,
  no_seasonality = TRUE
)

utils::write.csv(sim, file.path(out_dir, "timeseries.csv"), row.names = FALSE)
utils::write.csv(carrier_frequency, file.path(out_dir, "carrier_frequency.csv"), row.names = FALSE)
utils::write.csv(release_schedule, file.path(out_dir, "release_schedule.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
saveRDS(context, file.path(out_dir, "context.rds"))

cat("Saved outputs:\n")
cat(sprintf("  %s\n", file.path(out_dir, "timeseries.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "carrier_frequency.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "release_schedule.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "summary.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "context.rds")))

threshold_eir <- 0.05
threshold_pfpr <- 1

low_period <- ts$eir_annual_per_person < threshold_eir &
  ts$pfpr_2_10_lm < threshold_pfpr

rebound_days <- ts$timestep[ts$n_infections > 0 & dplyr::lag(low_period, default = FALSE)]

head(rebound_days)

end_time <- Sys.time()