#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "test-scripts/run_age_structured_eir10_WW_no_release.R"
}
script_path <- normalizePath(script_path, mustWork = TRUE)
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

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
    "malariasimulationGD is not installed. Run: Rscript test-scripts/install_local.R",
    call. = FALSE
  )
}

largest_remainder_counts <- function(counts, sample_size) {
  raw <- sample_size * counts / sum(counts)
  out <- floor(raw)
  remainder <- sample_size - sum(out)
  if (remainder > 0L) {
    add <- order(raw - out, decreasing = TRUE)[seq_len(remainder)]
    out[add] <- out[add] + 1L
  }
  as.integer(out)
}

deathrates_for_target_age_counts <- function(target_counts, age_width_days) {
  n_age <- length(target_counts)
  if (length(age_width_days) != n_age) {
    stop("age_width_days must have the same length as target_counts.",
         call. = FALSE)
  }
  if (any(target_counts <= 0L)) {
    stop("All target age groups must have positive counts.", call. = FALSE)
  }

  aging_rate <- 1 / age_width_days
  aging_rate[[n_age]] <- 0

  deathrates <- rep(1e-5, n_age)
  for (i in 2:n_age) {
    ratio <- target_counts[[i]] / target_counts[[i - 1L]]
    deathrates[[i]] <- aging_rate[[i - 1L]] / ratio - aging_rate[[i]]
  }
  if (any(!is.finite(deathrates)) || any(deathrates <= 0) || any(deathrates >= 1)) {
    stop("Could not derive valid daily death rates for the requested age structure.",
         call. = FALSE)
  }

  deathrates
}

make_native_WW_cube <- function() {
  genotype <- "WW"
  list(
    ih = array(1, dim = c(1L, 1L, 1L), dimnames = list(genotype, genotype, genotype)),
    tau = array(1, dim = c(1L, 1L, 1L), dimnames = list(genotype, genotype, genotype)),
    eta = matrix(1, nrow = 1L, ncol = 1L, dimnames = list(genotype, genotype)),
    b = setNames(1, genotype),
    c = setNames(1, genotype),
    phi = setNames(0.5, genotype),
    omega = setNames(1, genotype),
    xiF = setNames(1, genotype),
    xiM = setNames(1, genotype),
    s = setNames(1, genotype),
    genotypesID = genotype,
    wildType = genotype
  )
}

load_local_package(pkg_root)

set.seed(3997L)

human_population <- 1000L
init_EIR <- 10
timesteps <- 365L
year <- 365

out_dir <- file.path(pkg_root, "test-scripts", "output",
                     "age_structured_eir10_WW_no_release")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Regional age-structure data, reordered youngest to oldest because
# set_demography() expects increasing age groups.
regional_counts <- c(
  "0-9" = 259294L,
  "10-19" = 232677L,
  "20-29" = 189624L,
  "30-39" = 124181L,
  "40-49" = 62304L,
  "50-59" = 44182L,
  "60-69" = 32793L,
  "70-79" = 16059L,
  "80+" = 7648L
)
target_counts <- largest_remainder_counts(regional_counts, human_population)

age_min_years <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
age_max_years <- c(10, 20, 30, 40, 50, 60, 70, 80, 200)
age_min <- age_min_years * year
age_max <- age_max_years * year
age_width <- age_max - age_min
deathrates <- deathrates_for_target_age_counts(target_counts, age_width)

parameters <- malariasimulationGD::get_parameters(list(
  human_population = human_population,
  native_mosquito_backend = TRUE,
  individual_mosquitoes = FALSE,
  model_seasonality = FALSE,
  progress_bar = FALSE,
  cube = make_native_WW_cube(),
  age_group_rendering_min_ages = age_min,
  age_group_rendering_max_ages = age_max,
  prevalence_rendering_min_ages = age_min,
  prevalence_rendering_max_ages = age_max
))

parameters <- malariasimulationGD::set_demography(
  parameters,
  agegroups = age_max,
  timesteps = 0,
  deathrates = matrix(deathrates, nrow = 1L)
)
parameters <- malariasimulationGD::set_equilibrium(
  parameters,
  init_EIR = init_EIR,
  native_total_M = TRUE
)

sim <- malariasimulationGD::run_resumable_simulation(
  timesteps = timesteps,
  parameters = parameters
)
dat <- sim$data
female <- sim$mosquito_genotypes$female
male <- sim$mosquito_genotypes$male

age_columns <- paste0("n_age_", age_min, "_", age_max)
day1_age_counts <- as.integer(dat[1L, age_columns])
final_age_counts <- as.integer(dat[nrow(dat), age_columns])

age_summary <- data.frame(
  age_group = names(regional_counts),
  regional_count = as.integer(regional_counts),
  regional_proportion = as.numeric(regional_counts / sum(regional_counts)),
  target_count_for_1000 = as.integer(target_counts),
  day1_simulated_count = day1_age_counts,
  final_simulated_count = final_age_counts,
  daily_deathrate = deathrates,
  annual_deathrate = deathrates * year,
  stringsAsFactors = FALSE
)

mosquito_summary <- data.frame(
  timesteps = nrow(dat),
  init_EIR = init_EIR,
  solved_total_M = parameters$total_M,
  genotypes = paste(colnames(female), collapse = ","),
  final_female_WW = tail(female[, "WW"], 1),
  final_male_WW = tail(male[, "WW"], 1),
  release_schedule_present = !is.null(attr(dat, "mosquito_release_schedule")),
  stringsAsFactors = FALSE
)

utils::write.csv(dat, file.path(out_dir, "timeseries.csv"), row.names = FALSE)
utils::write.csv(age_summary, file.path(out_dir, "age_structure_summary.csv"),
                 row.names = FALSE)
utils::write.csv(mosquito_summary, file.path(out_dir, "mosquito_summary.csv"),
                 row.names = FALSE)
saveRDS(
  list(
    parameters = parameters,
    final_state = sim$state,
    mosquito_genotypes = sim$mosquito_genotypes,
    age_summary = age_summary,
    mosquito_summary = mosquito_summary
  ),
  file.path(out_dir, "context.rds")
)

print(age_summary)
print(mosquito_summary)
cat("Wrote outputs to:", out_dir, "\n")
