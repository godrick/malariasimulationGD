#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "test-scripts/run_minimal_release_example.R"
}
script_path <- normalizePath(script_path, mustWork = TRUE)
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (!requireNamespace("malariasimulationGD", quietly = TRUE)) {
  stop(
    "malariasimulationGD is not installed. Run: Rscript test-scripts/install_local.R",
    call. = FALSE
  )
}
if (!requireNamespace("MGDrivE", quietly = TRUE)) {
  stop(
    "MGDrivE is not installed. Run: Rscript test-scripts/install_local.R",
    call. = FALSE
  )
}

out_dir <- file.path(pkg_root, "test-scripts", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260514)

params <- malariasimulationGD::get_parameters(list(
  # Use the native count-based tau-leap mosquito engine. We explicitly set
  # individual_mosquitoes = FALSE so it is unambiguous that we are not
  # routing through the legacy event-based individual-mosquito backend.
  native_mosquito_backend = TRUE,
  individual_mosquitoes = FALSE,
  human_population = 50,
  total_M = 200,
  init_foim = 0,
  progress_bar = FALSE
))
params <- malariasimulationGD::parameterise_total_M(params, params$total_M)

cube <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
cube$releaseType <- "aa"
params$cube <- cube

params <- malariasimulationGD::set_releases(params, list(
  releasesStart = 20,
  releasesNumber = 1,
  releaseCount = 200,
  releaseSex = "M"
))

sim <- malariasimulationGD::run_resumable_simulation(60, parameters = params)
dat <- sim$data
geno <- sim$mosquito_genotypes
release_schedule <- attr(dat, "mosquito_release_schedule")

adult_aa <- geno$female[, "aa"] + geno$male[, "aa"]
adult_Aa <- geno$female[, "Aa"] + geno$male[, "Aa"]
first_adult_Aa <- which(adult_Aa > 0)
first_adult_Aa_day <- if (length(first_adult_Aa) > 0L) first_adult_Aa[[1L]] else NA_integer_

summary <- data.frame(
  timesteps = nrow(dat),
  release_day = release_schedule$timestep[[1L]],
  released_genotype = release_schedule$genotype[[1L]],
  released_count = release_schedule$count[[1L]],
  adult_aa_on_release_day = adult_aa[release_schedule$timestep[[1L]]],
  first_adult_Aa_day = first_adult_Aa_day,
  final_total_adult_females = tail(dat$Sm_gamb_count + dat$Pm_gamb_count + dat$Im_gamb_count, 1),
  stringsAsFactors = FALSE
)

utils::write.csv(dat, file.path(out_dir, "minimal_release_timeseries.csv"), row.names = FALSE)
utils::write.csv(release_schedule, file.path(out_dir, "minimal_release_schedule.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "minimal_release_summary.csv"), row.names = FALSE)

print(summary)
cat("Wrote outputs to:", out_dir, "\n")
