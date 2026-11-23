#!/usr/bin/env Rscript

script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
repo <- normalizePath(file.path(dirname(script), "../.."))
source(file.path(dirname(script), "comparison.R"))
mc_require()
args <- commandArgs(TRUE)
output <- if (length(args)) args[[1L]] else tempfile("mosquito-comparison-")
days <- if (length(args) > 1L) as.integer(args[[2L]]) else 120L
if (length(days) != 1L || is.na(days) || days < 90L) stop("Use at least 90 days.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output)
if (file.exists(file.path(output, "summary.csv"))) {
  stop("Output already contains a comparison; choose a new directory.", call. = FALSE)
}
set.seed(20261121L)
scenarios <- mc_scenarios()
results <- tight_results <- vector("list", length(scenarios))
for (i in seq_along(scenarios)) {
  cat(sprintf("[%02d/%02d] %s\n", i, length(scenarios), scenarios[[i]]$label))
  flush.console()
  results[[i]] <- mc_run_scenario(scenarios[[i]], days = days)
  tight_results[[i]] <- mc_run_scenario(scenarios[[i]], days = days, rtol = 1e-10, atol = 1e-12)
}
summary <- do.call(rbind, lapply(results, `[[`, "summary"))
tight <- do.call(rbind, lapply(tight_results, `[[`, "summary"))
convergence <- do.call(rbind, lapply(seq_along(results), function(i) {
  a <- results[[i]]$states
  b <- tight_results[[i]]$states
  data.frame(scenario = summary$scenario[i],
             mg2_max_abs_change = max(abs(a$mg2 - b$mg2)),
             mg4_max_abs_change = max(abs(a$mg4 - b$mg4)),
             standard_cross_engine_error = summary$max_compartment_abs_difference[i],
             tight_cross_engine_error = tight$max_compartment_abs_difference[i])
}))
write <- function(x, name) utils::write.csv(x, file.path(output, name), row.names = FALSE)
write(summary, "summary.csv")
write(tight, "summary_tight.csv")
write(convergence, "solver_convergence.csv")
bind <- function(field) do.call(rbind, lapply(results, function(x) {
  if (is.null(x[[field]])) return(NULL)
  transform(x[[field]], scenario = x$summary$scenario)
}))
write(bind("states"), "compartment_trajectories.csv")
write(bind("totals"), "stage_trajectories.csv")
write(bind("genotypes"), "genotype_trajectories.csv")
write(bind("rhs"), "derivative_comparison.csv")
write(bind("releases"), "release_audit.csv")
saveRDS(list(scenarios = scenarios, results = results, tight_results = tight_results),
        file.path(output, "comparison.rds"))
mc_plot_results(results, output)
revision <- system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE)
writeLines(c(paste("Source revision:", revision),
             "Execution: serial, deterministic, one run per scenario per solver precision.",
             "No human model, infection, stochastic simulation or household intervention model.",
             "Interventions are matched direct mosquito rate changes, not LLIN/IRS efficacy estimates.",
             "Rates on days 30-59 apply to intervals [30,31) through [59,60).",
             "Releases at day d are applied after integration to d and before recording d.",
             "Initial compartments and K independently checked against both equilibrium constructors.",
             "MGDrivE2 uses its own hazards and stoichiometry with deSolve LSODA.",
             "MGDrivE4 uses the installed native entomology engine and Boost Dopri5 solver.",
             "Main criterion, every compartment/time: abs(error) <= 1e-4 + 1e-5*max(abs(counts)).",
             "RHS criterion, matched scenarios: max scaled derivative difference < 1e-8.",
             "Standard rtol/atol = 1e-8/1e-10; tight = 1e-10/1e-12.",
             paste("malariasimulationGD library:", find.package("malariasimulationGD")),
             capture.output(sessionInfo())), file.path(output, "provenance.txt"))
print(summary[, c("scenario", "max_compartment_abs_difference", "max_stage_relative_difference",
                   "rhs_max_scaled_difference", "overlap_pass", "expectation_pass")], row.names = FALSE)
cat("\nResults:", output, "\n")
expected <- summary$expected == "overlap"
stopifnot(all(summary$expectation_pass), all(tight$expectation_pass),
          all(summary$rhs_max_scaled_difference[expected] < 1e-8),
          all(tight$rhs_max_scaled_difference[expected] < 1e-8))
