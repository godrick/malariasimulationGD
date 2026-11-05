#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else
  "test-scripts/four_node_linear/00_run.R"
example_dir <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(example_dir, "config.R"))
source(file.path(example_dir, "helpers.R"))

cfg <- four_node_config()
out_dir <- file.path(example_dir, "output")
for (arg in commandArgs(TRUE)) {
  if (arg == "--smoke") {
    cfg$burnin_days <- 30L
    cfg$days <- 30L
    cfg$release_day <- 10L
    cfg$population <- 100L
    next
  }
  if (!grepl("^--(out|library|days|burnin-days|release-day|population|seed)=.+$", arg)) {
    stop("Use --smoke or --out=PATH, --library=PATH, --days=N, --burnin-days=N, --release-day=N, --population=N, --seed=N.",
         call. = FALSE)
  }
  key <- sub("^--([^=]+)=.*$", "\\1", arg)
  value <- sub("^--[^=]+=", "", arg)
  if (key == "out") out_dir <- path.expand(value)
  else if (key == "library") .libPaths(c(normalizePath(path.expand(value), mustWork = TRUE), .libPaths()))
  else cfg[[gsub("-", "_", key)]] <- suppressWarnings(as.numeric(value))
}
validate_four_node_config(cfg)
total_days <- cfg$burnin_days + cfg$days
for (pkg in c("malariasimulationGD", "MGDrivE", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package ", pkg, ". See this example's README.md for installation instructions.", call. = FALSE)
  }
}
if (file.exists(file.path(out_dir, "context.rds"))) {
  stop("Output already exists; choose a fresh directory with --out=PATH.", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
set.seed(cfg$seed)
movement <- four_node_movement(cfg)
parameters <- four_node_parameters(cfg, movement)
# Audit the actual native capacity, including the genotype-aware lifecycle.
# These internal helpers only inspect parameters; they do not change the model.
capacity_at <- getFromNamespace("native_effective_carrying_capacity", "malariasimulationGD")
capacity <- vapply(parameters, capacity_at, numeric(1), species_i = 1L, timestep = 0L)
capacity_end <- vapply(parameters, capacity_at, numeric(1), species_i = 1L, timestep = total_days)
stopifnot(all(is.finite(capacity)), all(capacity > 0),
          max(capacity) - min(capacity) < 1e-8,
          isTRUE(all.equal(capacity, capacity_end)))
context <- list(config = cfg, movement = movement, capacity = capacity,
                initial_total_M = vapply(parameters, function(p) p$total_M, numeric(1)),
                package_version = as.character(utils::packageVersion("malariasimulationGD")),
                package_path = find.package("malariasimulationGD"),
                replicates = 1L, cores = 1L)
saveRDS(context, file.path(out_dir, "context.rds"))
saveRDS(parameters, file.path(out_dir, "parameters.rds"))
utils::write.csv(data.frame(site = cfg$sites, x_km = cfg$x_km, y_km = 0,
                            population = cfg$population, release = cfg$sites %in% cfg$release_sites,
                            native_carrying_capacity = capacity),
                 file.path(out_dir, "nodes.csv"), row.names = FALSE)
utils::write.csv(movement$distance, file.path(out_dir, "distance_km.csv"))
utils::write.csv(movement$probs, file.path(out_dir, "movement_probabilities.csv"))
cat(sprintf("Four equivalent nodes, %d burn-in + %d observed days, one replicate, one core, seed %d.\n",
            cfg$burnin_days, cfg$days, cfg$seed))
cat(sprintf("Release: observed day %d (simulation day %d), %d male HH at each of A and C.\n",
            cfg$release_day, cfg$burnin_days + cfg$release_day, cfg$release_count))
cat(sprintf("Movement rate: %.8f/day; achieved mean distance: %.3f km.\n",
            movement$rates[1L], movement$mean_distance_km))
cat(sprintf("Fixed native carrying capacity per site: %.3f.\n", capacity[1L]))
start <- Sys.time()
sim <- malariasimulationGD::run_metapop_simulation(
  timesteps = total_days, parameters = parameters,
  mixing_tt = 1L, export_mixing = list(diag(4L)), import_mixing = list(diag(4L)),
  p_captured_tt = 1L, p_captured = list(matrix(0, 4L, 4L)), p_success = 0,
  return_state = FALSE, render_output = TRUE, return_summary = TRUE
)
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
saveRDS(sim, file.path(out_dir, "simulation.rds"))
outputs <- four_node_extract(sim$data, cfg)
incidence <- four_node_incidence(outputs$daily, cfg)
utils::write.csv(outputs$carriers, file.path(out_dir, "carrier_frequency.csv"), row.names = FALSE)
utils::write.csv(outputs$daily, file.path(out_dir, "incidence_daily.csv"), row.names = FALSE)
utils::write.csv(incidence, file.path(out_dir, "incidence_windows.csv"), row.names = FALSE)
utils::write.csv(outputs$schedule, file.path(out_dir, "release_schedule.csv"), row.names = FALSE)
utils::write.csv(outputs$burnin_daily, file.path(out_dir, "burnin_incidence_daily.csv"), row.names = FALSE)
utils::write.csv(outputs$burnin_carriers, file.path(out_dir, "burnin_carrier_frequency.csv"), row.names = FALSE)
utils::write.csv(four_node_baseline_diagnostics(outputs, cfg),
                 file.path(out_dir, "baseline_diagnostics.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(cfg$sites, function(site) {
  d <- outputs$daily[outputs$daily$site == site, ]
  c <- outputs$carriers[outputs$carriers$site == site, ]
  first <- c$timestep[c$gd_carriers > 0]
  data.frame(site = site, days = cfg$days, burnin_days = cfg$burnin_days,
             total_simulation_days = total_days, population = cfg$population,
             seed = cfg$seed, replicates = 1L, cores = 1L,
             first_carrier_day = if (length(first)) min(first) else NA_integer_,
             final_carrier_frequency = tail(c$carrier_frequency, 1L),
             infections = sum(d$infections), person_years = sum(d$children) / 365,
             infection_incidence_ppy = 365 * sum(d$infections) / sum(d$children),
             elapsed_seconds = elapsed)
}))
utils::write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
four_node_plots(out_dir)
print(summary, row.names = FALSE)
cat(sprintf("Finished in %.1f simulation seconds. Outputs: %s\n", elapsed, out_dir))
