#!/usr/bin/env Rscript

# Standalone native mosquito lifecycle example, parameterized from the
# MGDrivE2 single-node lifecycle example.

# if (!requireNamespace("malariasimulationGD", quietly = TRUE)) {
#   stop("Install malariasimulationGD first, for example: Rscript test-scripts/install_local.R",
#        call. = FALSE)
# }
# if (!requireNamespace("MGDrivE", quietly = TRUE)) {
#   stop("Install MGDrivE first, for example: Rscript test-scripts/install_local.R",
#        call. = FALSE)
# }

library(MGDrivE)

cube <- MGDrivE::cubeMendelian()

NF <- 500
theta <- list(
  qE = 1 / 4, nE = 2,
  qL = 1 / 3, nL = 3,
  qP = 1 / 6, nP = 2,
  muE = 0.05, muL = 0.15, muP = 0.05,
  muF = 0.09, muM = 0.09,
  beta = 16,
  nu = 1 / (4 / 24)
)

parameters <- malariasimulationGD::get_mosquito_parameters(
  overrides = list(
    total_M = NF,
    native_mosquito_nE = theta$nE,
    native_mosquito_nL = theta$nL,
    native_mosquito_nP = theta$nP,
    del = theta$nE / theta$qE,
    dl = theta$nL / theta$qL,
    dpl = theta$nP / theta$qP,
    me = theta$muE,
    ml = theta$muL,
    mup = theta$muP,
    mum = theta$muF,
    beta = theta$beta,
    native_mosquito_nu = theta$nu,
    mosquito_tau_step = 0.1
  ),
  cube = cube
)

release_genotype <- if (!is.null(cube$releaseType)) cube$releaseType else tail(cube$genotypesID, 1)
release_times <- seq(from = 20L, length.out = 5L, by = 10L)
releases <- data.frame(
  timestep = release_times,
  node = "node1",
  species = "gamb",
  stage = "adult_female",
  genotype = release_genotype,
  count = 50,
  stringsAsFactors = FALSE
)

tmax <- 125L

ode_out <- malariasimulationGD::run_mosquito_simulation(
  timesteps = tmax,
  parameters = parameters,
  sampler = "ode",
  releases = releases
)

tau_out <- malariasimulationGD::run_mosquito_simulation(
  timesteps = tmax,
  parameters = parameters,
  sampler = "tau",
  releases = releases,
  seed = 123
)

adult_counts <- rbind(
  transform(
    subset(ode_out$counts, stage %in% c("adult_female", "adult_male")),
    sampler = "ode"
  ),
  transform(
    subset(tau_out$counts, stage %in% c("adult_female", "adult_male")),
    sampler = "tau"
  )
)

stopifnot(
  all(c("egg", "larva", "pupa", "adult_male", "adult_female") %in% unique(ode_out$counts$stage)),
  all(tau_out$state_counts$count >= 0),
  all(abs(tau_out$state_counts$count - round(tau_out$state_counts$count)) < 1e-8)
)

print(utils::head(adult_counts, 12))
cat("Final adult counts by sampler/stage/genotype:\n")
print(
  stats::aggregate(
    count ~ sampler + stage + genotype,
    data = subset(adult_counts, timestep == tmax),
    FUN = sum
  )
)

if (interactive() && requireNamespace("ggplot2", quietly = TRUE)) {
  p <- ggplot2::ggplot(adult_counts) +
    ggplot2::geom_line(ggplot2::aes(x = timestep, y = count, color = genotype)) +
    ggplot2::facet_grid(stage ~ sampler, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Standalone Native Mosquito Lifecycle")
  print(p)
}
