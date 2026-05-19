#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "test-scripts/install_local.R"
}
pkg_root <- normalizePath(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), ".."),
                          mustWork = TRUE)

repos <- c(
  mrcide = "https://mrc-ide.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
)
options(repos = repos)

runtime_packages <- c(
  "individual",
  "malariaEquilibrium",
  "malariaEquilibriumVivax",
  "Rcpp",
  "RcppArmadillo",
  "statmod",
  "MASS",
  "dqrng",
  "sitmo",
  "BH",
  "R6",
  "progress",
  "ggplot2",
  "dplyr",
  "scales",
  "patchwork",
  "MGDrivE"
)
local_installed <- rownames(utils::installed.packages())
missing_packages <- setdiff(runtime_packages, local_installed)
if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = repos)
}

install.packages(
  pkg_root,
  repos = NULL,
  type = "source"
)

cat("Installed malariasimulationGD into default library path.\n")
