#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
repo <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
args <- commandArgs(TRUE)
if (length(args) > 1L) {
  stop("Usage: Rscript test-scripts/tests/run_tests.R [test-file-filter]", call. = FALSE)
}
for (package in c("testthat", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Developer tests require the '", package, "' package.", call. = FALSE)
  }
}

set.seed(20261105)
pkgload::load_all(repo, quiet = TRUE, helpers = FALSE)
testthat::test_dir(
  file.path(dirname(script), "testthat"),
  filter = if (length(args)) args[[1L]] else NULL,
  package = "malariasimulationGD",
  reporter = "summary",
  stop_on_failure = TRUE
)
