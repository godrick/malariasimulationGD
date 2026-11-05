#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else
  "test-scripts/four_node_linear/01_plot.R"
example_dir <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(example_dir, "helpers.R"))
out_dir <- file.path(example_dir, "output")
for (arg in commandArgs(TRUE)) {
  if (!grepl("^--out=.+$", arg)) stop("Use --out=PATH.", call. = FALSE)
  out_dir <- path.expand(sub("^--out=", "", arg))
}
four_node_plots(normalizePath(out_dir, mustWork = TRUE))
cat("Plots saved in ", file.path(out_dir, "plots"), "\n", sep = "")
