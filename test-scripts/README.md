# Developer scenarios and scientific validation

This directory is the developer workspace for exploring scenarios, calibration,
scientific questions, and model validity. Its scripts, examples, local data,
plots, and generated outputs are excluded from `R CMD build` by `.Rbuildignore`
and do not ship in the final package. Keep generated outputs out of version
control and avoid overwriting existing scientific results.

Package regressions belong in `tests/testthat/` and must run independently of
this workspace. Tests specific to developer workflows live in `tests/` within
this directory; see [the developer test runner](tests/README.md).

Dependencies used only here belong in each workflow's documentation, rather
than the package `DESCRIPTION`. Install them in an isolated development library
when needed. `install_local.R` installs into the default library, so review that
choice before running it.

Examples:

- `install_local.R`: installs package dependencies and this local package into
  your default R library.
- `run_minimal_release_example.R`: runs a small generated-data mosquito genotype
  release example and writes outputs to `test-scripts/output`.
- `seven_node_example/`: runs a compact seven-village metapopulation example
  with movement, seasonality, a synthetic contact multiplier, and a homing-drive
  release. Use its `quick/` workflow first; it is the most helpful version for
  understanding the example.

From the package root:

```sh
Rscript test-scripts/run_minimal_release_example.R
Rscript test-scripts/seven_node_example/quick/00_run_quick.R
Rscript test-scripts/seven_node_example/quick/01_plot_quick.R
```

Install the current checkout into an isolated R library before running examples.
The minimal and seven-node examples use generated inputs; other developer
workflows can require optional local data or calibration artifacts. Consult
their READMEs before running them.
