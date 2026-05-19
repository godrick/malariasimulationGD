# msimGD Share Examples

This directory is intentionally small and self-contained. It is ignored by
`R CMD build` through the package `.Rbuildignore`, so it can hold local example
outputs without changing the package build.

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
Rscript test-scripts/install_local.R
Rscript test-scripts/run_minimal_release_example.R
Rscript test-scripts/seven_node_example/quick/00_run_quick.R
Rscript test-scripts/seven_node_example/quick/01_plot_quick.R
```

These examples do not read calibration files, saved outputs, CHIRPS files, or
analysis artifacts from the original working repository.
