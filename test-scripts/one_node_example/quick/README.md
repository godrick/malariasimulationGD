# Quick One-Node Run

Run from the package root:

```sh
Rscript test-scripts/one_node_example/quick/00_run_quick.R
Rscript test-scripts/one_node_example/quick/01_plot_quick.R
```

`00_run_quick.R` runs one deterministic native mosquito backend simulation with
one homing-drive release. `01_plot_quick.R` reads the generated CSV files and
creates PNG figures for adult drive carrier frequency, malaria incidence, and
child prevalence.

Outputs are written to:

```text
test-scripts/one_node_example/output/quick/
```
