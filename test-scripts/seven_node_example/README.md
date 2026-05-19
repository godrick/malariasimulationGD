# 7-node example

Self-contained seven-village metapopulation example for `malariasimulationGD`.
It includes mosquito movement, seasonality, a synthetic per-village contact
multiplier, and a homing-drive release in three villages.

Everything needed to run the example is under this folder. It does not read
saved outputs, calibration files, CHIRPS files, or private framework code from
the original research repository.

## Workflows

Use `quick/` first. It is the most helpful version of this example because it
shows the full seven-node release setup without the calibration and checkpoint
machinery. The `careful/` workflow is optional and mainly useful when you
specifically need a calibrated checkpoint and paired release/no-release arms.

| Workflow | What it does | When to use it |
|---|---|---|
| `quick/` | One release simulation from analytic equilibrium with a hardcoded `init_EIR`. | Recommended path for most users. |
| `careful/` | Calibrates `init_EIR` under the seasonal simulator, builds/promotes a baseline checkpoint, then runs paired release and no-release arms. | Optional; use only when you need the checkpoint workflow. |

Both workflows share the same `config/` and `lib/` files. You do not need to run
both.

## Install

From the package root:

```sh
Rscript test-scripts/install_local.R
```

This installs the package and example runtime dependencies into
your default R library.

## Run

```sh
# Quick path
Rscript test-scripts/seven_node_example/quick/00_run_quick.R
Rscript test-scripts/seven_node_example/quick/01_plot_quick.R

# Careful path
Rscript test-scripts/seven_node_example/careful/00_calibrate_init_eir.R
Rscript test-scripts/seven_node_example/careful/01_warmup_and_checkpoint.R
Rscript test-scripts/seven_node_example/careful/02_release_from_checkpoint.R
Rscript test-scripts/seven_node_example/careful/03_plot_careful.R
```

Outputs land in `output/quick/` or `output/careful/`. Those folders are
gitignored.

## What It Demonstrates

- A seven-node landscape with village-specific human population sizes.
- Between-village mosquito movement from `(mu, p_move, distance)`.
- Seasonality through the package's Fourier-rainfall mechanism; coefficients
  are already stored in `config/seasonality.R`.
- A synthetic contact multiplier applied through
  `malariasimulationGD::apply_node_contact_surface()`. This is retained for
  compatibility with the example metadata. New truth-generation work should
  provide household/slot contact values through `human_slot_contact_multiplier`
  directly.
- A homing-drive cube from `MGDrivE::cubeHomingDrive()` with explicit homing
  rates in `config/homing_drive.R`.
- A release scheduled in three villages and spread through the metapopulation.
- Incidence and PfPR summaries plotted at the readout day
  `release_day + horizon_day`.

## Layout

```text
seven_node_example/
  config/    example-specific landscape, movement, seasonality, contact, drive,
             and release settings
  lib/       reusable helper functions kept local to this example
  quick/     shortest release-simulation workflow
  careful/   calibrated checkpoint workflow with paired release/no-release arms
  output/    generated files, ignored by git
```
