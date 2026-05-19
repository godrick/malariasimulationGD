# Careful Workflow

Four-stage version of the same seven-node example. It is useful when you want a
seasonally calibrated baseline, a promoted checkpoint, and paired release versus
no-release simulations.

For most package users, `quick/` is the more helpful workflow. Run this careful
version only if you need the extra calibration and checkpoint steps.

## Run

From the package root:

```sh
Rscript test-scripts/seven_node_example/careful/00_calibrate_init_eir.R
Rscript test-scripts/seven_node_example/careful/01_warmup_and_checkpoint.R
Rscript test-scripts/seven_node_example/careful/02_release_from_checkpoint.R
Rscript test-scripts/seven_node_example/careful/03_plot_careful.R
```

## Stages

| Stage | Output |
|---|---|
| `00_calibrate_init_eir.R` | Finds an `init_EIR` whose realised seasonal PfPR is near the target. |
| `01_warmup_and_checkpoint.R` | Builds a checkpoint library, scores stationarity, promotes one snapshot, and stores context. |
| `02_release_from_checkpoint.R` | Runs paired release and no-release arms from the promoted checkpoint. |
| `03_plot_careful.R` | Renders release spread, incidence, and PfPR diagnostics. |

The calibration starts from a non-seasonal analytical anchor, then refines using
the same seasonal simulator settings used by the example. The checkpoint
stationarity checks follow the production audit-cell pattern, with
example-local tolerance choices documented in `lib/baseline_library.R`.

## Outputs

All outputs are written to `test-scripts/seven_node_example/output/careful/`:

- `calibrated_init_eir.rds` and `calibration_log.csv`.
- `baseline_checkpoint_library.rds` and `baseline_checkpoint.rds`.
- `stationarity_battery_summary.csv`, `multi_seed_validation.rds`, and
  `multi_seed_summary.csv`.
- `context.rds`, which records movement, contact multipliers, release timing,
  `horizon_day`, and `readout_day`.
- Release and no-release timeseries, carrier frequency, incidence, prevalence,
  summary, and figure files.
