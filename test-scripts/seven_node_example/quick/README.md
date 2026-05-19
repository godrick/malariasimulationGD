# Quick Workflow

Single release simulation plus one plotting script. This path skips calibration,
warmup, checkpoints, and counterfactual arms. It starts from
`set_equilibrium()` at a hardcoded `init_EIR`, releases on day 90, and runs to
the readout day `release_day + horizon_day` (`90 + 365 = 455` days).

From the package root:

```sh
Rscript test-scripts/seven_node_example/quick/00_run_quick.R
Rscript test-scripts/seven_node_example/quick/01_plot_quick.R
```

Outputs are written to `test-scripts/seven_node_example/output/quick/`:

- `timeseries.csv`: per-day, per-node simulation summary.
- `carrier_frequency.csv`: per-day, per-node adult drive-carrier frequency.
- `release_schedule.csv`: schedule consumed by the simulator.
- `nodes.csv` and `context.rds`: geometry, movement, timing, and contact inputs.
- `summary.csv`: one-row run summary with `release_day`, `horizon_day`, and
  `readout_day`.
- `figure_quick.png` and `figure_quick.pdf`: diagnostic plots.

Use `careful/` instead if you need a seasonally calibrated baseline, a promoted
checkpoint, or a paired no-release counterfactual.
