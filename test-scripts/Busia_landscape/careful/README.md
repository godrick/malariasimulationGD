# Busia 58-Node Careful Workflow

This folder contains the multi-stage Busia landscape workflow for the
`malariasimulationGD` simulation study. It uses the 58 Busia cluster
coordinates in `../Busia_data/busia_58_cluster_coordinates.csv`, the Busia
regional age structure in `../Busia_data/busia_age_structure.csv`, Busia
seasonality, and a homing-drive release design.

The careful workflow
separates baseline calibration, stationary checkpoint construction, release
simulation, and plotting so that the release and no-release arms begin from the
same promoted baseline state.

## Run Order

Run these scripts from the package root, in order:

```sh
Rscript test-scripts/Busia_landscape/careful/00_calibrate_init_eir.R
Rscript test-scripts/Busia_landscape/careful/01_warmup_and_checkpoint.R
Rscript test-scripts/Busia_landscape/careful/02_release_from_checkpoint.R
Rscript test-scripts/Busia_landscape/careful/03_plot_careful.R
```

By default, outputs are written to:

```text
test-scripts/Busia_landscape/output/careful/
```

You can redirect outputs for a trial run with:

```sh
BUSIA_CAREFUL_OUT_DIR=/private/tmp/busia_careful_test Rscript test-scripts/Busia_landscape/careful/00_calibrate_init_eir.R
```

Use the same `BUSIA_CAREFUL_OUT_DIR` value for all four stages if redirecting.

The default Busia population input is now a deterministic heterogeneous
58-node village vector with total population `44,875` residents and
cross-village sample SD approximately `60.3` residents. Because the total fixes
the arithmetic mean, the default mean is `44,875 / 58 = 773.7` residents per
village. You can override these inputs with `BUSIA_TOTAL_POPULATION` and
`BUSIA_POPULATION_SD`.

## Conceptual Stage Summary

| Script | Conceptual role | Main outputs |
|---|---|---|
| `00_calibrate_init_eir.R` | Finds the baseline `init_EIR` that gives the target seasonal mean PfPR2-10 in the Busia 58-node landscape. | `calibrated_init_eir.rds`, `calibration_log.csv`, `theta.rds` |
| `01_warmup_and_checkpoint.R` | Builds a stationary baseline checkpoint library, scores candidate snapshots, promotes one baseline state, and records workflow context. | `baseline_checkpoint_library.rds`, `baseline_checkpoint.rds`, stationarity and multi-seed summaries, `context.rds` |
| `02_release_from_checkpoint.R` | Starts from the promoted checkpoint and runs paired release and no-release simulations with identical baseline state and seed. | Release/no-release truth RDS files, timeseries CSVs, carrier frequency, prevalence, incidence, `summary.csv` |
| `03_plot_careful.R` | Reads the stage 02 summaries and renders the diagnostic figure for the release experiment. | `figure_careful.png`, `figure_careful.pdf` |

## `00_calibrate_init_eir.R`

This stage calibrates the starting annual entomological inoculation rate
(`init_EIR`) for the Busia baseline. The target is PfPR2-10 microscopy
prevalence, currently set by the shared Busia config as `TARGET_PFPR <- 0.39`.
Set `BUSIA_TARGET_PFPR` before stage 00 if you intentionally want a different
target.

Conceptually, it does two things:

1. It gets a fast analytical starting point using `malariaEquilibrium`.
2. It refines that starting point by running the actual Busia seasonal
   metapopulation simulator over a small EIR grid and measuring realised
   PfPR2-10.

The Busia landscape is built from the 58 coordinate rows, with latitude and
longitude converted to kilometer-scale `x`/`y` positions and a 58x58 distance
matrix. Human population and age structure are applied through the Busia
demography helpers. Mosquito infection-exposed/EIP staging is configured through
`theta$nEIP = 5L`.

The resulting calibrated `init_EIR`, target prevalence, realised prevalence,
Busia population metadata, and mosquito lifecycle settings are saved for later
stages.

## `01_warmup_and_checkpoint.R`

This stage turns the calibrated baseline into a reusable stationary baseline
state. It reads `calibrated_init_eir.rds`, rebuilds the same Busia 58-node setup,
then runs a long baseline burn-in and captures multiple candidate snapshots.

Conceptually, it answers: which saved baseline state is stable enough to use as
the common starting point for release and no-release comparisons?

It does this by:

1. Running a baseline burn-in.
2. Saving a library of snapshots spaced one year apart.
3. Re-running each snapshot for validation.
4. Scoring stationarity using year-over-year PfPR2-10 and mosquito abundance
   diagnostics.
5. Promoting the best snapshot.
6. Running a small multi-seed validation of the promoted snapshot.

The promoted checkpoint is saved as `baseline_checkpoint.rds`. The workflow
context saved in `context.rds` records the 58 Busia nodes, population vector,
movement settings, contact multipliers, release nodes, release timing,
seasonality, calibrated EIR, and mosquito lifecycle settings.

## `02_release_from_checkpoint.R`

This stage runs the actual paired intervention comparison. It reads the
promoted checkpoint and context from stage 01, then runs two arms:

- `release`: releases HH males into selected Busia nodes at the configured
  release day.
- `no_release`: runs the same baseline forward with no release.

Both arms use the same promoted checkpoint and same seed. This means differences
between the two arms are attributable to the release intervention rather than
different baseline warm-up histories.

The script writes raw truth outputs and analysis-ready CSVs:

- per-node release and no-release timeseries,
- drive carrier frequency by node and day,
- regional PfPR2-10 comparison,
- regional incidence comparison,
- one-row workflow summary.

It also performs sanity checks for the expected 58-node output shape, required
prevalence/incidence columns, genotype count attributes, and appearance of
H-carrier genotypes after release.

## `03_plot_careful.R`

This stage produces the final four-panel diagnostic figure. It reads `context.rds`
and the CSV outputs from stage 02, then plots:

- the Busia 58-node landscape with release nodes marked,
- drive carrier frequency over time by node,
- regional incidence for release versus no-release,
- regional PfPR2-10 for release versus no-release.

The figure is saved as both PNG and PDF. The script also prints the release and
no-release PfPR2-10 values closest to the configured readout day.

## Important Inputs

- `../Busia_data/busia_58_cluster_coordinates.csv`: 58 Busia cluster
  latitude/longitude coordinates.
- `../Busia_data/busia_age_structure.csv`: regional Busia age structure used to
  configure human demography.
- `../config/busia_landscape.R`: Busia landscape, demography, validation, and
  helper functions.
- `../config/seasonality.R`: Busia seasonal rainfall coefficients.
- `../config/movement.R`: mosquito movement assumptions.
- `../config/covariate.R`: per-node synthetic contact multiplier construction.
- `../config/homing_drive.R`: homing-drive inheritance cube.
- `../config/trial_design.R`: release timing, release size, and release-node
  selection settings.

## Notes

- The workflow is calibrated for a 58-node Busia landscape, not the original
  seven-node toy landscape.
- The current default total simulated population is `44,875` residents across
  the 58 Busia nodes, allocated heterogeneously with village-to-village SD
  approximately `60.3`. Set `BUSIA_TOTAL_POPULATION` before stage 00 if you
  intentionally want a different total population scale, and set
  `BUSIA_POPULATION_SD` if you want a different cross-village dispersion.
- Full stage 01 and stage 02 runs can be computationally expensive because they
  run many 58-node stochastic simulations.
