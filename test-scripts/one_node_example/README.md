# One-Node Native Gene-Drive Example

This is a minimal one-node `malariasimulationGD` example for quick native
mosquito backend sanity checks and visualization. It runs a single-node malaria
simulation with a homing-drive release and genotype fitness costs on the drive
cube.

This example intentionally has:

- one node only;
- no landscape;
- no mosquito movement;
- no human movement;
- no seasonality;
- no vector control.

It is not a spatial calibration workflow and does not use landscape geometry,
distance matrices, checkpoint machinery, or calibration helpers.

## Run

From the package root:

```sh
Rscript test-scripts/one_node_example/quick/00_run_quick.R
Rscript test-scripts/one_node_example/quick/01_plot_quick.R
```

The run script loads the working-tree package with `pkgload::load_all()` when
available, then builds a one-node native backend parameter set using
`set_equilibrium(..., native_total_M = TRUE)`.

## Scenario

- Human population: `1000`
- Initial annual EIR: `10`
- Native mosquito backend: `native_mosquito_backend = TRUE`
- Stochastic finite native mosquitoes: `individual_mosquitoes = TRUE`
- Mosquito movement: `move_probs = matrix(1, 1, 1)`, `move_rates = 0`
- Human movement: `human_mobility_enabled = FALSE`
- Seasonality: `model_seasonality = FALSE`
- Release: `5000` male `HH` mosquitoes on day `90`

For `individual_mosquitoes = TRUE`, `native_total_M` solves the continuous adult
mosquito abundance needed to match the requested annual EIR. The initialized
mosquito state is then integerized because the individual backend represents
finite mosquitoes. Therefore realised initial EIR may differ slightly from
requested `init_EIR`. The deterministic backend preserves the continuous target
exactly; the stochastic individual backend preserves integer mosquito counts and
reports the realised EIR.

## Fitness Costs

The homing-drive cube is built with `MGDrivE::cubeHomingDrive()` and then
applies multiplicative genotype fitness costs to `cube$s` from allele-level
values:

```r
allele_s <- c(W = 1, H = 0.8, B = 0.5, R = 1)
```

This gives:

```text
WW   WH   WR   WB   HH   HR   HB   RR   RB   BB
1.00 0.80 1.00 0.50 0.64 0.80 0.40 1.00 0.50 0.25
```

Drive carriers are adult mosquitoes with any genotype containing `H`.

## Outputs

Files are written under `test-scripts/one_node_example/output/quick/`:

- `timeseries.csv`: malaria and mosquito time series, including `n_infections`,
  `EIR_gamb`, and `pfpr_2_10_lm`.
- `carrier_frequency.csv`: adult drive carrier frequency and carrier counts,
  with female and male splits.
- `release_schedule.csv`: release schedule consumed by the simulator.
- `summary.csv`: one-row run summary.
- `context.rds`: scenario metadata, fitness costs, and genotype labels.

The plotting script writes:

- `adult_drive_carrier_frequency.png`
- `incidence_over_time.png`
- `child_prevalence_over_time.png`
- `one_node_quick_summary.png` when `patchwork` is available

## Layout

```text
one_node_example/
  config/homing_drive.R  one-node homing-drive cube and fitness-cost helper
  quick/                 runnable quick simulation and plotting scripts
  output/quick/          generated CSVs, RDS context, and figures
```
