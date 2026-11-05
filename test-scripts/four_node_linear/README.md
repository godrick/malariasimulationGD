# Four sites on a line: A–B–C–D

Run one stochastic replicate with **1,825 release-free burn-in days**, followed
by **1,825 observed days (five 365-day years)**. Release **10,000 male HH
mosquitoes at each of A and C on observed day 183** (six months, rounded up;
absolute simulation day **2,008**). The full continuous simulation is 3,650 days.
The burn-in and observation share all human and mosquito state and the same RNG
stream; the model is not reinitialized at observation day 1.

Sites are at 0, 2, 4 and 6 km, with **10,000 residents each** to reduce stochastic
variation relative to the original 1,500-resident example. All baseline
parameter lists are identical, including contact multipliers (explicitly 1),
biting rates, human exposure-heterogeneity distribution, bednet coverage,
initial EIR, carrying capacity, and mosquito life history. A check rejects any
difference other than release schedules. Realized populations and incidence
can differ because individuals are stochastic and the movement geometry has
edge and interior sites. There are no human trips, cross-site human biting, or
border captures. Mosquitoes can move between every pair of sites.

## Run

Use an installed `malariasimulationGD` representing the current working tree,
plus `MGDrivE` and `ggplot2`. From the repository root:

```sh
Rscript test-scripts/four_node_linear/00_run.R
```

The script runs the model, checks outputs, and creates all four requested plots
automatically in `test-scripts/four_node_linear/output/plots/`:

- `carrier_B_D.png` / `.pdf`: GD carrier frequency at B and D.
- `carrier_A_C.png` / `.pdf`: GD carrier frequency at A and C.
- `incidence_B_D.png` / `.pdf`: child infection incidence at B and D.
- `incidence_A_C.png` / `.pdf`: child infection incidence at A and C.

`overview.png` combines the four plots. Time zero is the end of the burn-in;
only the five observed years are plotted. The dashed line marks the release day.
Regenerate plots without rerunning the simulation:

```sh
Rscript test-scripts/four_node_linear/01_plot.R
```

For an isolated installation, avoiding changes to the default R library:

```sh
mkdir -p /tmp/four_node_library
R CMD INSTALL --library=/tmp/four_node_library .
Rscript test-scripts/four_node_linear/00_run.R --library=/tmp/four_node_library
```

The run refuses to overwrite existing output. Select a fresh destination with
`--out=PATH` (also supported by `01_plot.R`). A short end-to-end check uses
30 burn-in days plus 30 observed days, observed release day 10, and 100 residents
per site (the reduced population is for smoke testing only):

```sh
Rscript test-scripts/four_node_linear/00_run.R --smoke --out=/tmp/four_node_smoke
```

Optional run overrides: `--days=N` (observed days), `--burnin-days=N`,
`--release-day=N` (relative to observation), `--population=N`, and `--seed=N`.
`config.R` holds the remaining scenario settings. The default seed
is 20261029; one replicate uses one core.

## Scientific settings and interpretation

Settings were taken from the local Godrick calibration files
`config/p20_scientific_defaults.R`, `config/p20_preset.R`,
`R/epidemiology_parameters.R`, and `diagnostics/config.R`. They are declared
here explicitly; this example does not require that optional local bundle,
external rainfall data, or existing calibration artifacts.

- Native stochastic count-based mosquitoes: `native_mosquito_backend = TRUE`,
  `individual_mosquitoes = TRUE`, tau step 0.1 day.
- Analytic initialization uses `set_equilibrium(..., native_total_M = TRUE)`
  and Godrick's accepted `init_EIR = 18.4730728587`. This is an initialization
  parameter borrowed from that calibration, **not a recalibration of this
  four-site scenario**. Five release-free years of coupled human/mosquito
  burn-in precede observation, followed by six observed months before release.
  `baseline_diagnostics.csv` reports annual burn-in incidence and incidence
  during the observed pre-release period. These diagnostics describe remaining
  drift and stochastic variation; a fixed burn-in is not a formal proof of
  equilibrium.
- Constant carrying capacity: no rainfall seasonality or time-dependent
  rescaling; all sites share identical baseline parameters and population.
  Native genotype-aware capacity is calculated once from this equilibrium,
  audited for equality across sites and time, and recorded in `nodes.csv`.
  Adult mosquito abundance can still vary with genetics and stochastic events.
- Mosquito biology: mean egg/larva/pupa durations 6/21/2 days, 2/3/2 stages;
  mortality 0.05/0.15/0.05 per day, adult mortality 0.132 per day; fecundity
  16, mating rate 6 per day, mean EIP 11 days with six stages.
- Human baseline: `b0 = 0.55`, `d1 = 0.160527`, zero routine treatment,
  lognormal exposure heterogeneity variance 1.67 with five groups.
- Stationary PBO bednet use 68.4%; `dn0 = 0.437`, `rn = 0.555`, `rnm = 0.24`,
  efficacy and retention time constants `1e12` days, initial net time 1,
  timestep-1 distribution skipped. There are no subsequent campaigns.
- Homing drive: cleavage 1, homing 0.95, resistant repair 0.05 in both sexes.
  Allele fecundities W=1, H=0.8, R=1, B=0.5 multiply within a genotype.
  Mosquito-to-human transmission is 0 for H carriers and 0.55 otherwise;
  mosquito infection susceptibility is 1 for every genotype.
- Movement uses `P[i,j] ∝ exp(-distance[i,j] / 2 km)` for `i != j`, with zero
  diagonal and origin rows normalized to one. All sites have equal destination
  attractiveness. The **2 km value is a decay length**, not a target mean move
  distance. Godrick's target mean of 2 km is the lower feasible bound on this
  geometry, so a finite distance-decay kernel cannot achieve it while keeping
  longer moves positive. The achieved mean distance is saved in `context.rds`.
- Godrick's `p_move = 0.02` is preserved using the rate
  `0.02 * 0.132 / (1 - 0.02) ≈ 0.002694/day`. It represents the probability
  of moving before death under competing hazards, not a 2% daily move rate.

GD carrier frequency is the number of **adult males plus females carrying H**
(WH, HH, HR, HB), divided by all adults at the site. It is a carrier fraction,
not an allele frequency; it is missing if there are no adults.

The primary incidence endpoint is **new infections**, not clinical cases or
first-ever infections. Rendered age bounds are 182–3649 days, matching the
Godrick calibration's six-month-to-under-ten-year band. Nonoverlapping 30-day
windows show `365 * sum(new infections) / sum(daily children in band)`, in
infections per child-year. The final shorter window uses its actual observed
person-time. Clinical counts and rates are also saved for reference.

## Saved data and checks

`simulation.rds` preserves all 3,650 raw daily outputs per node and genotype attributes;
`parameters.rds`, `context.rds`, `sessionInfo.txt`, `nodes.csv`,
`distance_km.csv`, and `movement_probabilities.csv` preserve inputs.
`carrier_frequency.csv`, `incidence_daily.csv`, `incidence_windows.csv`,
`release_schedule.csv`, and `summary.csv` contain derived outputs and runtime.
The main carrier/incidence tables contain only the observed period, with
`timestep` relative to observation and `simulation_day` where applicable.
The consumed release schedule records observed day 183 and simulation day
2,008. `burnin_incidence_daily.csv` and `burnin_carrier_frequency.csv` retain
the excluded burn-in; `baseline_diagnostics.csv` records release-free baseline
incidence and mean adult abundance for each site and burn-in year.

Checks cover identical baseline site parameters, row-normalized movement and
node order, fixed equal capacity,
complete daily age-specific outputs, nonnegative genotype counts, carrier
frequency bounds, no H carriers before release, the simulator's consumed
schedule, and the day-of-release HH male pulse at A and C. Package equations,
defaults, and interfaces are unchanged. Generated output is ignored by Git.
