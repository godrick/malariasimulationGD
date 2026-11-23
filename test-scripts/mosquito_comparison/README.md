# Deterministic mosquito comparison: MGDrivE2 and native malariasimulationGD

This developer experiment tests whether the two installed implementations give
overlapping mosquito-only trajectories after matching biological parameters,
all initial compartments, density dependence, releases and intervention clocks.
It does not modify package equations or defaults. It is excluded from package CI.

## Run

Use an installed `malariasimulationGD` built from the intended checkout. Prefer
an isolated temporary R library; do not use the default-library installer.
Additional workflow dependencies are `MGDrivE2` (tested with 2.1.1), `MGDrivE`
(tested with 1.6.2), `deSolve`, and `ggplot2` for plots. `testthat` is needed for
the developer regression tests. These are not new package dependencies.

From the repository root, selecting your library through `R_LIBS_USER` if needed:

```sh
Rscript test-scripts/mosquito_comparison/run_comparison.R /tmp/mosquito-comparison-results 120
Rscript -e 'testthat::test_dir("test-scripts/tests/testthat", filter = "mosquito-deterministic-comparison", reporter = "summary", stop_on_failure = TRUE)'
```

Choose a fresh results directory: the runner refuses to overwrite an existing
`summary.csv`. An omitted directory creates a temporary directory. The second
argument is the final day (minimum 90, default 120). The suite runs serially,
once per scenario at each of two solver precisions. No remote data are needed.

## What is matched

- One node, one species and 500 equilibrium adult females, with three genotypes.
- Mean egg/larval/pupal durations 4/3/6 days and 2/3/2 Erlang substages.
  **Correct mapping:** `del = 1/qE`, `dl = 1/qL`, `dpl = 1/qP`. The progression
  rates within the native engine are then `nE*qE`, `nL*qL`, `nP*qP`.
- Egg/larval/pupal mortality 0.05/0.15/0.05 per day; both adult sexes 0.09 per day;
  daily egg production 16; unmated-female mating rate 6 per day.
- Logistic density dependence; the actual K and complete equilibrium vectors
  are independently calculated by both packages and checked before each run.
- Exact state mapping by aquatic substage, genotype, sex, and mother's/mate's
  genotype. Each arm starts from the same complete vector, not just equal totals.
- Identical inheritance cubes and viable-offspring, fecundity, adult mortality,
  emergence and mating-preference modifiers. Wild-type traits remain neutral.
- Matching continuous-state hazard thresholds of `1e-12` in the comparison
  engines. This is an explicit numerical configuration, not a package default
  change. The native public wrapper uses a mating threshold of `1e-8`; developer
  tests separately compare its release trajectories with the comparison driver.

The MGDrivE2 arm evaluates that package's actual `spn_hazards()` and `spn_S()`
with `deSolve::ode(method = "lsoda")`. The native arm uses the installed C++
entomology engine and its Boost Dormand-Prince solver. No lifecycle equations
are copied into the comparison. Direct use of the engines allows explicit,
identical solver tolerances and mosquito-rate interventions without adding a
new public intervention API.

Both arms integrate one day, apply releases at the next boundary, and record
the state after the release. Intervention rates on days 30 through 59 apply
over `[30,60)`. Each rate-change boundary restarts the integrator with the full
current state. K is held at its original calibrated value unless K itself is
the intervention; it is never silently recalibrated to restore the population.
Seasonality uses identical daily piecewise-constant K in both arms.

## Scenarios

The 19 matched scenarios are:

1. Unperturbed equilibrium.
2. Initial removal of 70% of adults, leaving aquatic stages intact.
3. Single release of 250 aa adult males on day 20.
4. Five releases of 250 aa males, on days 20/30/40/50/60.
5. Release of 100 aa females already mated to AA males.
6. Mixed-sex release: 100 aa males followed by 100 aa females; females mate
   immediately using each arm's local male pool and preferences.
7. 2,000 aa eggs into egg substage 2.
8. 1,000 aa larvae into larval substage 2.
9. 250 aa pupae into pupal substage 1.
10. Initially 200 unmated AA females and no other mosquitoes; 200 AA males
    arrive on day 20.
11. Five releases of 500 aa sterile males. Females mated to aa males produce
    no viable offspring (`tau[, "aa", ] = 0`).
12. Male and female releases with aa fitness costs: adult lifespan modifier 0.8,
    fecundity modifier 0.6, female/male emergence modifiers 0.7/0.9, and positive,
    asymmetric mating preferences.
13. An illustrative homing-drive cube: heterozygotes transmit allele a with
    probability 0.95 in either sex. Release 100 Aa males on day 20. This is a
    synthetic inheritance test, not a calibrated gene-drive construct.
14. Adult mortality multiplied by 3 for both sexes over `[30,60)`.
15. Baseline larval mortality multiplied by 3 over `[30,60)`. Because mortality
    is `muL*(1+L/K)`, this scales the full larval death hazard.
16. Habitat K multiplied by 0.4 over `[30,60)`.
17. Egg production multiplied by 0.4 over `[30,60)`.
18. Repeated male releases plus the adult-mortality and habitat interventions.
19. Daily habitat K = K0 * (1 + 0.6*sin(2*pi*day/60)).

All single releases occur on day 20. A twentieth, deliberately mismatched
negative control increases only the native egg-production rate by 5%; the
comparison must reject overlap for this control.

The rate interventions represent mechanistic perturbations, not household
LLIN/IRS coverage, repellency or efficacy calculations. The standalone model
does not implement those human-derived intervention schedules. The experiment
also excludes infection, human dynamics, stochastic runs and spatial movement.

## Acceptance and diagnostics

The criterion is declared before running and applied to every compartment at
every daily boundary:

`abs(MG4 - MG2) <= 1e-4 + 1e-5 * max(abs(MG2), abs(MG4))`.

Standard solver tolerances are rtol/atol = `1e-8/1e-10`; every scenario is
repeated at `1e-10/1e-12`. Derivatives are also compared at the same state before
every day; matched scenarios must have maximum scaled derivative difference
below `1e-8`. Every release's injected count and compartment allocation are
audited. Solver convergence, nonnegative finite output and the deliberate
mismatch control provide independent checks against superficial agreement.

`summary.csv` and `summary_tight.csv` contain pass/fail results and maximum
errors. `compartment_trajectories.csv` retains all compartments; the stage and
genotype CSVs give aggregated trajectories. `release_audit.csv`,
`derivative_comparison.csv`, `solver_convergence.csv`, `comparison.rds` and
`provenance.txt` retain the evidence, configuration and installed versions.
PNG plots show adult females, all life stages, genotype trajectories, and the
error relative to the acceptance limit. The runner saves evidence before
raising an error if any expected outcome fails. Generated outputs are not
source files and should not be committed.

For the chosen illustrative parameters, lower egg production and sterile-male
releases increase adult females slightly in both implementations. This is
consistent with overcompensation from reduced larval crowding and must not be
read as a general intervention effect. Tests require a nontrivial, matched
response without assuming a monotonic decline in adults.
