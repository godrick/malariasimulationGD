# Code and model guide

Use this guide to locate behavior before changing it. Detailed parameter/output
definitions remain in the roxygen sources and installed help. Read
[CONTRIBUTING.md](CONTRIBUTING.md) for setup and validation, and
[AGENTS.md](AGENTS.md) for the scientific guardrails.

## Entry points and source map

| Responsibility | Main entry or source |
| --- | --- |
| Construct parameters and apply overrides | `get_parameters()` in [R/parameters.R](R/parameters.R) |
| Parameterize equilibrium and initialize individual state | `set_equilibrium()` in [R/compatibility.R](R/compatibility.R), variables in [R/variables.R](R/variables.R), and auxiliary human initialization in [R/stationary_human_initializer.R](R/stationary_human_initializer.R) |
| Run one population | `run_simulation()` and `run_resumable_simulation()` in [R/model.R](R/model.R) |
| Run interacting populations | `run_metapop_simulation()` in [R/model.R](R/model.R) |
| Build ordered daily processes | `create_processes()` in [R/processes.R](R/processes.R) |
| Scheduled individual events and listeners | [R/events.R](R/events.R) and the intervention modules |
| Mosquito backend selection and shared native state | [R/native_mosquito_backend.R](R/native_mosquito_backend.R) |
| Legacy mosquito models/solvers | [R/compartmental.R](R/compartmental.R), [R/mosquito_biology.R](R/mosquito_biology.R) |
| Native deterministic/stochastic engines | [src/det_engine_epi.cpp](src/det_engine_epi.cpp), [src/stoch_engine_epi.cpp](src/stoch_engine_epi.cpp), [src/engine_common.h](src/engine_common.h) |
| Human infection, progression, treatment | [R/human_infection.R](R/human_infection.R), [R/disease_progression.R](R/disease_progression.R) |
| Bites and exposure | [R/biting_process.R](R/biting_process.R) |
| Releases and consumed schedules | [R/mosquito_releases.R](R/mosquito_releases.R) |
| Overnight trips and human exposure/infectivity history | [R/human_mobility.R](R/human_mobility.R), [R/human_exposure_lag.R](R/human_exposure_lag.R), [R/human_infectivity_lag.R](R/human_infectivity_lag.R) |
| Rendering and age-specific outputs | [R/render.R](R/render.R), [R/output_parameters.R](R/output_parameters.R) |
| Mosquito lifecycle without the human infection model | `run_mosquito_simulation()` in [R/mosquito_standalone.R](R/mosquito_standalone.R) |

`run_simulation()` returns the data frame from `run_resumable_simulation()`.
The resumable function additionally returns state. Metapopulation output is
organized by node, with optional state when `return_state = TRUE`.

## Follow a run

`R/model.R` chooses the backend, validates configuration, initializes individual
variables and events, builds mosquito models/solvers and lag buffers, and creates
renderers. For a resumed run it restores the saved objects before continuing.
It then delegates the simulation loop to the `individual` package and assembles
data, optional genotype histories, and checkpoint state.

For one population, inspect `create_processes()` to trace the process order.
Its principal stages are:

1. Decay immunity and, for vivax, process hypnozoite batches.
2. Apply supported mosquito releases; legacy individual mosquitoes have their
   own emergence path.
3. When explicit mobility is enabled, update trips and the human infectivity
   context before biting.
4. Calculate biting/exposure and infection/progression hazards, then resolve
   the competing hazards.
5. Step the mosquito solvers at their designated position.
6. Apply configured PEV/PMC processes, render outputs, and process the remaining
   configured vector-control/progress processes.
7. Process human mortality at the end of the process list.

This is a source-reading map, not a replacement for the event scheduler's
semantics. Event listeners are initialized separately. Native metapopulations
use `native_interleave_metapop_processes()` so processes across nodes share the
backend in the required order. Do not reorder processes or move solver stepping
merely to simplify code: it can change exposure, releases, deaths, and RNG draws.

## Backend state and boundaries

The four flag combinations are listed in [the README](README.md#choose-the-mosquito-backend-deliberately).
In the native backend, `individual_mosquitoes = TRUE` selects count-based
stochastic dynamics; it does not create legacy individual mosquitoes. The human
infection model remains individual-based even with deterministic mosquito ODEs.

The native adapter builds indices for mosquito life stages, genotypes, species,
and nodes, then passes them to C++ through the generated Rcpp interface. Node
solvers can share one backend, so stepping/copying them as independent state
objects can be incorrect. Check native stage-layout and movement tests before
altering that ownership or stepping behavior.

`R/RcppExports.R` and `src/RcppExports.cpp` are generated interfaces. R-to-C++
index conventions must be verified at the boundary; do not assume an R index
can be passed to an arbitrary native array without conversion. Changing stage
counts also changes the state layout and can invalidate checkpoints.

## Movement and node ordering

Use a consistent node order for parameter lists, movement/mixing matrices,
contact surfaces, release targets, and output arrays. Validate both dimensions
and the meaning of each axis.

Explicit overnight human mobility currently requires:

- At least two nodes and the native metapopulation mosquito backend.
- A numeric, finite, nonnegative `human_move_probs` matrix with one row/column
  per node; each row sums to one. Rows are trip origins/home nodes, columns are
  destinations, and diagonal probabilities mean staying home when eligible.
- Identical values for any repeated non-NULL human movement matrices.
- Identity `export_mixing` and `import_mixing` matrices and zero `p_captured`
  matrices. Use `p_success = 0` for the no-capture setup.
- A positive `de` and the documented fixed/geometric trip-duration settings.
  `human_move_rates` is not supported by this explicit mobility path.

Travel changes location and exposure, not ownership of the resident's record.
Human infection, immunity, and intervention state stays associated with the
home record. Exposure and infectivity have separate lag buffers. A one-night
trip contributes away exposure for that step and returns on the following
step; it is not permanent migration.

Mosquito movement uses its own `move_probs`/`move_rates` configuration. Treat
that interface separately from human trip probabilities. The detailed
restrictions are checked by `validate_explicit_human_mobility_metapop()` and
related validators in `R/native_mosquito_backend.R`.

## Releases and outputs

Construct releases through `set_releases()` and inspect the consumed schedule
on output, together with the resulting sex/genotype histories. Preserve the
release day, interval, number, count, sex, genotype, and target-node semantics.
Early/late lifecycle effects are not established by schedule validation alone.

Output columns are not all in the same units. In particular, an annual
per-person initialization EIR is not directly comparable to a per-timestep
population-level rendered EIR. Read the relevant parameter and run-function
help before computing rates, prevalence, age ranges, or resident/present-node
summaries. Keep denominators and observation windows explicit in analyses.

## Continue a checkpoint

Using `parameters` from the README example and the selected installed library:

```r
set.seed(20261105)
first <- run_resumable_simulation(5L, parameters = parameters)
continued <- run_resumable_simulation(
  10L,
  parameters = parameters,
  initial_state = first$state,
  restore_random_state = TRUE
)
stopifnot(nrow(first$data) == 5L, nrow(continued$data) == 5L)
```

Here `10L` is the final simulation endpoint, so the second call returns days
6-10 rather than ten additional days. Preserve the complete state object;
rendered data frames alone are not checkpoints. Restore RNG state when a
continuation should follow its saved stochastic trajectory.

For strict continuous/resumed comparisons, follow
[test-resume.R](tests/testthat/test-resume.R): it initializes the shared lazy
correlation object before resetting the seed for the two arms. Explicit
mobility additionally saves per-human lag histories; combined resume coverage
is in [test-human-mobility-stage8.R](tests/testthat/test-human-mobility-stage8.R).
Changing parameters during resume is only supported where the relevant API and
tests establish it; changing native state dimensions is not a compatible edit.

## Choose regression coverage

| Behavior | Existing tests under `tests/testthat/` |
| --- | --- |
| Native samplers, stage layouts, mosquito movement | `test-native-mosquito-backend.R` |
| Legacy genotype/release compatibility | `test-mosquito-genotype-cube.R` |
| Mobility validation and trip/return timing | `test-human-mobility-stage1.R`, `test-human-mobility-stage2.R` |
| Per-human exposure, delayed infection, infectivity | `test-human-mobility-stage3.R` through `test-human-mobility-stage5.R` |
| Native handoff, releases with mobility, solver order | `test-human-mobility-stage6.R` |
| Documented mobility configuration and carried state | `test-human-mobility-stage7.R` |
| Combined mobility, diagnostics, and resume | `test-human-mobility-stage8.R` |
| Checkpoints and simulation entry points | `test-resume.R`, `test-simulation-e2e.R` |
| Vivax treatment recipient/drug alignment | `test-vivax-treatment-regression.R` |
| Standalone mosquito lifecycle | `test-mosquito-standalone.R` |

Historical stage numbers in test filenames describe the implementation sequence,
not supported feature versions. Use small synthetic fixtures to capture package
guarantees; keep long scenario studies in the developer workspace. See
[the baseline reports](tests/BASELINE-REPAIRS.md) for recorded validation and
known limits, rather than treating old pass counts as permanent expectations.
