# Development baseline — 2026-11-05

## Purpose and provenance

This is the stage 1 baseline before introducing formatting, linting, GitHub
workflows, or refactoring. It records observed behavior and existing problems;
passing these checks does not establish scientific validity for every scenario.

- Package: `malariasimulationGD` 2.0.2.
- Source revision: `1f4e8bd58dd10a60b264909a859aea0afd409a79`.
- The working tree was clean before baseline work began.
- Checks used disposable copies of tracked files. Ignored local files and
  external integration fixtures were deliberately not included, matching what
  another contributor can obtain from the repository.
- Builds, libraries, logs, serialized test results, and example outputs stayed
  in a temporary directory. No dependency was installed or upgraded.
- This report is the only repository change. Implementation, tests, generated
  interfaces, documentation, defaults, equations, and RNG behavior are unchanged.
- After validation, all 358 pre-existing tracked files matched their recorded
  SHA-256 hashes.

## Environment and dependencies

The observed environment was R 4.6.1 on macOS Tahoe 26.6.2, ARM64, with the R
BLAS/LAPACK implementation. This is a single-platform baseline, not a supported
R-version or operating-system matrix.

| Dependency/tool | Installed version |
| --- | --- |
| individual | 0.1.19 |
| malariaEquilibrium | 1.0.1 |
| malariaEquilibriumVivax | 1.0.1 |
| Rcpp | 1.1.2 |
| RcppArmadillo | 15.6.0.1 |
| statmod | 1.5.2 |
| MASS | 7.3.65 |
| dqrng | 0.4.1 |
| sitmo | 2.0.2 |
| BH | 1.90.0.1 |
| R6 | 2.6.1 |
| progress | 1.2.3 |
| testthat | 3.3.2 |
| mockery | 0.4.5 |
| pkgload | 1.5.3 |
| MGDrivE | 1.6.2 |
| knitr | 1.52 |
| rmarkdown | 2.32 |
| digest (used optionally, currently undeclared) | 0.6.39 |

All declared `Imports` and `LinkingTo` packages were available. Missing declared
`Suggests` were `reshape2`, `DiagrammeR`, and `cowplot`. External integration
tests also encountered missing `customMGDrive2`, `nloptr`, and fixture files.
These are coverage limitations, not successful checks of those integrations.

`src/Makevars` requests C++14. This R installation reports its C++14 configuration
as defunct and compiled the package using Apple clang 21.0.0 with
`-std=gnu++20`. The package check explicitly notes that the C++14 request is
ignored. Successful compilation here does **not** verify C++14 compatibility.

## Executed checks

| Check | Result | Elapsed time |
| --- | --- | --- |
| Focused native-backend source tests | 154 passing expectations; no failures, errors, warnings, or skips | 4.456 s |
| Complete source test suite | 1,875 passing expectations; no failures or errors; one warning; nine reported skips | 169.072 s |
| Source package build without manuals/vignette outputs | Succeeded | Not separately timed |
| Strict package check | Stopped at dependency checking: one error for missing suggested packages | Not separately timed |
| Installed-package check with missing suggestions allowed and vignette execution disabled | One error, three warnings, two notes; details below | 84.898 s |
| Minimal release example using the freshly built package | Succeeded; package location verified against the temporary check library | 0.332 s |

The source-test harness set `set.seed(20261105)`, retained existing per-test
seeds, and disabled early stopping so all available tests could run. `CI` was
unset: the large falciparum and vivax equilibrium tests ran and passed. The
source result contained 372 test-result rows across 53 files; eight additional
files skipped at file scope. The ninth skip was an empty test.

The installed-package check used `CI=true` and
`_R_CHECK_FORCE_SUGGESTS_=false`. Existing `skip_on_ci()` calls consequently
skipped eleven equilibrium tests that had run in the source suite. It reported
1,875 passes, one test error, two test warnings, and 21 skips. Source and installed
test totals should not be treated as identical coverage.

Times describe single local runs with compilation already cached for the
focused source test. Some checks overlapped. These are workflow-duration
observations, not performance benchmarks or speedup claims.

## Existing failures and warnings

1. **Installed-package test failure.**
   `test-human-mobility-stage6.R:375` reads
   `../../R/native_mosquito_backend.R` to inspect constructor calls. That source
   path is available in a source checkout but not in the installed package test
   context. The resulting file-open error fails package checking, although the
   same test passes in the source suite. This is a test-portability problem;
   it does not demonstrate a change in mobility behavior. A future fix must
   preserve the constraint being tested rather than simply remove the assertion.

2. **Dependency declaration warning.**
   `R/stationary_human_initializer.R` optionally calls `digest`, but `digest`
   is not declared in `DESCRIPTION`. The checker also reports declared imports
   without namespace imports: `BH`, `dqrng`, `R6`, `RcppArmadillo`, and `sitmo`.
   Review their R and compiled-code roles before changing dependency fields.

3. **R code analysis note.**
   Reported unresolved names are `setNames`, `drug`, `parameter_draws_pf`, and
   `parameter_draws_pv`. These require individual review; the report does not
   assume that every static-analysis finding is a runtime bug.

4. **Source-test warning.**
   The age-output test in `test-output.R:4` reports
   `coercing argument of type 'list' to logical`. The installed suite also warns
   when the source-inspection test cannot open its input file.

5. **Compiler diagnostics.**
   `src/solver.h:14` emits a diagnostic-pop warning without a matching push.
   `src/stoch_engine_epi.cpp:1077` contains an unused `nState` variable.
   The C++14/C++20 discrepancy above is a separate compatibility limitation.

6. **Vignette limitations.**
   Building without vignette outputs leaves no `inst/doc`, producing vignette
   warnings. An exploratory check using only `--no-build-vignettes` also
   executed vignette code. Its long metapopulation vignette run was interrupted;
   the resulting incomplete-run error is not classified as a package defect.
   The bounded check explicitly adds `--no-vignettes`. A complete vignette build
   and uninterrupted vignette execution remain unverified.

7. **Dependency repository access.**
   Check output includes unavailable CRAN/Bioconductor index warnings in the
   restricted environment. Installed dependencies were used; online dependency
   resolution and a fresh-machine installation were not verified.

The strict-check error for missing suggested packages is environmental. The
installed-package source-path failure is independently actionable even when
missing suggestions are allowed. Neither was fixed during baseline work.

## Small regression scenarios for subsequent changes

These existing tests provide useful starting points before the full suite.
Times and expectation counts below come from the complete source run.

| Test file (under `tests/testthat/`) | Main coverage | Passing expectations | Time |
| --- | --- | --- | --- |
| `test-native-mosquito-backend.R` | Native deterministic/stochastic paths, one/two nodes, initialization, mosquito movement, releases | 154 | 3.391 s |
| `test-mosquito-genotype-cube.R` | Legacy individual-mosquito genotype compatibility, releases, trait validation | 101 | 2.044 s |
| `test-human-mobility-stage1.R` | Backend restrictions, movement dimensions, mixing constraints | 18 | 0.221 s |
| `test-human-mobility-stage2.R` | Identity/nonidentity movement, trip duration, return timing, resume | 63 | 0.357 s |
| `test-human-mobility-stage8.R` | Identity RNG equivalence, destination exposure, diagnostic allocation, checkpoint state | 62 | 1.626 s |
| `test-resume.R` | Legacy resume equivalence and adding interventions | 60 | 7.402 s |
| `test-simulation-e2e.R` | Simulation entry points, metapopulation integration, rendering/state | 24 | 0.361 s |
| `test-mosquito-standalone.R` | Native lifecycle, ODE/tau samplers, releases, checkpoint equivalence | 21 | 0.277 s |
| `test-four-node-linear-example.R` | Burn-in/release indexing and output extraction | 25 | 0.086 s |

Keep backend flags explicit in new scenarios. In the native tests,
`individual_mosquitoes = FALSE` selects deterministic dynamics and `TRUE`
selects stochastic dynamics. This does not make the native backend
interchangeable with the legacy individual-mosquito backend.

### Recorded minimal release output

`test-scripts/run_minimal_release_example.R` ran against the package installed
from this snapshot, with its existing seed `20260514`, 50 humans, initial
`total_M = 200`, `native_mosquito_backend = TRUE`, and
`individual_mosquitoes = FALSE`. There was one male `aa` release on day 20.

| Output | Observed value |
| --- | --- |
| Timesteps | 60 |
| Consumed release day | 20 |
| Consumed released genotype | `aa` |
| Consumed release count | 200 |
| Adult `aa` on release day | 175.268199016194 |
| First reported adult `Aa` day | 20 |
| Final total adult females | 200.000017654131 |

These are observations of this configuration, not new biological acceptance
criteria. Preserve existing scientific tolerances and compare configurations
and dependency versions before interpreting differences. Full trajectory
reference files were not added to the repository.

## Repeating the baseline

Run from a disposable source copy so compilation and simulation outputs do not
alter a contributor's working tree. Keep raw logs outside the repository.

```sh
Rscript -e 'set.seed(20261105); testthat::set_max_fails(Inf); testthat::test_local(".", reporter = "summary")'
R CMD build --no-manual --no-build-vignettes .
CI=true _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --no-build-vignettes --no-vignettes malariasimulationGD_2.0.2.tar.gz
```

For a strict dependency check, omit `_R_CHECK_FORCE_SUGGESTS_=false`. For the
example, install this same source build into an isolated temporary R library
and put that library first in `R_LIBS`; do not use an unrelated preinstalled
package or modify the default library.

For a focused test, use the file stem as the filter, for example:

```sh
Rscript -e 'set.seed(20261105); testthat::test_local(".", filter = "^native-mosquito-backend$")'
```

## Follow-up boundaries

- Keep the installed-package test repair separate from formatting work.
- Resolve dependency and static-analysis findings individually before requiring
  a clean package check in GitHub branch protection.
- Establish a supported R/compiler/platform matrix; the declared R minimum and
  C++14 compatibility have not been validated here.
- Complete optional integration fixtures and vignette validation separately.
- The seven-node workflow, Busia workflows, calibration runs, profiling,
  replicated scientific validation, coverage, and linting were not run as part
  of this baseline. Small tests and the minimal example establish the initial
  comparison points without regenerating landscape outputs.
