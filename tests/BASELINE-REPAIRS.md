# Baseline repairs — 2026-11-05

This records the follow-up to [the original baseline](BASELINE.md), which is
retained unchanged. Validation used disposable copies of the working tree on
R 4.6.1/macOS ARM64. No default R library, dependency installation, branch,
commit, or existing scientific output was changed.

## Package and developer boundaries

`test-scripts/` is the developer workspace for scenarios, calibration,
scientific validity checks, examples, and generated outputs. None of it is
required by the distributed package or its regression suite.

Nine workflow-specific test files moved from `tests/testthat/` to
`test-scripts/tests/testthat/`:

- `test-four-node-linear-example.R`
- `test-equilibrium-initialization-diagnostics.R`
- The seven `test-malaria-sim-epi-*.R` files

Their assertions were retained; path resolution and skip messages were adapted
to the developer location. The separate runner loads the current checkout:

```sh
Rscript test-scripts/tests/run_tests.R '^four-node-linear-example$'
Rscript test-scripts/tests/run_tests.R
```

Optional external project paths can be configured with
`MSIMGD_INTEGRATION_SCRIPTS_ROOT` and `MSIMGD_INTEGRATION_REPO_ROOT`; the runner's
README describes the expected layout. External fixtures remain optional and
are not fetched automatically.

`AGENTS.md` and the developer README now state this boundary explicitly.
`.Rbuildignore` continues to exclude the entire developer directory and also
excludes these baseline reports. Archive inspection confirmed that developer
scripts, developer tests, and baseline reports do not ship.

## Repairs

### Installed-package mobility regression

The mobility test now captures actual calls to both native engine constructors
instead of opening an R source file. It verifies that C++ human movement is
disabled and the human movement rates are zero for both samplers. The real
constructors still run. No backend or movement constraints changed.

### Dependencies and static analysis

- Declared optional `digest` in `Suggests`; retained and tested the existing
  fallback when it is unavailable. The availability-dependent hash algorithms
  are unchanged.
- Removed `reshape2`, `DiagrammeR`, and `cowplot` from `Suggests` after auditing
  the tracked package implementation, tests, and shipped documentation for use.
- Retained header dependencies `BH`, `dqrng`, `sitmo`, and `RcppArmadillo` in
  `LinkingTo`, removing their unused R-level `Imports` entries. The existing
  `dqrng >= 0.4` requirement is preserved in `LinkingTo`.
- Retained the used runtime dependency `R6` and generated its namespace import
  from roxygen source. Roxygen2 8.1.0 regenerated `NAMESPACE`; the corresponding
  tool-version metadata was updated. Rcpp interface files remained unchanged.
- Qualified `stats::setNames()` and package dataset references explicitly.

The installed-package checks for dependency declarations, loading with stated
dependencies, and possible R code problems now all pass.

### Output warning and test determinism

`set_epi_outputs()` now uses typed logical results from `vapply()` instead of
coercing a list inside `any()`. Its existing acceptance/rejection rules were
preserved. Successful and rejected configurations are tested without warnings;
stricter validation of partially incompatible output requests is a separate
API decision.

Expanded testing exposed an order-dependent MDA fixture: its assertions assumed
all treatments succeeded, but it had omitted the efficacy mock already used in
the neighboring tests. The missing mock and an explicit seed were added. The
original assertions, biological parameters, and production sampling were not
weakened or changed.

### Vivax treatment correctness

New tests reproduced failures in the original implementation before repair.
The affected path had three related problems:

1. Appending a Bitset via `c()` did not preserve the expected named result field.
2. Liver-stage success indices were applied to the original target even though
   their efficacy draws corresponded to blood-stage treatment successes.
3. Clinical liver-drug updates referenced undefined `drug`; blood-only vivax
   treatment also attempted to read an absent liver-treatment result.

The result now retains the liver-recipient Bitset and its matching drug vector.
For example, original recipients `2, 4, 5, 7` receiving drugs `2, 1, 2, 1`, with
blood-stage successes `2, 5, 7` and liver successes at positions `2, 3`, update
people `5, 7` with drugs `2, 1`. Early treatment failures, no successes, and
blood-only treatment are covered as well.

This preserves the existing assumption that these liver efficacy draws apply
to the blood-stage-success subset; it does not redesign drug efficacy biology.
The sampling calls, their order, and their probability vectors within treatment
were retained. Affected vivax clinical/MDA treatment outcomes are intentionally
corrected, and subsequent trajectories can consequently differ. Other model
defaults, equations, release conventions, mobility constraints, checkpoint
formats, and RNG implementations were not changed.

### Compiler diagnostics and documentation

Added the matching diagnostic push in `src/solver.h` and removed the unused
`nState` local in `stoch_epi_step_native()`. No C++ signature, equation, or
integration setting changed. `src/Makevars` still requests C++14.

All 15 existing user-facing vignettes built successfully with their declared
chunk settings, using the available Pandoc 3.10 executable. No software was
installed to enable that build. Intentionally unevaluated vignette chunks are
not claimed as executed scientific validation.

## Validation results

| Validation | Result |
| --- | --- |
| Focused mobility/output/dataset and vivax regressions | Passed |
| Expanded infection, MDA, cache, and vivax tests | Passed after correcting the test fixture and repair |
| Complete source suite, seed `20261105`, CI unset | 1,892 passing expectations; zero failures, errors, or warnings; one pre-existing empty-test skip; 178.973 s |
| Installed-package check, strict suggested dependencies, `CI=true` | Zero errors or warnings; two compatibility notes; 82.889 s |
| Installed test suite within that check | 1,920 passes; zero failures or warnings; 12 skips |
| Developer test runner on the tracked-files snapshot | 25 passes; eight optional integration/fixture skips |
| Package build including all 15 HTML vignettes | Passed; 303.971 s |
| Archive audit | Developer files/reports absent; all 15 built HTML vignettes present |
| Minimal release example against the fresh temporary installation | Passed; every summary field matches the original baseline, numeric absolute tolerance `1e-10` |

The installed suite's 12 skips consist of eleven existing `skip_on_ci()`
equilibrium tests and the existing empty listener test. The eleven equilibrium
tests ran in the complete source suite. Source/installed reporters and coverage
differ, so raw expectation totals are not a measure of a numerical change.

The minimal example retained its original seed `20260514`, native backend,
`individual_mosquitoes = FALSE`, 50 humans, 60 timesteps, and day-20 male release.
Its release schedule and numerical summary remain unchanged. Full trajectories
or new golden-data artifacts were not added.

Times are single local validation observations, not benchmark comparisons.
Logs, compiled objects, package archives, installed libraries, rendered
vignettes, and simulation outputs remain in temporary directories.

## Remaining compatibility decisions and limits

The two package-check notes are:

- R 4.6 ignores the obsolete C++14 request and uses C++20 in this environment.
  C++14 compatibility has not been demonstrated here; the requested standard
  was preserved as instructed.
- The declared R minimum of 2.10 predates versioned `LinkingTo` support, which
  requires R >= 3.0.2. This note becomes visible when the preserved minimum
  dqrng version is placed with its header dependency. The real supported R
  minimum must be established from code/dependency requirements and a tested
  version matrix rather than guessed to silence the note.

Choose the supported R/compiler matrix in a separate compatibility change.
No cross-platform or older-R claims follow from these local results.

The pre-existing empty listener test remains a coverage gap. The external
developer integrations with missing packages/fixtures were not validated.
Busia workflows, large calibration/replicated scientific studies, profiling,
PDF manual production, coverage, and linting were not run in this repair phase.

## Repeating validation

Run these commands in a disposable checkout/copy. Configure an available Pandoc
for rmarkdown and use an isolated library for installations and examples.

```sh
Rscript -e 'set.seed(20261105); testthat::set_max_fails(Inf); testthat::test_local(".", reporter = "summary")'
Rscript test-scripts/tests/run_tests.R '^four-node-linear-example$'
R CMD build --no-manual .
CI=true R CMD check --no-manual --no-build-vignettes --no-vignettes malariasimulationGD_2.0.2.tar.gz
```

Here `R CMD build` performs the vignette build. The final check avoids repeating
those simulations while still checking the shipped vignette artifacts. Unlike
the initial baseline workaround, `_R_CHECK_FORCE_SUGGESTS_=false` is not needed
in this environment after dependency cleanup.
