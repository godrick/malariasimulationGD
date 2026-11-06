# Developer workflow tests

These tests validate developer scenarios, calibration helpers, and external
integration workflows. They are excluded from package builds with the rest of
`test-scripts/` and are not part of `testthat::test_local()` or `R CMD check`.
Package behavior must have small, self-contained regressions in
`tests/testthat/`; those tests must not source this workspace.

From the repository root:

```sh
Rscript test-scripts/tests/run_tests.R '^four-node-linear-example$'
Rscript test-scripts/tests/run_tests.R
```

The runner loads this checkout and uses a fixed initial seed. Optional workflows
skip explicitly when their fixtures or dependencies are unavailable. Running
every developer test is opt-in: external integrations can be more expensive
than the package suite.

The `malaria-sim-epi-*` tests retain their original external project layout.
Set `MSIMGD_INTEGRATION_SCRIPTS_ROOT` to the directory containing `msimGD/`,
`landscape-movement/`, and the other required integration helpers. Set
`MSIMGD_INTEGRATION_REPO_ROOT` to the parent project containing its `DESCRIPTION`
if it differs from the parent of that scripts directory. The defaults preserve
the old locations relative to the package checkout. No external data or packages
are downloaded by this runner.
