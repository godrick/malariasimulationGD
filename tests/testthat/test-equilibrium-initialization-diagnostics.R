# The diagnostic helpers live outside the built package. Pure helper tests run
# locally; packaged checks skip explicitly when example sources are not bundled.
eqd_repo <- normalizePath(file.path(testthat::test_path(), "../.."), mustWork = TRUE)
eqd_root <- file.path(eqd_repo, "test-scripts/Godrick_calibration")
eqd_dir <- file.path(eqd_root, "diagnostics")
if (!file.exists(file.path(eqd_dir, "config.R"))) {
  testthat::skip("Godrick diagnostic sources are excluded from package builds.")
}
eqd_test_env <- new.env(parent = globalenv())
for (f in c("config.R", "simulation.R", "analysis.R", "legacy_equilibrium.R")) sys.source(file.path(eqd_dir, f), eqd_test_env)

test_that("diagnostic arguments and source fingerprint are validated", {
  e <- eqd_test_env
  expect_identical(e$eqd_options(character())$mode, "smoke")
  for (arg in c("--mode=full", "--cores=0", "--cores=3", "--cores=1.5", "--days=31", "--out=")) {
    expect_error(e$eqd_options(arg))
  }
  expect_equal(nrow(e$eqd_scenarios("screen")), 12L)
  expect_identical(e$eqd_scenarios("confirm")$id, "native_D")
  expect_identical(e$eqd_source_fingerprint(eqd_repo), e$eqd_source_fingerprint(eqd_repo))
})

test_that("constant exposure buffers resist attempted feedback and validate restart", {
  b <- eqd_test_env$eqd_constant_buffer(2)
  b$save(200, 1L)
  expect_identical(b$get(-100L), 2)
  expect_identical(b$get(9999L), 2)
  expect_identical(b$writes, 1L)
  expect_silent(b$restore_state(1L, b$save_state()))
  expect_error(b$restore_state(1L, list(value = 3)), "changed")
  expect_error(eqd_test_env$eqd_constant_buffer(-1))
})

test_that("independent RNG streams are reproducible without changing caller RNG", {
  set.seed(912L)
  old <- .Random.seed
  a <- eqd_test_env$eqd_streams(3L)
  expect_identical(.Random.seed, old)
  expect_identical(a, eqd_test_env$eqd_streams(3L))
  expect_false(identical(a[[1L]], a[[2L]]))
  expect_false(identical(a[[2L]], a[[3L]]))
})

test_that("stationarity bounds distinguish practical stability, drift and uncertainty", {
  e <- eqd_test_env
  expect_equal(e$eqd_margin("state_A__all", 0.4), 0.01)
  expect_equal(e$eqd_margin("ib_mean__all", 0), 0.1)
  expect_equal(e$eqd_margin("icm_mean__all", 0), 0.01)
  expect_equal(e$eqd_margin("mosquito_Im_gamb_count", 0), 1)
  expect_identical(e$eqd_interval(rep(0.4, 10), rep(0.405, 10), "state_A__all")$status, "pass")
  expect_identical(e$eqd_interval(rep(0.4, 10), rep(0.42, 10), "state_A__all")$status, "fail")
  expect_identical(e$eqd_interval(rep(0.4, 3), c(0.38, 0.4, 0.42), "state_A__all")$status, "inconclusive")
  expect_identical(e$eqd_interval(c(NA, 0.4, 0.4), rep(0.4, 3), "state_A__all")$status, "inconclusive")
  expect_true(e$eqd_interval(rep(0, 3), rep(0, 3), "state_Tr__all")$structural_zero)
  expect_false(e$eqd_interval(rep(0, 3), rep(0, 3), "ivm_q10__age_20_plus")$structural_zero)
  a <- expand.grid(replicate = 1:10, year = 27:31)
  a$days <- 365L
  a$state_A__all <- 0.4
  expect_identical(e$eqd_verdict(e$eqd_gate(a)$status), "pass")
  a$state_A__all[a$year == 31L] <- 0.43
  g <- e$eqd_gate(a)
  expect_true(all(g$status[g$phase == "late_burnin"] == "pass"))
  expect_true(all(g$status[g$phase == "observation"] == "fail"))
  expect_identical(e$eqd_verdict(e$eqd_gate(a[a$replicate != 1L, ])$status), "inconclusive")
  a$state_A__all <- 0.4
  a$days[1L] <- 30L
  expect_identical(e$eqd_verdict(e$eqd_gate(a)$status), "inconclusive")
})

test_that("Godrick scenarios preserve fixed biology and vary only stated controls", {
  skip_if_not(file.exists(file.path(eqd_root, "config/p20_scientific_defaults.R")), "Optional Godrick fixture is absent.")
  skip_if_not_installed("malariaEquilibrium")
  e <- eqd_test_env
  a <- e$eqd_parameters(eqd_repo, eqd_root, "A", population = 100L)
  d <- e$eqd_parameters(eqd_repo, eqd_root, "D", population = 100L)
  expect_false(a$enable_heterogeneity)
  expect_false(a$bednets)
  expect_true(d$enable_heterogeneity)
  expect_true(d$bednets)
  expect_equal(d$initial_bednet_coverage, 0.684)
  for (p in list(a, d)) {
    expect_true(p$native_mosquito_backend)
    expect_true(p$individual_mosquitoes)
    expect_true(p$native_total_M_equilibrium)
    expect_false(p$model_seasonality)
    expect_null(p$releases)
    expect_null(p$cube)
    expect_null(p$move_probs)
    expect_equal(p$b0, 0.55)
    expect_equal(p$dem, 11)
    expect_equal(p$native_mosquito_nEIP, 6L)
    expect_equal(e$eqd_ns("equilibrium_treatment_coverage")(p), 0)
  }
  legacy <- e$eqd_parameters(eqd_repo, eqd_root, "D", "legacy", population = 100L)
  expect_false(legacy$native_mosquito_backend)
  expect_false(legacy$individual_mosquitoes)
  deterministic <- e$eqd_parameters(eqd_repo, eqd_root, "D", "native_deterministic", population = 100L)
  expect_true(deterministic$native_mosquito_backend)
  expect_false(deterministic$individual_mosquitoes)
})

test_that("late summaries use five burn-in years and preserve missing observations", {
  a <- data.frame(scenario = "native_D", replicate = 1L, year = 1:31,
    ib_mean__all = 1:31, state_A__all = 0.4)
  a$state_A__all[28L] <- NA_real_
  x <- eqd_test_env$eqd_phase_summary(a)
  expect_equal(x$last_five_years[x$metric == "ib_mean__all"], mean(26:30))
  expect_true(is.na(x$last_five_years[x$metric == "state_A__all"]))
})

test_that("observing joint human state does not change production simulation", {
  skip_if_not(file.exists(file.path(eqd_root, "config/p20_scientific_defaults.R")), "Optional Godrick fixture is absent.")
  e <- eqd_test_env
  p <- e$eqd_parameters(eqd_repo, eqd_root, "D", population = 100L)
  rng <- e$eqd_streams(1L)[[1L]]
  assign(".Random.seed", rng, .GlobalEnv)
  production <- suppressWarnings(e$eqd_ns("run_resumable_simulation")(31L, p))
  observed <- suppressWarnings(e$eqd_run(p, 31L, rng))
  expect_equal(observed$data, production$data, tolerance = 0)
  expect_true(all(rowSums(observed$data[paste0(e$eqd_states, "_count")]) == 100L))
  expect_equal(nrow(observed$initial_profile), 100L)
  audit <- e$eqd_initial_audit(observed, p)
  expect_equal(audit$continuous_vector_EIR, p$init_EIR, tolerance = 1e-8)
  expect_true(all(audit[grep("initialization_error", names(audit))] == 0))
  expect_true(all(e$eqd_immunities %in% names(observed$initial_profile)))
})

test_that("fixed exposure drives real infection and remains independent of mosquito feedback", {
  skip_if_not(file.exists(file.path(eqd_root, "config/p20_scientific_defaults.R")), "Optional Godrick fixture is absent.")
  e <- eqd_test_env
  p <- e$eqd_parameters(eqd_repo, eqd_root, "A", population = 100L)
  rng <- e$eqd_streams(1L, 17L)[[1L]]
  fixed <- suppressWarnings(e$eqd_run(p, 31L, rng, "fixed"))
  expect_true(all(fixed$data$EIR_gamb == p$init_EIR * p$human_population / 365))
  expect_gt(sum(fixed$data$n_bitten), 0)
  expect_gt(sum(fixed$buffer_audit$ignored_writes), 0)
  expect_gt(sum(fixed$buffer_audit$reads), 0)
  expect_equal(fixed$buffer_audit$weighted, fixed$buffer_audit$raw)
  expect_true(all(rowSums(fixed$data[paste0(e$eqd_states, "_count")]) == 100L))
  annual <- e$eqd_annual(fixed, p)
  expect_equal(annual$days, 31L)
  expect_equal(annual$state_S__all, mean(fixed$data$S_count / 100))
  expect_true("ib_q90__age_2_10" %in% names(annual))
  fixed$details$ib_q90__age_2_10[1L] <- NA_real_
  expect_true(is.na(e$eqd_annual(fixed, p)$ib_q90__age_2_10))
})

test_that("ICL human equilibrium matches its fixed-exposure ODE", {
  skip_if_not_installed("MGDrivE2")
  r <- eqd_test_env$eqd_icl_reference()
  expect_true(all(r$pass))
  expect_lt(max(r$residual), 1e-10)
})

test_that("legacy 5000-human experiment derives abundance and disables interventions", {
  skip_if_not(file.exists(file.path(eqd_root, "config/p20_scientific_defaults.R")), "Optional Godrick fixture is absent.")
  e <- eqd_test_env
  p <- e$eqd_legacy_parameters(eqd_repo, eqd_root)
  expect_identical(p$human_population, 5000L)
  expect_false(p$native_mosquito_backend)
  expect_false(p$individual_mosquitoes)
  expect_true(p$enable_heterogeneity)
  expect_equal(p$sigma_squared, 1.67)
  for (name in c("bednets", "mda", "smc", "pev", "tbv")) expect_false(p[[name]])
  expect_length(p$clinical_treatment_coverages, 0L)
  expect_equal(e$eqd_ns("equilibrium_treatment_coverage")(p), 0)
  expect_equal(p$total_M, e$eqd_ns("equilibrium_total_M")(p, p$init_EIR), tolerance = 1e-10)
  smaller <- e$eqd_parameters(eqd_repo, eqd_root, case = "B", backend = "legacy", population = 2500L)
  expect_equal(p$total_M, 2 * smaller$total_M, tolerance = 1e-8)
})

test_that("50-year legacy gate detects slower drift and incomplete windows", {
  e <- eqd_test_env
  a <- expand.grid(replicate = 1:3, year = 1:50)
  a$days <- 365L
  a$state_A__all <- 0.4
  a$ib_mean__all <- 50
  expect_identical(e$eqd_verdict(e$eqd_legacy_gate(a)$status), "pass")
  a$state_A__all[a$year %in% 46:50] <- 0.42
  g <- e$eqd_legacy_gate(a)
  expect_true(all(g$status[g$metric == "state_A__all" & g$phase == "years_36_40_to_46_50"] == "fail"))
  expect_identical(e$eqd_verdict(e$eqd_legacy_gate(a[a$replicate != 1L, ])$status), "inconclusive")
  a$state_A__all <- 0.4
  a$days[a$year == 37L & a$replicate == 2L] <- 30L
  g <- e$eqd_legacy_gate(a)
  expect_true(all(g$status[g$phase == "years_36_40_to_46_50"] == "inconclusive"))
  a$state_A__all[a$year == 48L & a$replicate == 1L] <- NA_real_
  expect_true(is.na(e$eqd_legacy_windows(a)$years_46_50[1L]))
})

test_that("legacy replicate count is validated and all 50 streams enter assessment", {
  e <- eqd_test_env
  expect_identical(e$eqd_legacy_options(character())$reps, 3L)
  cfg <- e$eqd_legacy_options(c("--mode=run", "--reps=50"))
  expect_identical(cfg$reps, 50L)
  expect_identical(cfg$mode, "run")
  for (arg in c("--reps=0", "--reps=-1", "--reps=1.5", "--reps=Inf", "--reps=abc", "--reps=")) expect_error(e$eqd_legacy_options(arg))
  expect_error(e$eqd_legacy_options(c("--reps=3", "--reps=50")))
  a <- expand.grid(replicate = 1:50, year = 1:50)
  a$days <- 365L
  a$state_A__all <- 0.4
  a$state_A__all[a$year >= 46L & a$replicate <= 3L] <- 0.42
  g <- e$eqd_legacy_gate(a, reps = 50L)
  expect_true(all(g$n == 50L))
  expect_identical(e$eqd_verdict(g$status), "pass")
  expect_identical(e$eqd_verdict(e$eqd_legacy_gate(a, reps = 3L)$status), "fail")
  expect_identical(e$eqd_verdict(e$eqd_legacy_gate(a[a$replicate != 50L, ], reps = 50L)$status), "inconclusive")
})

test_that("ensemble intervals and cumulative means reflect replicate variation", {
  e <- eqd_test_env
  x <- e$eqd_legacy_mean_interval(c(1, 2, 3))
  expect_equal(unname(x[c("n", "mean", "se")]), c(3, 2, 1 / sqrt(3)))
  expect_lt(x[["ci_low"]], 2)
  expect_gt(x[["ci_high"]], 2)
  expect_true(is.na(e$eqd_legacy_mean_interval(1)[["se"]]))
  expect_true(is.na(e$eqd_legacy_mean_interval(c(NA, 1, 2))[["mean"]]))
  a <- expand.grid(replicate = 1:3, year = 1:50)
  metrics <- c(paste0("state_", e$eqd_states, "__all"), paste0(e$eqd_immunities, "_mean__all"),
    "prevalence_lm__182_5475", "mosquito_adult_female")
  for (metric in metrics) a[[metric]] <- a$replicate
  c <- e$eqd_legacy_convergence(a)
  expect_equal(c$mean[c$metric == "state_S__all"], c(1, 1.5, 2))
  expect_true(all(c$final_ensemble_mean == 2))
  m <- e$eqd_legacy_ensemble(a)
  expect_true(all(m$n == 3 & m$mean == 2))
})
