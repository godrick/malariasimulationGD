for (package in c("MGDrivE2", "MGDrivE", "deSolve", "malariasimulationGD")) {
  testthat::skip_if_not_installed(package)
}
mosquito_comparison_env <- new.env(parent = globalenv())
sys.source(file.path(developer_test_repo(), "test-scripts/mosquito_comparison/comparison.R"),
           envir = mosquito_comparison_env)

test_that("independent mosquito equilibria and actual derivatives agree", {
  e <- mosquito_comparison_env
  result <- e$mc_run_scenario(e$mc_scenarios()[[1L]], days = 3L)
  expect_lt(result$summary$initial_max_abs_difference, 1e-7)
  expect_lt(result$summary$K_abs_difference, 1e-7)
  expect_lt(result$summary$rhs_max_scaled_difference, 1e-8)
  expect_true(result$summary$overlap_pass)
  expect_equal(result$parameters$del, 1 / result$theta$qE)
  expect_equal(result$parameters$dl, 1 / result$theta$qL)
  expect_equal(result$parameters$dpl, 1 / result$theta$qP)
  female <- subset(result$totals, stage == "adult_female")
  expect_equal(female$mg2, rep(500, 4L), tolerance = 1e-7)
  expect_equal(female$mg4, rep(500, 4L), tolerance = 1e-7)
})

test_that("release comparison agrees with the public native standalone runner", {
  e <- mosquito_comparison_env
  for (id in c("mated_female", "mixed_sex", "egg_release", "larval_release", "pupal_release")) {
    scenarios <- e$mc_scenarios()
    scenario <- scenarios[[match(id, vapply(scenarios, `[[`, character(1), "id"))]]
    scenario$releases$timestep <- 2L
    result <- e$mc_run_scenario(scenario, days = 6L)
    public <- malariasimulationGD::run_mosquito_simulation(
      6L, result$parameters, sampler = "ode", releases = scenario$releases
    )
    expected <- public$counts[, c("timestep", "stage", "genotype", "count")]
    names(expected)[1L] <- "day"
    joined <- merge(result$genotypes, expected, by = c("day", "stage", "genotype"))
    expect_equal(nrow(joined), nrow(result$genotypes), info = id)
    expect_true(max(abs(joined$mg4 - joined$count)) < 1e-4, info = id)
    expect_true(result$summary$overlap_pass, info = id)
    expect_equal(result$releases$added_mg2, result$releases$requested)
    expect_equal(result$releases$added_mg4, result$releases$requested)
    expect_lt(max(result$releases$max_allocation_difference), 1e-6)
  }
})

test_that("direct mortality and habitat interventions have matched nontrivial responses", {
  e <- mosquito_comparison_env
  scenarios <- e$mc_scenarios()
  for (id in c("adult_mortality", "larvicide", "habitat", "fecundity")) {
    scenario <- scenarios[[match(id, vapply(scenarios, `[[`, character(1), "id"))]]
    result <- e$mc_run_scenario(scenario, days = 40L)
    female <- subset(result$totals, stage == "adult_female" & day == 40L)
    expect_true(result$summary$overlap_pass, info = id)
    expect_true(result$summary$rhs_max_scaled_difference < 1e-8, info = id)
    # Do not assume reducing egg input must reduce adults: crowding can
    # compensate or overcompensate in this Erlang larval model.
    expect_true(abs(female$mg2 - 500) > 1, info = id)
    expect_true(abs(female$mg4 - 500) > 1, info = id)
  }
})

test_that("the comparison detects a deliberately mismatched egg production rate", {
  e <- mosquito_comparison_env
  scenarios <- e$mc_scenarios()
  scenario <- scenarios[[match("control_beta_mismatch",
                               vapply(scenarios, `[[`, character(1), "id"))]]
  result <- e$mc_run_scenario(scenario, days = 10L)
  expect_false(result$summary$overlap_pass)
  expect_true(result$summary$expectation_pass)
  expect_gt(result$summary$max_tolerance_ratio, 1)
  expect_gt(result$summary$rhs_max_abs_difference, 1)
})
