vivax_treatment_regression_parameters <- function() {
  parameters <- get_parameters(list(human_population = 8L), parasite = "vivax")
  parameters <- set_drugs(parameters, list(CQ_PQ_params_vivax, CQ_TQ_params_vivax))
  parameters <- set_clinical_treatment(parameters, 1L, 1L, 0.5)
  set_clinical_treatment(parameters, 2L, 1L, 0.5)
}

test_that("liver treatment keeps drug IDs aligned after blood-stage failures", {
  set.seed(20261105)
  parameters <- vivax_treatment_regression_parameters()
  target <- individual::Bitset$new(8L)$insert(c(2L, 4L, 5L, 7L))
  draws <- mockery::mock(c(1L, 3L, 4L), c(2L, 3L))
  local_mocked_bindings(bernoulli_multi_p = draws)

  treated <- calculate_successful_treatments(
    parameters, target, c(2L, 1L, 2L, 1L), 5L, NullRender$new(5L), ""
  )

  expect_equal(treated$successfully_treated$to_vector(), c(2L, 5L, 7L))
  expect_equal(treated$successfully_treated_hypnozoites$to_vector(), c(5L, 7L))
  expect_equal(treated$drugs_hypnozoites, c(2L, 1L))
  mockery::expect_args(draws, 2L, parameters$drug_hypnozoite_efficacy[c(2L, 2L, 1L)])
  expect_length(mockery::mock_args(draws), 2L)
})

test_that("liver treatment remains aligned after early treatment failure", {
  set.seed(20261105)
  parameters <- vivax_treatment_regression_parameters()
  parameters$antimalarial_resistance <- TRUE
  target <- individual::Bitset$new(8L)$insert(c(2L, 4L, 5L, 7L))
  draws <- mockery::mock(c(1L, 3L, 4L), c(2L, 3L), integer(0), 1L)
  local_mocked_bindings(
    bernoulli_multi_p = draws,
    get_antimalarial_resistance_parameters = function(...) list(
      artemisinin_resistance_proportion = rep(1, 3L),
      early_treatment_failure_probability = c(1, 0, 0),
      slow_parasite_clearance_probability = rep(0, 3L),
      dt_slow_parasite_clearance = rep(parameters$dt, 3L)
    )
  )

  treated <- calculate_successful_treatments(
    parameters, target, c(2L, 1L, 2L, 1L), 5L, NullRender$new(5L), ""
  )

  expect_equal(treated$successfully_treated$to_vector(), c(5L, 7L))
  expect_equal(treated$successfully_treated_hypnozoites$to_vector(), 5L)
  expect_equal(treated$drugs_hypnozoites, 2L)
  mockery::expect_args(draws, 4L, parameters$drug_hypnozoite_efficacy[c(2L, 1L)])
  expect_length(mockery::mock_args(draws), 4L)
})

test_that("clinical vivax treatment updates liver drugs for the correct people", {
  set.seed(20261105)
  parameters <- vivax_treatment_regression_parameters()
  target <- individual::Bitset$new(8L)$insert(c(2L, 4L, 5L, 7L))
  variables <- lapply(c("state", "infectivity", "progression_rates", "drug",
                        "drug_time", "hypnozoites", "ls_drug", "ls_drug_time"),
                      function(name) list(queue_update = mockery::mock()))
  names(variables) <- c("state", "infectivity", "progression_rates", "drug",
                       "drug_time", "hypnozoites", "ls_drug", "ls_drug_time")
  local_mocked_bindings(
    sample_bitset = function(...) target,
    bernoulli_multi_p = mockery::mock(c(1L, 3L, 4L), c(2L, 3L))
  )
  calculate <- calculate_treated
  mockery::stub(calculate, "sample.int", function(...) c(2L, 1L, 2L, 1L))

  treated <- calculate(variables, target, parameters, 5L, NullRender$new(5L))

  expect_equal(treated$to_vector(), c(2L, 5L, 7L))
  expect_bitset_update(variables$hypnozoites$queue_update, 0, c(5L, 7L))
  expect_bitset_update(variables$ls_drug$queue_update, c(2L, 1L), c(5L, 7L))
  expect_bitset_update(variables$ls_drug_time$queue_update, 5L, c(5L, 7L))
})

test_that("blood-only vivax treatment does not require a liver-treatment result", {
  set.seed(20261105)
  parameters <- get_parameters(list(human_population = 8L), parasite = "vivax")
  parameters <- set_drugs(parameters, list(CQ_params_vivax))
  parameters <- set_clinical_treatment(parameters, 1L, 1L, 1)
  target <- individual::Bitset$new(8L)$insert(c(2L, 7L))
  variables <- lapply(c("state", "infectivity", "progression_rates", "drug", "drug_time"),
                      function(name) list(queue_update = mockery::mock()))
  names(variables) <- c("state", "infectivity", "progression_rates", "drug", "drug_time")
  local_mocked_bindings(
    sample_bitset = function(...) target,
    bernoulli_multi_p = function(...) c(1L, 2L)
  )

  treated <- calculate_treated(variables, target, parameters, 5L, NullRender$new(5L))

  expect_equal(treated$to_vector(), c(2L, 7L))
})

test_that("empty blood or liver treatment successes retain empty recipient mappings", {
  parameters <- vivax_treatment_regression_parameters()
  target <- individual::Bitset$new(8L)$insert(c(2L, 7L))
  for (blood_successes in list(integer(0), c(1L, 2L))) {
    set.seed(20261105)
    draws <- mockery::mock(blood_successes, integer(0))
    calculate <- calculate_successful_treatments
    mockery::stub(calculate, "bernoulli_multi_p", draws)

    treated <- calculate(parameters, target, c(2L, 1L), 5L, NullRender$new(5L), "")

    expect_equal(treated$successfully_treated_hypnozoites$size(), 0)
    expect_length(treated$drugs_hypnozoites, 0L)
    expect_length(mockery::mock_args(draws), 2L)
  }
})
