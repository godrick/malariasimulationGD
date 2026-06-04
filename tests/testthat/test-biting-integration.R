test_that('biting_process integrates mosquito effects and human infection', {
  population <- 4
  timestep <- 5
  parameters <- get_parameters(
    list(human_population = population)
  )

  renderer <- individual::Render$new(5)
  events <- mockery::mock()
  age <- c(20, 24, 5, 39) * 365
  variables <- list(birth = individual::DoubleVariable$new((-age + timestep)))
  lagged_foim <- LaggedValue$new(1, 1)
  lagged_eir <- LaggedValue$new(1, 1)
  models <- parameterise_mosquito_models(parameters, timestep)
  solvers <- parameterise_solvers(models, parameters)

  infection_outcome <- CompetingOutcome$new(
    targeted_process = function(timestep, target){
      infection_process_resolved_hazard(timestep, target, variables, renderer, parameters)
    },
    size = parameters$human_population
    )

  biting_process <- create_biting_process(
    renderer,
    solvers,
    models,
    variables,
    events,
    parameters,
    lagged_foim,
    lagged_eir,
    infection_outcome=infection_outcome
  )
  
  bitten <- list(
    bitten_humans = individual::Bitset$new(parameters$human_population),
    n_bites_per_person = numeric(0),
    transmission_multiplier = 1
  )
  bites_mock <- mockery::mock(bitten, cycle = T)
  infection_mock <- mockery::mock()

  mockery::stub(biting_process, 'simulate_bites', bites_mock)
  mockery::stub(biting_process, 'simulate_infection', infection_mock)
  biting_process(timestep)
  
  mockery::expect_args(
    bites_mock,
    1,
    renderer,
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    lagged_foim,
    lagged_eir,
    NULL,
    1,
    lagged_eir,
    NULL,
    NULL
  )

  mockery::expect_args(
    infection_mock,
    1,
    variables,
    events,
    bitten$bitten_humans,
    bitten$n_bites_per_person,
    age,
    parameters,
    timestep,
    renderer,
    infection_outcome,
    bitten$transmission_multiplier
  )
})

test_that('biting_process passes delayed human bites when mobility is enabled', {
  population <- 2
  timestep <- 5
  parameters <- get_parameters(
    list(
      human_population = population,
      human_mobility_enabled = TRUE,
      de = 1
    )
  )

  renderer <- individual::Render$new(5)
  events <- mockery::mock()
  age <- c(20, 24) * 365
  variables <- list(birth = individual::DoubleVariable$new((-age + timestep)))
  lagged_foim <- LaggedValue$new(1, 1)
  lagged_eir <- LaggedValue$new(1, 1)
  models <- list()
  solvers <- list()
  infection_outcome <- CompetingOutcome$new(
    targeted_process = function(timestep, target) invisible(NULL),
    size = parameters$human_population
  )
  human_exposure_lag_context <- new.env(parent = emptyenv())
  human_exposure_lag_context$weighted_active <- FALSE
  human_exposure_lag_context$buffers <- list(HumanExposureLagBuffer$new(
    max_lag = 2,
    default_exposure = rep(0, population)
  ))
  human_exposure_lag_context$buffers[[1L]]$save(timestep - parameters$de, c(2, 0))

  biting_process <- create_biting_process(
    renderer,
    solvers,
    models,
    variables,
    events,
    parameters,
    lagged_foim,
    lagged_eir,
    infection_outcome = infection_outcome,
    human_exposure_lag_context = human_exposure_lag_context
  )

  bitten <- list(
    bitten_humans = individual::Bitset$new(parameters$human_population),
    n_bites_per_person = numeric(0),
    transmission_multiplier = 99
  )
  bites_mock <- mockery::mock(bitten, cycle = TRUE)
  infection_mock <- mockery::mock()

  mockery::stub(biting_process, 'simulate_bites', bites_mock)
  mockery::stub(biting_process, 'simulate_infection', infection_mock)
  biting_process(timestep)

  mockery::expect_args(
    bites_mock,
    1,
    renderer,
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    lagged_foim,
    lagged_eir,
    NULL,
    1,
    lagged_eir,
    human_exposure_lag_context,
    NULL
  )

  infection_args <- mockery::mock_args(infection_mock)[[1L]]
  expect_identical(infection_args[[1L]], variables)
  expect_identical(infection_args[[3L]], bitten$bitten_humans)
  expect_identical(infection_args[[9L]], infection_outcome)
  expect_equal(infection_args[[10L]], bitten$transmission_multiplier)
  expect_null(infection_args$infection_exposure)
})

human_mobility_biting_integration_context <- function(parameters, variables, timestep, bite_rates, weighted_rates = NULL) {
  lagged_eir <- list(list(LaggedValue$new(2, 0)))
  context <- create_human_exposure_lag_context(
    parameters = list(parameters),
    variables = list(variables),
    solvers = list(NULL),
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_eir
  )
  context$buffers[[1L]]$save(
    timestep - parameters$de,
    exposure = bite_rates,
    weighted_exposure = weighted_rates,
    species_exposure = matrix(bite_rates, nrow = 1L),
    species_weighted_exposure = if (is.null(weighted_rates)) NULL else matrix(weighted_rates, nrow = 1L)
  )
  context
}

test_that('simulate_bites reconstructs delayed falciparum bitten_humans under mobility', {
  population <- 3L
  timestep <- 5L
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(list(
    human_population = population,
    human_mobility_enabled = TRUE,
    de = 1
  ))
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  age <- rep(20 * 365, population)
  context <- human_mobility_biting_integration_context(
    parameters,
    variables,
    timestep,
    bite_rates = c(0, 2, 1)
  )
  infectivity_context <- create_human_infectivity_lag_context(list(parameters), list(variables))

  pois_mock <- mockery::mock(2)
  sample_mock <- mockery::mock(c(2L, 3L))
  mockery::stub(simulate_bites, 'rpois', pois_mock)
  mockery::stub(simulate_bites, 'fast_weighted_sample', sample_mock)

  models <- parameterise_mosquito_models(parameters, timestep)
  solvers <- parameterise_solvers(models, parameters)
  bitten <- simulate_bites(
    renderer,
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    LaggedValue$new(2, 0),
    list(LaggedValue$new(2, 0)),
    human_exposure_lag_context = context,
    human_infectivity_lag_context = infectivity_context
  )

  expect_equal(bitten$bitten_humans$to_vector(), c(2L, 3L))
  expect_equal(bitten$n_bites_per_person, numeric(0))
  mockery::expect_args(pois_mock, 1, 1, 3)
  mockery::expect_args(sample_mock, 1, 2, c(0, 2, 1))
})

test_that('simulate_bites reconstructs delayed vivax n_bites_per_person under mobility', {
  population <- 3L
  timestep <- 5L
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(
    list(
      human_population = population,
      human_mobility_enabled = TRUE,
      de = 1
    ),
    parasite = "vivax"
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  age <- rep(20 * 365, population)
  context <- human_mobility_biting_integration_context(
    parameters,
    variables,
    timestep,
    bite_rates = c(0, 2, 1)
  )
  infectivity_context <- create_human_infectivity_lag_context(list(parameters), list(variables))

  pois_mock <- mockery::mock(3)
  sample_mock <- mockery::mock(c(2L, 2L, 3L))
  mockery::stub(simulate_bites, 'rpois', pois_mock)
  mockery::stub(simulate_bites, 'fast_weighted_sample', sample_mock)

  models <- parameterise_mosquito_models(parameters, timestep)
  solvers <- parameterise_solvers(models, parameters)
  bitten <- simulate_bites(
    renderer,
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    LaggedValue$new(2, 0),
    list(LaggedValue$new(2, 0)),
    human_exposure_lag_context = context,
    human_infectivity_lag_context = infectivity_context
  )

  expect_equal(bitten$bitten_humans$to_vector(), c(2L, 3L))
  expect_equal(bitten$n_bites_per_person, c(0L, 2L, 1L))
  mockery::expect_args(pois_mock, 1, 1, 3)
  mockery::expect_args(sample_mock, 1, 3, c(0, 2, 1))
})

test_that('simulate_bites integrates eir calculation and mosquito side effects', {
  population <- 4
  timestep <- 5
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(
    list(human_population = population,
         individual_mosquitoes = TRUE)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  infectivity <- c(.6, 0, .2, .3)
  age <- c(20, 24, 5, 39) * 365

  variables$zeta <- individual::DoubleVariable$new((c(.2, .3, .5, .9)))
  variables$infectivity <- individual::DoubleVariable$new(infectivity)
  variables$mosquito_state <- individual::CategoricalVariable$new(
    c('Sm', 'Pm', 'Im', 'NonExistent'),
    c(rep('Im', 10), rep('Sm', 15), rep('NonExistent', 75))
  )
  variables$species <- individual::CategoricalVariable$new(
    c('gamb'),
    rep('gamb', 100)
  )

  lambda_mock <- mockery::mock(c(.5, .5, .5, .5))
  mosquito_effects_mock <- mockery::mock()
  eqs_update <- mockery::mock()
  sample_mock <- mockery::mock(c(2, 3))
  pois_mock <- mockery::mock(2)

  mockery::stub(simulate_bites, 'biting_effects_individual', mosquito_effects_mock)
  mockery::stub(simulate_bites, 'rpois', pois_mock)
  mockery::stub(simulate_bites, 'fast_weighted_sample', sample_mock)
  mockery::stub(simulate_bites, 'effective_biting_rates', lambda_mock)
  mockery::stub(simulate_bites, 'aquatic_mosquito_model_update', eqs_update)
  models <- parameterise_mosquito_models(parameters, timestep)
  solvers <- parameterise_solvers(models, parameters)
  lagged_foim <- LaggedValue$new(12.5, .001)
  lagged_eir <- list(LaggedValue$new(12, 10))
  bitten <- simulate_bites(
    renderer,
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    lagged_foim,
    lagged_eir
  )

  expect_equal(bitten$bitten_humans$to_vector(), c(2, 3))

  f <- parameters$blood_meal_rates[[1]]

  effects_args <- mockery::mock_args(mosquito_effects_mock)

  expect_equal(effects_args[[1]][[1]], variables)
  expect_equal(effects_args[[1]][[3]], events)
  expect_equal(effects_args[[1]][[4]], 1)
  expect_equal(effects_args[[1]][[5]]$to_vector(), 11:25)
  expect_equal(effects_args[[1]][[6]]$to_vector(), c(1:10, 11:25))
  expect_equal(effects_args[[1]][[7]], parameters$mum)
  expect_equal(effects_args[[1]][[8]], parameters)
  expect_equal(effects_args[[1]][[9]], timestep)

  mockery::expect_args(eqs_update, 1, models[[1]]$.model, 25, f, parameters$mum)
  mockery::expect_args(
    pois_mock,
    1,
    1,
    10 * mean(unique_biting_rate(age, parameters))
  )
  mockery::expect_args(
    sample_mock,
    1,
    2,
    c(.5, .5, .5, .5)
  )
})
