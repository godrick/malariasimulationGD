human_mobility_stage4_variables <- function(
  population,
  states = rep("S", population),
  hypnozoites = rep(0L, population)
) {
  list(
    birth = individual::DoubleVariable$new(rep(-1000, population)),
    state = individual::CategoricalVariable$new(
      categories = c("S", "A", "U", "D", "Tr"),
      initial_values = states
    ),
    drug = individual::IntegerVariable$new(rep(0L, population)),
    drug_time = individual::DoubleVariable$new(rep(-1, population)),
    last_eff_pev_timestep = individual::IntegerVariable$new(rep(-1L, population)),
    pev_profile = individual::IntegerVariable$new(rep(-1L, population)),
    ib = individual::DoubleVariable$new(rep(0, population)),
    last_boosted_ib = individual::DoubleVariable$new(rep(-1, population)),
    hypnozoites = individual::IntegerVariable$new(hypnozoites)
  )
}

human_mobility_stage4_outcome <- function(population) {
  CompetingOutcome$new(
    targeted_process = function(timestep, target) invisible(NULL),
    size = population
  )
}

human_mobility_stage4_context <- function(exposure, weighted_exposure = NULL, max_lag = 3) {
  context <- new.env(parent = emptyenv())
  context$weighted_active <- !is.null(weighted_exposure)
  context$buffers <- list(HumanExposureLagBuffer$new(
    max_lag = max_lag,
    default_exposure = rep(0, length(exposure)),
    default_weighted_exposure = if (is.null(weighted_exposure)) NULL else rep(0, length(exposure))
  ))
  context$buffers[[1L]]$save(3, exposure, weighted_exposure)
  context
}

test_that("delayed destination exposure drives home-node infection rates", {
  population <- 2L
  timestep <- 5L
  parameters <- get_parameters(list(
    human_population = population,
    human_mobility_enabled = TRUE,
    de = 2
  ))
  variables <- human_mobility_stage4_variables(population)
  outcome <- human_mobility_stage4_outcome(population)
  context <- human_mobility_stage4_context(c(0, 4))
  input <- human_exposure_lag_get_infection_input(context, 1L, timestep, parameters)

  simulate_infection(
    variables = variables,
    events = create_events(parameters),
    bitten_humans = individual::Bitset$new(population),
    n_bites_per_person = numeric(0),
    age = rep(1000, population),
    parameters = parameters,
    timestep = timestep,
    renderer = NullRender$new(timestep),
    infection_outcome = outcome,
    transmission_multiplier = input$transmission_multiplier,
    infection_exposure = input$infection_exposure
  )

  expect_equal(outcome$target$to_vector(), 2L)
  expect_equal(
    outcome$rates,
    input$infection_exposure[[2L]] * blood_immunity(variables$ib$get_values(2L), parameters)
  )
})

test_that("current high exposure does not affect infection before the delay", {
  population <- 1L
  timestep <- 5L
  parameters <- get_parameters(list(
    human_population = population,
    human_mobility_enabled = TRUE,
    de = 2
  ))
  variables <- human_mobility_stage4_variables(population)
  outcome <- human_mobility_stage4_outcome(population)

  context <- new.env(parent = emptyenv())
  context$weighted_active <- FALSE
  context$buffers <- list(HumanExposureLagBuffer$new(max_lag = 3, default_exposure = 0))
  context$buffers[[1L]]$save(timestep, 10)
  input <- human_exposure_lag_get_infection_input(context, 1L, timestep, parameters)

  expect_equal(input$infection_exposure, 0)

  simulate_infection(
    variables = variables,
    events = create_events(parameters),
    bitten_humans = individual::Bitset$new(population),
    n_bites_per_person = numeric(0),
    age = 1000,
    parameters = parameters,
    timestep = timestep,
    renderer = NullRender$new(timestep),
    infection_outcome = outcome,
    transmission_multiplier = input$transmission_multiplier,
    infection_exposure = input$infection_exposure
  )

  expect_null(outcome$rates)
})

test_that("weighted delayed exposure supplies multiplier with zero-denominator fallback", {
  parameters <- get_parameters(list(
    human_population = 2L,
    human_mobility_enabled = TRUE,
    de = 2
  ))
  context <- human_mobility_stage4_context(c(0, 4), c(9, 2))

  input <- human_exposure_lag_get_infection_input(context, 1L, 5L, parameters)

  expect_equal(input$infection_exposure, c(0, 4))
  expect_equal(input$weighted_exposure, c(9, 2))
  expect_equal(input$transmission_multiplier, c(1, 0.5))
})

test_that("infection from mobility exposure updates the home-node human record", {
  population <- 1L
  timestep <- 5L
  parameters <- get_parameters(list(
    human_population = population,
    human_mobility_enabled = TRUE,
    de = 2
  ))
  variables <- human_mobility_stage4_variables(population)
  infection_outcome <- CompetingOutcome$new(
    targeted_process = function(timestep, target) {
      variables$state$queue_update("A", target)
    },
    size = population
  )

  mockery::stub(calculate_falciparum_infections, "blood_immunity", mockery::mock(1))

  simulate_infection(
    variables = variables,
    events = create_events(parameters),
    bitten_humans = individual::Bitset$new(population),
    n_bites_per_person = numeric(0),
    age = 1000,
    parameters = parameters,
    timestep = timestep,
    renderer = NullRender$new(timestep),
    infection_outcome = infection_outcome,
    transmission_multiplier = 1,
    infection_exposure = 10
  )
  CompetingHazard$new(
    outcomes = list(infection_outcome),
    size = population,
    rng = function(n) rep(0.5, n)
  )$resolve(timestep)
  variables$state$.update()

  expect_equal(variables$state$get_values(), "A")
})

test_that("vivax relapse rates are preserved while bite exposure uses delayed intensity", {
  population <- 2L
  timestep <- 5L
  parameters <- get_parameters(
    list(
      human_population = population,
      human_mobility_enabled = TRUE,
      de = 2,
      b = 0.25,
      f = 0.1
    ),
    parasite = "vivax"
  )
  variables <- human_mobility_stage4_variables(population, hypnozoites = c(3L, 0L))
  outcome <- human_mobility_stage4_outcome(population)

  calculate_vivax_infections(
    variables = variables,
    bitten_humans = individual::Bitset$new(population),
    n_bites_per_person = numeric(0),
    parameters = parameters,
    renderer = NullRender$new(timestep),
    timestep = timestep,
    infection_outcome = outcome,
    transmission_multiplier = 1,
    infection_exposure = c(0, 2)
  )

  expect_equal(outcome$target$to_vector(), c(1L, 2L))
  expect_equal(outcome$rates, c(0.3, 0.5), tolerance = 1e-12)
  expect_equal(outcome$args$relative_rates, 1)
})

test_that("explicit mobility rejects nonpositive de only when enabled", {
  disabled <- list(
    get_parameters(list(de = 0)),
    get_parameters(list(de = 0))
  )
  expect_silent(human_mobility_validate_parameters(disabled, n_nodes = 2L))

  enabled <- disabled
  enabled[[1L]]$human_mobility_enabled <- TRUE
  enabled[[1L]]$human_move_probs <- diag(2)

  expect_error(
    human_mobility_validate_parameters(enabled, n_nodes = 2L),
    "greater than 0"
  )
})

test_that("diag human_move_probs preserves home-node delayed infection input", {
  parameters <- list(
    get_parameters(list(human_population = 2L, human_mobility_enabled = TRUE, de = 1)),
    get_parameters(list(human_population = 3L, human_mobility_enabled = TRUE, de = 1))
  )
  for (i in seq_along(parameters)) {
    parameters[[i]]$human_move_probs <- diag(2)
    parameters[[i]]$human_mobility_node_index <- as.integer(i)
  }
  variables <- lapply(parameters, create_variables)
  lagged_eir <- lapply(c(0, 0), function(default) list(LaggedValue$new(3, default)))
  context <- create_human_exposure_lag_context(
    parameters = parameters,
    variables = variables,
    solvers = list(NULL, NULL),
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_eir
  )

  human_exposure_lag_record_node(context, node_index = 1L, timestep = 3L, exposure = 5)
  human_exposure_lag_record_node(context, node_index = 2L, timestep = 3L, exposure = 7)

  input_1 <- human_exposure_lag_get_infection_input(context, 1L, 4L, parameters[[1L]])
  input_2 <- human_exposure_lag_get_infection_input(context, 2L, 4L, parameters[[2L]])

  expect_equal(input_1$infection_exposure, rep(5, 2))
  expect_equal(input_2$infection_exposure, rep(7, 3))
})
