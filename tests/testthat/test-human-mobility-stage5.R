human_mobility_stage5_params <- function(
  delay_gam = 1,
  population_1 = 2L,
  population_2 = 2L,
  overrides_1 = list(),
  overrides_2 = list()
) {
  make_params <- function(population, overrides) {
    values <- list(
      human_population = population,
      human_mobility_enabled = TRUE,
      human_move_probs = diag(2),
      native_mosquito_backend = TRUE,
      individual_mosquitoes = FALSE,
      total_M = 20,
      delay_gam = delay_gam,
      progress_bar = FALSE
    )
    values[names(overrides)] <- overrides
    get_parameters(values)
  }

  parameters <- list(make_params(population_1, overrides_1), make_params(population_2, overrides_2))
  for (i in seq_along(parameters)) {
    parameters[[i]]$human_mobility_node_index <- as.integer(i)
  }
  parameters
}

human_mobility_stage5_set_node <- function(
  variables,
  infectivity,
  zeta = rep(1, length(infectivity)),
  birth = rep(-1000, length(infectivity)),
  current_node = NULL,
  states = rep("S", length(infectivity)),
  tbv_vaccinated = rep(-1, length(infectivity))
) {
  variables$infectivity <- individual::DoubleVariable$new(infectivity)
  variables$zeta <- individual::DoubleVariable$new(zeta)
  variables$birth <- individual::DoubleVariable$new(birth)
  variables$state <- individual::CategoricalVariable$new(
    categories = c("S", "A", "U", "D", "Tr"),
    initial_values = states
  )
  variables$tbv_vaccinated <- individual::DoubleVariable$new(tbv_vaccinated)
  if (!is.null(current_node)) {
    variables$current_node <- individual::IntegerVariable$new(as.integer(current_node))
  }
  variables
}

human_mobility_stage5_context <- function(parameters, variables) {
  create_human_infectivity_lag_context(
    parameters = parameters,
    variables = variables
  )
}

human_mobility_stage5_run <- function(parameters, timesteps = 3L, return_state = FALSE, initial_state = NULL) {
  run_metapop_simulation(
    timesteps = timesteps,
    parameters = parameters,
    mixing_tt = 1,
    export_mixing = list(diag(2)),
    import_mixing = list(diag(2)),
    p_captured_tt = 1,
    p_captured = list(matrix(0, nrow = 2, ncol = 2)),
    p_success = 0,
    return_state = return_state,
    initial_state = initial_state,
    restore_random_state = TRUE
  )
}

test_that("human infectivity lag buffers use defaults, interpolate, and clear to zero", {
  buffer <- HumanInfectivityLagBuffer$new(
    max_lag = 3.5,
    default_infectivity = c(0.1, 0.2)
  )

  expect_equal(buffer$get(0), c(0.1, 0.2))

  buffer$save(1, c(0.2, 0.4))
  buffer$save(3, c(0.6, 1.0))

  expect_equal(buffer$get(1), c(0.2, 0.4))
  expect_equal(buffer$get(2), c(0.4, 0.7))
  expect_warning(
    expect_equal(buffer$get(4), c(0.1, 0.2)),
    "after the latest saved timestep"
  )

  buffer$clear(2L)
  expect_equal(buffer$get(0), c(0.1, 0))
  expect_equal(buffer$get(1), c(0.2, 0))
  expect_equal(buffer$get(3), c(0.6, 0))
})

test_that("disabled mobility creates no human infectivity lag context", {
  parameters <- human_mobility_stage5_params()
  parameters <- lapply(parameters, function(p) {
    p$human_mobility_enabled <- FALSE
    p
  })

  expect_null(create_human_infectivity_lag_context(
    parameters = parameters,
    variables = lapply(parameters, create_variables)
  ))
})

test_that("diag human_move_probs reproduces scalar FOIM reservoir input with equal weights", {
  parameters <- human_mobility_stage5_params(population_1 = 2L, population_2 = 3L)
  variables <- lapply(parameters, create_variables)
  variables[[1L]] <- human_mobility_stage5_set_node(variables[[1L]], c(0.2, 0.4))
  variables[[2L]] <- human_mobility_stage5_set_node(variables[[2L]], c(0.1, 0.3, 0.5))
  context <- human_mobility_stage5_context(parameters, variables)

  human_infectivity_lag_record_node(context, 1L, timestep = 1L)
  human_infectivity_lag_record_node(context, 2L, timestep = 1L)

  reservoir <- human_infectivity_lag_get_reservoir(context, timestep = 2L)

  expect_equal(reservoir[[1L]], mean(c(0.2, 0.4)))
  expect_equal(reservoir[[2L]], mean(c(0.1, 0.3, 0.5)))
})

test_that("travelling infectious humans contribute to their current destination reservoir", {
  parameters <- human_mobility_stage5_params(population_1 = 2L, population_2 = 2L)
  variables <- lapply(parameters, create_variables)
  variables[[1L]] <- human_mobility_stage5_set_node(
    variables[[1L]],
    infectivity = c(1, 0),
    current_node = c(2L, 1L)
  )
  variables[[2L]] <- human_mobility_stage5_set_node(
    variables[[2L]],
    infectivity = c(0, 0),
    current_node = c(2L, 2L)
  )
  context <- human_mobility_stage5_context(parameters, variables)

  human_infectivity_lag_record_node(context, 1L, timestep = 1L)
  human_infectivity_lag_record_node(context, 2L, timestep = 1L)
  reservoir <- human_infectivity_lag_get_reservoir(context, timestep = 2L)

  expect_equal(reservoir[[1L]], 0)
  expect_equal(reservoir[[2L]], 1 / 3)
})

test_that("location is current at lookup timestep, not lagged with infectivity", {
  parameters <- human_mobility_stage5_params(population_1 = 1L, population_2 = 1L)
  variables <- lapply(parameters, create_variables)
  variables[[1L]] <- human_mobility_stage5_set_node(
    variables[[1L]],
    infectivity = 1,
    current_node = 1L
  )
  variables[[2L]] <- human_mobility_stage5_set_node(
    variables[[2L]],
    infectivity = 0,
    current_node = 2L
  )
  context <- human_mobility_stage5_context(parameters, variables)

  human_infectivity_lag_record_node(context, 1L, timestep = 1L)
  human_infectivity_lag_record_node(context, 2L, timestep = 1L)

  variables[[1L]]$current_node$queue_update(2L, 1L)
  variables[[1L]]$current_node$.update()
  reservoir <- human_infectivity_lag_get_reservoir(context, timestep = 2L)

  expect_equal(reservoir[[1L]], 0)
  expect_equal(reservoir[[2L]], 0.5)
})

test_that("new infectivity contributes only after delay_gam has elapsed", {
  parameters <- human_mobility_stage5_params(population_1 = 1L, population_2 = 1L)
  variables <- lapply(parameters, create_variables)
  variables[[1L]] <- human_mobility_stage5_set_node(variables[[1L]], infectivity = 0)
  variables[[2L]] <- human_mobility_stage5_set_node(variables[[2L]], infectivity = 0)
  context <- human_mobility_stage5_context(parameters, variables)

  human_infectivity_lag_record_node(context, 1L, timestep = 1L)
  human_infectivity_lag_record_node(context, 2L, timestep = 1L)
  variables[[1L]]$infectivity$queue_update(1, 1L)
  variables[[1L]]$infectivity$.update()
  human_infectivity_lag_record_node(context, 1L, timestep = 2L)
  human_infectivity_lag_record_node(context, 2L, timestep = 2L)

  expect_equal(human_infectivity_lag_get_reservoir(context, timestep = 2L)[[1L]], 0)
  expect_equal(human_infectivity_lag_get_reservoir(context, timestep = 3L)[[1L]], 1)
})

test_that("delay_gam zero uses same-timestep per-human infectivity", {
  parameters <- human_mobility_stage5_params(
    delay_gam = 0,
    population_1 = 1L,
    population_2 = 1L
  )
  variables <- lapply(parameters, create_variables)
  variables[[1L]] <- human_mobility_stage5_set_node(variables[[1L]], infectivity = 0.7)
  variables[[2L]] <- human_mobility_stage5_set_node(variables[[2L]], infectivity = 0)
  context <- human_mobility_stage5_context(parameters, variables)

  human_infectivity_lag_record_node(context, 1L, timestep = 5L)
  human_infectivity_lag_record_node(context, 2L, timestep = 5L)

  expect_equal(human_infectivity_lag_get_reservoir(context, timestep = 5L)[[1L]], 0.7)
})

test_that("empty destination nodes have zero reservoir input", {
  parameters <- human_mobility_stage5_params(population_1 = 1L, population_2 = 1L)
  variables <- lapply(parameters, create_variables)
  variables[[1L]] <- human_mobility_stage5_set_node(
    variables[[1L]],
    infectivity = 1,
    current_node = 1L
  )
  variables[[2L]] <- human_mobility_stage5_set_node(
    variables[[2L]],
    infectivity = 0,
    current_node = 1L
  )
  context <- human_mobility_stage5_context(parameters, variables)

  human_infectivity_lag_record_node(context, 1L, timestep = 1L)
  human_infectivity_lag_record_node(context, 2L, timestep = 1L)

  expect_equal(human_infectivity_lag_get_reservoir(context, timestep = 2L)[[2L]], 0)
})

test_that("reset_target clears replaced human infectivity lag history to zero", {
  parameters <- human_mobility_stage5_params(population_1 = 3L)
  variables <- lapply(parameters, create_variables)
  variables[[1L]] <- human_mobility_stage5_set_node(
    variables[[1L]],
    infectivity = c(0.1, 0.2, 0.3)
  )
  context <- human_mobility_stage5_context(parameters, variables)
  context$buffers[[1L]]$save(1, c(0.4, 0.5, 0.6))
  context$buffers[[1L]]$save(2, c(0.7, 0.8, 0.9))

  target <- individual::Bitset$new(3)$insert(2L)
  reset_target(
    variables[[1L]],
    create_events(parameters[[1L]]),
    target,
    "S",
    parameters[[1L]],
    timestep = 3L,
    human_infectivity_lag_context = context
  )

  expect_equal(context$buffers[[1L]]$get(0), c(0.1, 0, 0.3))
  expect_equal(context$buffers[[1L]]$get(1), c(0.4, 0, 0.6))
  expect_equal(context$buffers[[1L]]$get(2), c(0.7, 0, 0.9))
})

test_that("TBV-adjusted infectivity is recorded in the per-human lag buffer", {
  parameters <- human_mobility_stage5_params(population_1 = 5L, population_2 = 1L)
  parameters[[1L]]$tbv <- TRUE
  variables <- lapply(parameters, create_variables)
  infectivity <- c(0, 0.1, 0.15, 0.5, 0.3)
  variables[[1L]] <- human_mobility_stage5_set_node(
    variables[[1L]],
    infectivity = infectivity,
    states = c("S", "U", "A", "D", "Tr"),
    tbv_vaccinated = c(-1, -1, 50, 50, 50)
  )
  context <- human_mobility_stage5_context(parameters, variables)

  human_infectivity_lag_record_node(context, 1L, timestep = 55L)

  expect_equal(
    context$buffers[[1L]]$get(55L),
    account_for_tbv(55L, infectivity, variables[[1L]], parameters[[1L]])
  )
})

test_that("human infectivity lag state is resumable and older states are tolerated", {
  parameters <- human_mobility_stage5_params(
    population_1 = 2L,
    population_2 = 2L,
    overrides_1 = list(average_age = 1e12),
    overrides_2 = list(average_age = 1e12)
  )

  set.seed(123)
  checkpoint <- human_mobility_stage5_run(parameters, timesteps = 2L, return_state = TRUE)
  expect_true("human_infectivity_lag" %in% names(checkpoint$state$malariasimulationGD))

  resumed <- human_mobility_stage5_run(
    parameters,
    timesteps = 4L,
    return_state = TRUE,
    initial_state = checkpoint$state
  )
  for (node in seq_along(checkpoint$state$malariasimulationGD$human_infectivity_lag)) {
    checkpoint_buffer <- checkpoint$state$malariasimulationGD$human_infectivity_lag[[node]]
    resumed_buffer <- resumed$state$malariasimulationGD$human_infectivity_lag[[node]]
    restored_rows <- match(checkpoint_buffer$timesteps, resumed_buffer$timesteps)
    expect_false(anyNA(restored_rows))
    expect_equal(
      resumed_buffer$infectivity[restored_rows, , drop = FALSE],
      checkpoint_buffer$infectivity,
      ignore_attr = TRUE
    )
  }

  old_state <- checkpoint$state
  old_state$malariasimulationGD$human_infectivity_lag <- NULL
  expect_silent(human_mobility_stage5_run(
    parameters,
    timesteps = 4L,
    return_state = TRUE,
    initial_state = old_state
  ))
})
