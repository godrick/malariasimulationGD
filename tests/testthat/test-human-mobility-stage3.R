human_mobility_stage3_matrix <- function(stay_home = TRUE) {
  if (isTRUE(stay_home)) {
    return(diag(2))
  }
  matrix(c(0, 1, 0, 1), nrow = 2, byrow = TRUE)
}

human_mobility_stage3_params <- function(
  overrides_1 = list(),
  overrides_2 = list(),
  population_1 = 2L,
  population_2 = 3L
) {
  base <- list(
    total_M = 20,
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    progress_bar = FALSE
  )
  make_params <- function(population, overrides) {
    values <- c(base, list(human_population = population))
    values[names(overrides)] <- overrides
    get_parameters(values)
  }
  parameters <- list(make_params(population_1, overrides_1), make_params(population_2, overrides_2))
  for (i in seq_along(parameters)) {
    parameters[[i]]$human_mobility_enabled <- TRUE
    if (is.null(parameters[[i]]$human_move_probs)) {
      parameters[[i]]$human_move_probs <- human_mobility_stage3_matrix()
    }
    parameters[[i]]$human_mobility_node_index <- as.integer(i)
  }
  parameters
}

human_mobility_stage3_lagged <- function(defaults) {
  lapply(defaults, function(default) list(LaggedValue$new(3, default)))
}

run_human_mobility_stage3 <- function(parameters, timesteps = 3L, return_state = FALSE, initial_state = NULL) {
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

test_that("human exposure lag buffers use per-human defaults and interpolate", {
  buffer <- HumanExposureLagBuffer$new(
    max_lag = 3.5,
    default_exposure = c(10, 20),
    default_weighted_exposure = c(30, 40)
  )

  expect_equal(buffer$get(0), c(10, 20))
  expect_equal(buffer$get(0, weighted = TRUE), c(30, 40))

  buffer$save(1, c(2, 4), c(6, 8))
  buffer$save(3, c(6, 10), c(10, 16))

  expect_equal(buffer$get(1), c(2, 4))
  expect_equal(buffer$get(2), c(4, 7))
  expect_equal(buffer$get(2, weighted = TRUE), c(8, 12))
  expect_warning(
    expect_equal(buffer$get(4), c(10, 20)),
    "after the latest saved timestep"
  )
})

test_that("diag human_move_probs records allocated home-node unlagged exposure for residents", {
  parameters <- human_mobility_stage3_params()
  variables <- lapply(parameters, create_variables)
  variables[[1L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[1L]]$human_population))
  variables[[2L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[2L]]$human_population))
  variables[[1L]]$zeta <- individual::DoubleVariable$new(rep(1, parameters[[1L]]$human_population))
  variables[[2L]]$zeta <- individual::DoubleVariable$new(rep(1, parameters[[2L]]$human_population))
  lagged_eir <- human_mobility_stage3_lagged(c(11, 22))

  context <- create_human_exposure_lag_context(
    parameters = parameters,
    variables = variables,
    solvers = list(NULL, NULL),
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_eir
  )

  human_exposure_lag_record_node(context, node_index = 1L, timestep = 1L, exposure = 5)
  expect_equal(context$buffers[[1]]$get(1), rep(11 / 2, 2))

  human_exposure_lag_record_node(context, node_index = 2L, timestep = 1L, exposure = 7)

  expect_equal(context$buffers[[1]]$get(1), rep(5 / 2, 2))
  expect_equal(context$buffers[[2]]$get(1), rep(7 / 3, 3))
})

test_that("recording allocates destination exposure across explicit travellers and residents", {
  parameters <- human_mobility_stage3_params()
  variables <- lapply(parameters, create_variables)
  variables[[1L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[1L]]$human_population))
  variables[[2L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[2L]]$human_population))
  variables[[1L]]$zeta <- individual::DoubleVariable$new(rep(1, parameters[[1L]]$human_population))
  variables[[2L]]$zeta <- individual::DoubleVariable$new(rep(1, parameters[[2L]]$human_population))
  lagged_eir <- human_mobility_stage3_lagged(c(11, 22))

  variables[[1]]$current_node$queue_update(2L, 1L)
  variables[[1]]$current_node$.update()

  context <- create_human_exposure_lag_context(
    parameters = parameters,
    variables = variables,
    solvers = list(NULL, NULL),
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_eir
  )

  human_exposure_lag_record_node(context, node_index = 1L, timestep = 1L, exposure = 5)
  human_exposure_lag_record_node(context, node_index = 2L, timestep = 1L, exposure = 7)

  expect_equal(context$buffers[[1]]$get(1), c(7 / 4, 5))
  expect_equal(context$buffers[[2]]$get(1), rep(7 / 4, 3))
})

test_that("reset_target clears replaced slots to home-node exposure defaults", {
  parameters <- human_mobility_stage3_params(population_1 = 3L)
  variables <- lapply(parameters, create_variables)
  variables[[1L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[1L]]$human_population))
  variables[[1L]]$zeta <- individual::DoubleVariable$new(rep(1, parameters[[1L]]$human_population))
  variables[[2L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[2L]]$human_population))
  variables[[2L]]$zeta <- individual::DoubleVariable$new(rep(1, parameters[[2L]]$human_population))
  lagged_eir <- human_mobility_stage3_lagged(c(11, 22))
  context <- create_human_exposure_lag_context(
    parameters = parameters,
    variables = variables,
    solvers = list(NULL, NULL),
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_eir
  )
  context$buffers[[1]]$save(1, c(3, 4, 5))
  context$buffers[[1]]$save(2, c(6, 7, 8))

  target <- individual::Bitset$new(3)$insert(2L)
  reset_target(
    variables[[1]],
    create_events(parameters[[1]]),
    target,
    "S",
    parameters[[1]],
    timestep = 3L,
    human_exposure_lag_context = context
  )

  expect_equal(context$buffers[[1]]$get(1), c(3, 11 / 3, 5))
  expect_equal(context$buffers[[1]]$get(2), c(6, 11 / 3, 8))
})

test_that("weighted exposure is optional and scalar per human when active", {
  inactive_parameters <- human_mobility_stage3_params()
  inactive_variables <- lapply(inactive_parameters, create_variables)
  lagged_eir <- human_mobility_stage3_lagged(c(11, 22))
  inactive_context <- create_human_exposure_lag_context(
    parameters = inactive_parameters,
    variables = inactive_variables,
    solvers = list(NULL, NULL),
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_eir
  )
  expect_null(inactive_context$buffers[[1]]$get(0, weighted = TRUE))

  active_parameters <- human_mobility_stage3_params()
  active_parameters <- lapply(active_parameters, function(p) {
    p$vector_infectivity_g_by_species <- list(gamb = 2)
    p
  })
  active_variables <- lapply(active_parameters, create_variables)
  active_variables[[1L]]$birth <- individual::DoubleVariable$new(rep(-1e12, active_parameters[[1L]]$human_population))
  active_variables[[1L]]$zeta <- individual::DoubleVariable$new(rep(1, active_parameters[[1L]]$human_population))
  active_variables[[2L]]$birth <- individual::DoubleVariable$new(rep(-1e12, active_parameters[[2L]]$human_population))
  active_variables[[2L]]$zeta <- individual::DoubleVariable$new(rep(1, active_parameters[[2L]]$human_population))
  lagged_transmission_eir <- human_mobility_stage3_lagged(c(101, 202))
  active_context <- create_human_exposure_lag_context(
    parameters = active_parameters,
    variables = active_variables,
    solvers = list(NULL, NULL),
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_transmission_eir
  )

  expect_equal(active_context$buffers[[1]]$get(0, weighted = TRUE), rep(101 / 2, 2))
  human_exposure_lag_record_node(active_context, node_index = 1L, timestep = 1L, exposure = 5, weighted_exposure = 15)
  human_exposure_lag_record_node(active_context, node_index = 2L, timestep = 1L, exposure = 7, weighted_exposure = 25)

  expect_equal(active_context$buffers[[1]]$get(1), rep(5 / 2, 2))
  expect_equal(active_context$buffers[[1]]$get(1, weighted = TRUE), rep(15 / 2, 2))
})

test_that("human exposure allocation uses biting weights and preserves destination totals", {
  parameters <- human_mobility_stage3_params(population_1 = 2L, population_2 = 1L)
  variables <- lapply(parameters, create_variables)
  variables[[1L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[1L]]$human_population))
  variables[[2L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[2L]]$human_population))
  variables[[1L]]$zeta <- individual::DoubleVariable$new(c(1, 3))
  variables[[2L]]$zeta <- individual::DoubleVariable$new(1)

  context <- create_human_exposure_lag_context(
    parameters = parameters,
    variables = variables,
    solvers = list(NULL, NULL),
    lagged_eir = human_mobility_stage3_lagged(c(0, 0)),
    lagged_transmission_eir = human_mobility_stage3_lagged(c(0, 0))
  )

  human_exposure_lag_record_node(context, node_index = 1L, timestep = 1L, exposure = 8)
  human_exposure_lag_record_node(context, node_index = 2L, timestep = 1L, exposure = 9)

  node_1_exposure <- context$buffers[[1L]]$get(1L)
  node_2_exposure <- context$buffers[[2L]]$get(1L)
  expect_equal(node_1_exposure, c(2, 6), tolerance = 1e-12)
  expect_equal(node_2_exposure, 9, tolerance = 1e-12)
  expect_equal(sum(node_1_exposure), 8, tolerance = 1e-12)
  expect_equal(sum(node_2_exposure), 9, tolerance = 1e-12)
  expect_lt(max(node_1_exposure), 8)
})

test_that("weighted human exposure uses the same biting allocation shares", {
  parameters <- human_mobility_stage3_params(population_1 = 2L, population_2 = 1L)
  parameters <- lapply(parameters, function(p) {
    p$vector_infectivity_g_by_species <- list(gamb = 2)
    p
  })
  variables <- lapply(parameters, create_variables)
  variables[[1L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[1L]]$human_population))
  variables[[2L]]$birth <- individual::DoubleVariable$new(rep(-1e12, parameters[[2L]]$human_population))
  variables[[1L]]$zeta <- individual::DoubleVariable$new(c(1, 3))
  variables[[2L]]$zeta <- individual::DoubleVariable$new(1)

  context <- create_human_exposure_lag_context(
    parameters = parameters,
    variables = variables,
    solvers = list(NULL, NULL),
    lagged_eir = human_mobility_stage3_lagged(c(0, 0)),
    lagged_transmission_eir = human_mobility_stage3_lagged(c(0, 0))
  )

  human_exposure_lag_record_node(context, node_index = 1L, timestep = 1L, exposure = 8, weighted_exposure = 20)
  human_exposure_lag_record_node(context, node_index = 2L, timestep = 1L, exposure = 9, weighted_exposure = 45)

  expect_equal(context$buffers[[1L]]$get(1L), c(2, 6), tolerance = 1e-12)
  expect_equal(context$buffers[[1L]]$get(1L, weighted = TRUE), c(5, 15), tolerance = 1e-12)
  expect_equal(context$buffers[[2L]]$get(1L), 9, tolerance = 1e-12)
  expect_equal(context$buffers[[2L]]$get(1L, weighted = TRUE), 45, tolerance = 1e-12)
})

test_that("disabled mobility creates no human exposure lag context", {
  parameters <- lapply(
    human_mobility_stage3_params(),
    function(p) {
      p$human_mobility_enabled <- FALSE
      p
    }
  )

  expect_null(create_human_exposure_lag_context(
    parameters = parameters,
    variables = lapply(parameters, create_variables),
    solvers = list(NULL, NULL),
    lagged_eir = human_mobility_stage3_lagged(c(11, 22)),
    lagged_transmission_eir = human_mobility_stage3_lagged(c(11, 22))
  ))
})

test_that("human exposure lag state is resumable and older states are tolerated", {
  parameters <- human_mobility_stage3_params(list(
    human_move_probs = human_mobility_stage3_matrix(stay_home = FALSE),
    human_trip_duration_type = "fixed",
    human_trip_duration_mean = 1,
    average_age = 1e12
  ), list(
    human_move_probs = human_mobility_stage3_matrix(stay_home = FALSE),
    average_age = 1e12
  ))

  set.seed(123)
  checkpoint <- run_human_mobility_stage3(parameters, timesteps = 2, return_state = TRUE)
  resumed <- run_human_mobility_stage3(
    parameters,
    timesteps = 4,
    return_state = TRUE,
    initial_state = checkpoint$state
  )

  for (node in seq_along(checkpoint$state$malariasimulationGD$human_exposure_lag)) {
    checkpoint_buffer <- checkpoint$state$malariasimulationGD$human_exposure_lag[[node]]
    resumed_buffer <- resumed$state$malariasimulationGD$human_exposure_lag[[node]]
    restored_rows <- match(checkpoint_buffer$timesteps, resumed_buffer$timesteps)
    expect_false(anyNA(restored_rows))
    expect_equal(
      resumed_buffer$exposure[restored_rows, , drop = FALSE],
      checkpoint_buffer$exposure,
      ignore_attr = TRUE
    )
  }

  old_state <- checkpoint$state
  old_state$malariasimulationGD$human_exposure_lag <- NULL
  expect_silent(run_human_mobility_stage3(
    parameters,
    timesteps = 4,
    return_state = TRUE,
    initial_state = old_state
  ))
})
