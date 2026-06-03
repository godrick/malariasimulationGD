human_mobility_stage2_matrix <- function(duration_test = TRUE) {
  if (isTRUE(duration_test)) {
    return(matrix(c(0, 1, 0, 1), nrow = 2, byrow = TRUE))
  }
  matrix(c(0.75, 0.25, 0, 1), nrow = 2, byrow = TRUE)
}

human_mobility_stage2_params <- function(
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
  list(make_params(population_1, overrides_1), make_params(population_2, overrides_2))
}

run_human_mobility_stage2 <- function(parameters, timesteps = 3L, return_state = FALSE) {
  run_metapop_simulation(
    timesteps = timesteps,
    parameters = parameters,
    mixing_tt = 1,
    export_mixing = list(diag(2)),
    import_mixing = list(diag(2)),
    p_captured_tt = 1,
    p_captured = list(matrix(0, nrow = 2, ncol = 2)),
    p_success = 0,
    return_state = return_state
  )
}

test_that("mobility variables are added only when enabled and initialize at home", {
  disabled <- get_parameters(list(human_population = 4))
  disabled_vars <- create_variables(disabled)
  expect_true(is.null(disabled_vars$home_node))
  expect_true(is.null(disabled_vars$current_node))

  enabled <- get_parameters(list(
    human_population = 4,
    human_mobility_enabled = TRUE
  ))
  enabled$human_mobility_node_index <- 2L
  enabled_vars <- create_variables(enabled)
  expect_equal(enabled_vars$home_node$get_values(), rep(2L, 4))
  expect_equal(enabled_vars$current_node$get_values(), rep(2L, 4))
  expect_equal(enabled_vars$travel_destination$get_values(), rep(2L, 4))
  expect_equal(enabled_vars$travel_remaining_nights$get_values(), rep(0L, 4))
  expect_equal(enabled_vars$is_travelling$get_values(), rep(0L, 4))
})

test_that("mobility disabled runs without mobility outputs", {
  output <- run_human_mobility_stage2(human_mobility_stage2_params(), timesteps = 1)

  expect_false("humans_present" %in% names(output[[1]]))
  expect_false("visitors_present" %in% names(output[[1]]))
  expect_false("residents_away" %in% names(output[[1]]))
  expect_false("trips_started" %in% names(output[[1]]))
})

test_that("stochastic resample initialization starts mobility state at home", {
  parameters <- get_parameters(list(
    human_population = 2L,
    human_initialization = "stochastic_resample",
    human_initialization_burnin_timesteps = 1L,
    human_mobility_enabled = TRUE,
    progress_bar = FALSE
  ))
  parameters$init_EIR <- 1
  parameters$init_foim <- 0
  parameters$human_mobility_node_index <- 2L

  sampled <- list(
    state = c("S", "A"),
    birth = c(-10L, -25L),
    last_boosted_ica = c(-1, -1),
    icm = c(0.1, 0.2),
    ica = c(1, 2),
    zeta = c(1, 1.5),
    zeta_group = c("1", "2"),
    infectivity = c(0, 0.05),
    progression_rates = c(0, 1 / 200),
    drug = c(0L, 0L),
    drug_time = c(-1L, -1L),
    last_pev_timestep = c(-1L, -1L),
    last_eff_pev_timestep = c(-1L, -1L),
    pev_profile = c(-1L, -1L),
    tbv_vaccinated = c(-1, -1),
    net_time = c(-1L, -1L),
    spray_time = c(-1L, -1L),
    last_boosted_ib = c(-1, -1),
    last_boosted_iva = c(-1, -1),
    last_boosted_id = c(-1, -1),
    ivm = c(0.01, 0.02),
    ib = c(1, 2),
    iva = c(1, 2),
    id = c(1, 2)
  )

  create_variables_stub <- create_variables
  mockery::stub(
    create_variables_stub,
    "stationary_human_initializer_resample",
    function(library, size, parameters) sampled
  )
  mockery::stub(
    create_variables_stub,
    "get_stationary_human_initialization_library",
    function(parameters) list(timesteps = 20L, variables = sampled)
  )

  variables <- create_variables_stub(parameters)

  expect_equal(variables$home_node$get_values(), c(2L, 2L))
  expect_equal(variables$current_node$get_values(), c(2L, 2L))
  expect_equal(variables$travel_remaining_nights$get_values(), c(0L, 0L))
})

test_that("reset_target resets mobility state to permanent home", {
  parameters <- get_parameters(list(
    human_population = 4,
    human_mobility_enabled = TRUE
  ))
  parameters$human_mobility_node_index <- 1L
  variables <- create_variables(parameters)
  events <- create_events(parameters)
  target <- individual::Bitset$new(4)$insert(c(2, 4))

  variables$current_node$queue_update(2L, target)
  variables$travel_destination$queue_update(2L, target)
  variables$travel_remaining_nights$queue_update(3L, target)
  variables$is_travelling$queue_update(1L, target)
  for (variable in variables[c("current_node", "travel_destination", "travel_remaining_nights", "is_travelling")]) {
    variable$.update()
  }

  reset_target(variables, events, target, "S", parameters, timestep = 5L)
  for (variable in variables[c("current_node", "travel_destination", "travel_remaining_nights", "is_travelling")]) {
    variable$.update()
  }

  expect_equal(variables$home_node$get_values(), rep(1L, 4))
  expect_equal(variables$current_node$get_values(c(2, 4)), c(1L, 1L))
  expect_equal(variables$travel_destination$get_values(c(2, 4)), c(1L, 1L))
  expect_equal(variables$travel_remaining_nights$get_values(c(2, 4)), c(0L, 0L))
  expect_equal(variables$is_travelling$get_values(c(2, 4)), c(0L, 0L))
})

test_that("mobility process is inserted after mosquito releases and before biting", {
  parameters <- human_mobility_stage2_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage2_matrix()
  ))
  parameters <- lapply(seq_along(parameters), function(i) {
    parameters[[i]]$human_mobility_enabled <- TRUE
    parameters[[i]]$human_move_probs <- human_mobility_stage2_matrix()
    parameters[[i]]$human_mobility_node_index <- as.integer(i)
    parameters[[i]]
  })
  variables <- lapply(parameters, create_variables)
  renderers <- lapply(seq_along(parameters), function(.) individual::Render$new(1))
  context <- create_human_mobility_context(parameters, variables, timesteps = 1)
  processes <- create_processes(
    renderer = renderers[[1]],
    variables = variables[[1]],
    events = create_events(parameters[[1]]),
    parameters = parameters[[1]],
    models = list(),
    solvers = list(),
    lagged_eir = list(),
    lagged_infectivity = list(),
    lagged_transmission_eir = list(),
    human_mobility_context = context
  )

  process_names <- names(processes)
  expect_lt(match("mosquito_release_process", process_names), match("human_mobility_process", process_names))
  expect_lt(match("human_mobility_process", process_names), match("biting_process", process_names))

  disabled <- get_parameters(list(
    human_population = 4,
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE
  ))
  disabled_processes <- create_processes(
    renderer = individual::Render$new(1),
    variables = create_variables(disabled),
    events = create_events(disabled),
    parameters = disabled,
    models = list(),
    solvers = list(),
    lagged_eir = list(),
    lagged_infectivity = list(),
    lagged_transmission_eir = list()
  )
  expect_false("human_mobility_process" %in% names(disabled_processes))
})

test_that("fixed L=1 is away for current timestep and home next timestep", {
  parameters <- human_mobility_stage2_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage2_matrix(),
    human_trip_duration_type = "fixed",
    human_trip_duration_mean = 1
  ))

  output <- run_human_mobility_stage2(parameters, timesteps = 3)

  expect_equal(output[[1]]$residents_away, c(2, 0, 2))
  expect_equal(output[[1]]$trips_started, c(2, 0, 2))
  expect_equal(output[[1]]$humans_present, c(0, 2, 0))
  expect_equal(output[[2]]$visitors_present, c(2, 0, 2))
})

test_that("fixed L=2 is away for current and next timestep, home at t+2", {
  parameters <- human_mobility_stage2_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage2_matrix(),
    human_trip_duration_type = "fixed",
    human_trip_duration_mean = 2
  ))

  output <- run_human_mobility_stage2(parameters, timesteps = 3)

  expect_equal(output[[1]]$residents_away, c(2, 2, 0))
  expect_equal(output[[1]]$trips_started, c(2, 0, 0))
  expect_equal(output[[1]]$humans_present, c(0, 0, 2))
  expect_equal(output[[2]]$visitors_present, c(2, 2, 0))
})

test_that("traveller with one remaining night stays away this timestep", {
  parameters <- human_mobility_stage2_params(
    list(
      human_mobility_enabled = TRUE,
      human_move_probs = human_mobility_stage2_matrix()
    ),
    population_1 = 2L,
    population_2 = 1L
  )
  parameters <- lapply(seq_along(parameters), function(i) {
    parameters[[i]]$human_mobility_enabled <- TRUE
    parameters[[i]]$human_move_probs <- human_mobility_stage2_matrix()
    parameters[[i]]$human_mobility_node_index <- as.integer(i)
    parameters[[i]]
  })
  variables <- lapply(parameters, create_variables)
  variables[[1]]$current_node$queue_update(2L, 1L)
  variables[[1]]$travel_destination$queue_update(2L, 1L)
  variables[[1]]$travel_remaining_nights$queue_update(1L, 1L)
  variables[[1]]$is_travelling$queue_update(1L, 1L)
  for (variable in variables[[1]][c("current_node", "travel_destination", "travel_remaining_nights", "is_travelling")]) {
    variable$.update()
  }

  renderer <- individual::Render$new(1)
  context <- create_human_mobility_context(parameters, variables, timesteps = 1)
  process <- create_human_mobility_process(context, node_index = 1L, renderer = renderer)
  process(1L)

  expect_equal(variables[[1]]$current_node$get_values(), c(2L, 2L))
  expect_equal(renderer$to_dataframe()$residents_away, 2)

  for (variable in variables[[1]][c("current_node", "travel_destination", "travel_remaining_nights", "is_travelling")]) {
    variable$.update()
  }
  expect_equal(variables[[1]]$current_node$get_values(), c(1L, 1L))
  expect_equal(variables[[1]]$travel_remaining_nights$get_values(), c(-1L, -1L))
})

test_that("geometric mean one always produces one-night trips", {
  parameters <- human_mobility_stage2_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage2_matrix(),
    human_trip_duration_type = "geometric",
    human_trip_duration_mean = 1
  ))

  output <- run_human_mobility_stage2(parameters, timesteps = 2)

  expect_equal(output[[1]]$residents_away, c(2, 0))
  expect_equal(output[[1]]$trips_started, c(2, 0))
})

test_that("shared matrix supplied once is consumed by all nodes", {
  parameters <- human_mobility_stage2_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage2_matrix(),
    human_trip_duration_type = "fixed",
    human_trip_duration_mean = 1
  ))

  output <- run_human_mobility_stage2(parameters, timesteps = 1)

  expect_true("humans_present" %in% names(output[[1]]))
  expect_true("humans_present" %in% names(output[[2]]))
  expect_equal(output[[2]]$visitors_present, 2)
})

test_that("identical matrices supplied on multiple nodes are accepted", {
  parameters <- human_mobility_stage2_params(
    list(
      human_mobility_enabled = TRUE,
      human_move_probs = human_mobility_stage2_matrix(),
      human_trip_duration_type = "fixed",
      human_trip_duration_mean = 1
    ),
    list(human_move_probs = human_mobility_stage2_matrix())
  )

  output <- run_human_mobility_stage2(parameters, timesteps = 1)

  expect_equal(output[[1]]$residents_away, 2)
  expect_equal(output[[2]]$visitors_present, 2)
})

test_that("mobility diagnostics are optional output-list attributes", {
  parameters <- human_mobility_stage2_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage2_matrix(),
    human_mobility_store_diagnostics = TRUE
  ))

  output <- run_human_mobility_stage2(parameters, timesteps = 2)

  expect_false("OD_started_trips" %in% names(output[[1]]))
  expect_equal(dim(attr(output, "OD_started_trips")), c(2L, 2L, 2L))
  expect_equal(dim(attr(output, "OD_active_overnight_stays")), c(2L, 2L, 2L))
  expect_equal(dim(attr(output, "mean_remaining_trip_duration")), c(2L, 2L))

  parameters[[1]]$human_mobility_store_diagnostics <- FALSE
  output_without_diagnostics <- run_human_mobility_stage2(parameters, timesteps = 1)
  expect_null(attr(output_without_diagnostics, "OD_started_trips"))
})

test_that("mobility outputs are resumable", {
  parameters <- human_mobility_stage2_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage2_matrix(),
    human_trip_duration_type = "fixed",
    human_trip_duration_mean = 1,
    human_mobility_store_diagnostics = TRUE,
    average_age = 1e12
  ))

  set.seed(123)
  full <- run_human_mobility_stage2(parameters, timesteps = 4, return_state = TRUE)
  set.seed(123)
  checkpoint <- run_human_mobility_stage2(parameters, timesteps = 2, return_state = TRUE)
  set.seed(456)
  resumed <- run_metapop_simulation(
    timesteps = 4,
    parameters = parameters,
    mixing_tt = 1,
    export_mixing = list(diag(2)),
    import_mixing = list(diag(2)),
    p_captured_tt = 1,
    p_captured = list(matrix(0, nrow = 2, ncol = 2)),
    p_success = 0,
    initial_state = checkpoint$state,
    restore_random_state = TRUE,
    return_state = TRUE
  )

  mobility_cols <- c("humans_present", "visitors_present", "residents_away", "trips_started")
  expect_equal(
    resumed$data[[1]][mobility_cols],
    full$data[[1]][3:4, mobility_cols, drop = FALSE],
    ignore_attr = TRUE
  )
  expect_equal(
    resumed$data[[2]][mobility_cols],
    full$data[[2]][3:4, mobility_cols, drop = FALSE],
    ignore_attr = TRUE
  )
  expect_equal(
    attr(resumed$data, "OD_started_trips"),
    attr(full$data, "OD_started_trips")[3:4, , , drop = FALSE],
    ignore_attr = TRUE
  )
})
