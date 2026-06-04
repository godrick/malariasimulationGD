human_mobility_stage7_matrix <- function() {
  matrix(c(0.8, 0.2, 0.1, 0.9), nrow = 2, byrow = TRUE)
}

human_mobility_stage7_params <- function(
  overrides_1 = list(),
  overrides_2 = list(),
  population_1 = 3L,
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

run_human_mobility_stage7 <- function(
  parameters,
  timesteps = 2L,
  export_mixing = list(diag(2)),
  import_mixing = list(diag(2)),
  p_captured = list(matrix(0, nrow = 2, ncol = 2)),
  p_success = 0
) {
  run_metapop_simulation(
    timesteps = timesteps,
    parameters = parameters,
    mixing_tt = 1,
    export_mixing = export_mixing,
    import_mixing = import_mixing,
    p_captured_tt = 1,
    p_captured = p_captured,
    p_success = p_success
  )
}

test_that("documentation-style 2-node explicit mobility setup runs", {
  human_move_probs <- human_mobility_stage7_matrix()
  parameters <- human_mobility_stage7_params(
    list(
      human_mobility_enabled = TRUE,
      human_move_probs = human_move_probs,
      human_trip_duration_type = "fixed",
      human_trip_duration_mean = 1,
      human_mobility_store_diagnostics = TRUE
    ),
    list(
      human_trip_duration_type = "fixed",
      human_trip_duration_mean = 1
    )
  )

  output <- run_human_mobility_stage7(parameters, timesteps = 2L)

  mobility_outputs <- c("humans_present", "visitors_present", "residents_away", "trips_started")
  expect_true(all(mobility_outputs %in% names(output[[1L]])))
  expect_true(all(mobility_outputs %in% names(output[[2L]])))
  expect_equal(dim(attr(output, "OD_started_trips")), c(2L, 2L, 2L))
  expect_equal(dim(attr(output, "OD_active_overnight_stays")), c(2L, 2L, 2L))
  expect_equal(dim(attr(output, "mean_remaining_trip_duration")), c(2L, 2L))
})

test_that("human-carried intervention state remains on the home-node record when current_node changes", {
  parameters <- human_mobility_stage7_params(
    list(
      human_mobility_enabled = TRUE,
      human_move_probs = matrix(c(0, 1, 0, 1), nrow = 2, byrow = TRUE),
      human_trip_duration_type = "fixed",
      human_trip_duration_mean = 1
    ),
    list(
      human_mobility_enabled = TRUE,
      human_move_probs = matrix(c(0, 1, 0, 1), nrow = 2, byrow = TRUE)
    ),
    population_1 = 2L,
    population_2 = 2L
  )
  for (i in seq_along(parameters)) {
    parameters[[i]]$human_mobility_node_index <- as.integer(i)
    parameters[[i]]$human_mobility_n_nodes <- 2L
  }

  variables <- lapply(parameters, create_variables)
  variables[[1L]]$net_time <- individual::IntegerVariable$new(c(10L, -1L))
  variables[[1L]]$spray_time <- individual::IntegerVariable$new(c(20L, -1L))
  variables[[1L]]$drug <- individual::IntegerVariable$new(c(1L, 0L))
  variables[[1L]]$drug_time <- individual::IntegerVariable$new(c(30L, -1L))
  variables[[1L]]$last_pev_timestep <- individual::IntegerVariable$new(c(40L, -1L))
  variables[[1L]]$last_eff_pev_timestep <- individual::IntegerVariable$new(c(50L, -1L))
  variables[[1L]]$pev_profile <- individual::IntegerVariable$new(c(1L, -1L))
  variables[[1L]]$tbv_vaccinated <- individual::DoubleVariable$new(c(60, -1))
  variables[[1L]]$zeta <- individual::DoubleVariable$new(c(1.5, 0.5))
  variables[[1L]]$human_slot_contact_multiplier <- individual::DoubleVariable$new(c(2, 0.75))

  context <- create_human_mobility_context(parameters, variables, timesteps = 1L)
  process <- create_human_mobility_process(
    context = context,
    node_index = 1L,
    renderer = NullRender$new(1)
  )

  process(1L)

  expect_equal(variables[[1L]]$current_node$get_values(), c(2L, 2L))
  expect_equal(variables[[1L]]$net_time$get_values(), c(10L, -1L))
  expect_equal(variables[[1L]]$spray_time$get_values(), c(20L, -1L))
  expect_equal(variables[[1L]]$drug$get_values(), c(1L, 0L))
  expect_equal(variables[[1L]]$drug_time$get_values(), c(30L, -1L))
  expect_equal(variables[[1L]]$last_pev_timestep$get_values(), c(40L, -1L))
  expect_equal(variables[[1L]]$last_eff_pev_timestep$get_values(), c(50L, -1L))
  expect_equal(variables[[1L]]$pev_profile$get_values(), c(1L, -1L))
  expect_equal(variables[[1L]]$tbv_vaccinated$get_values(), c(60, -1))
  expect_equal(variables[[1L]]$zeta$get_values(), c(1.5, 0.5))
  expect_equal(variables[[1L]]$human_slot_contact_multiplier$get_values(), c(2, 0.75))
})

test_that("default mobility output remains minimal and diagnostics require opt-in", {
  parameters <- human_mobility_stage7_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage7_matrix(),
    human_trip_duration_type = "fixed",
    human_trip_duration_mean = 1
  ))

  output <- run_human_mobility_stage7(parameters, timesteps = 1L)

  mobility_outputs <- c("humans_present", "visitors_present", "residents_away", "trips_started")
  diagnostics <- c("OD_started_trips", "OD_active_overnight_stays", "mean_remaining_trip_duration")
  expect_true(all(mobility_outputs %in% names(output[[1L]])))
  expect_false(any(diagnostics %in% names(output[[1L]])))
  expect_null(attr(output, "OD_started_trips"))
  expect_null(attr(output, "OD_active_overnight_stays"))
  expect_null(attr(output, "mean_remaining_trip_duration"))
})

test_that("documented explicit mobility compatibility restrictions remain clear", {
  parameters <- human_mobility_stage7_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage7_matrix()
  ))

  non_native <- parameters
  non_native[[1L]]$native_mosquito_backend <- FALSE
  non_native[[2L]]$native_mosquito_backend <- FALSE
  expect_error(
    run_human_mobility_stage7(non_native, timesteps = 1L),
    "native metapop mosquito backend"
  )

  non_identity <- list(matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE))
  expect_error(
    run_human_mobility_stage7(parameters, timesteps = 1L, export_mixing = non_identity),
    "export_mixing"
  )
  expect_error(
    run_human_mobility_stage7(parameters, timesteps = 1L, import_mixing = non_identity),
    "import_mixing"
  )

  p_captured <- matrix(0, nrow = 2, ncol = 2)
  p_captured[1, 2] <- 0.1
  expect_error(
    run_human_mobility_stage7(parameters, timesteps = 1L, p_captured = list(p_captured)),
    "p_captured"
  )

  rates <- parameters
  rates[[1L]]$human_move_rates <- c(1, 1)
  expect_error(
    human_mobility_validate_parameters(rates, n_nodes = 2L),
    "human_move_rates"
  )

  unsupported_mode <- parameters
  unsupported_mode[[1L]]$human_mobility_mode <- "native"
  expect_error(
    human_mobility_validate_parameters(unsupported_mode, n_nodes = 2L),
    "human_mobility_mode"
  )

  disabled_de_zero <- human_mobility_stage7_params(
    list(de = 0),
    list(de = 0)
  )
  expect_error(
    human_mobility_validate_parameters(disabled_de_zero, n_nodes = 2L),
    NA
  )

  enabled_de_zero <- parameters
  enabled_de_zero[[1L]]$de <- 0
  expect_error(
    human_mobility_validate_parameters(enabled_de_zero, n_nodes = 2L),
    "greater than 0"
  )
})
