human_mobility_stage8_cube <- function() {
  G <- 2L
  g <- c("WW", "HH")
  ih <- array(0, dim = c(G, G, G), dimnames = list(g, g, g))
  ih["WW", "WW", "WW"] <- 1
  ih["WW", "HH", "HH"] <- 1
  ih["HH", "WW", "HH"] <- 1
  ih["HH", "HH", "HH"] <- 1

  list(
    ih = ih,
    tau = array(1, dim = c(G, G, G), dimnames = list(g, g, g)),
    eta = matrix(1, nrow = G, ncol = G, dimnames = list(g, g)),
    b = setNames(c(1, 1), g),
    c = setNames(c(1, 1), g),
    phi = setNames(c(0.5, 0.5), g),
    omega = setNames(c(1, 1), g),
    xiF = setNames(c(1, 1), g),
    xiM = setNames(c(1, 1), g),
    s = setNames(c(1, 1), g),
    genotypesID = g,
    wildType = "WW"
  )
}

human_mobility_stage8_matrix <- function(identity = FALSE) {
  if (isTRUE(identity)) {
    return(diag(2))
  }
  matrix(c(0.25, 0.75, 0.2, 0.8), nrow = 2, byrow = TRUE)
}

human_mobility_stage8_params <- function(
  overrides_1 = list(),
  overrides_2 = list(),
  population_1 = 4L,
  population_2 = 4L,
  mobility = TRUE,
  move_probs = human_mobility_stage8_matrix()
) {
  make_params <- function(population, overrides) {
    values <- list(
      human_population = population,
      total_M = 20,
      native_mosquito_backend = TRUE,
      individual_mosquitoes = FALSE,
      average_age = 1e12,
      progress_bar = FALSE,
      de = 1,
      delay_gam = 1
    )
    if (isTRUE(mobility)) {
      values$human_mobility_enabled <- TRUE
      values$human_move_probs <- move_probs
      values$human_trip_duration_type <- "fixed"
      values$human_trip_duration_mean <- 1
    }
    values[names(overrides)] <- overrides
    get_parameters(values)
  }

  parameters <- list(make_params(population_1, overrides_1), make_params(population_2, overrides_2))
  for (i in seq_along(parameters)) {
    parameters[[i]]$human_mobility_node_index <- as.integer(i)
    parameters[[i]]$human_mobility_n_nodes <- length(parameters)
  }
  parameters
}

human_mobility_stage8_run <- function(parameters, timesteps = 4L, return_state = FALSE, initial_state = NULL) {
  run_metapop_simulation(
    timesteps = timesteps,
    parameters = parameters,
    mixing_tt = 1,
    export_mixing = list(diag(length(parameters))),
    import_mixing = list(diag(length(parameters))),
    p_captured_tt = 1,
    p_captured = list(matrix(0, nrow = length(parameters), ncol = length(parameters))),
    p_success = 0,
    return_state = return_state,
    initial_state = initial_state,
    restore_random_state = TRUE
  )
}

human_mobility_stage8_variable_state <- function(state, variable) {
  matches <- grep(paste0("^", variable, "$"), names(state$individual$variables))
  state$individual$variables[matches]
}

test_that("mobility diagnostics allocate only when requested and rendering is enabled", {
  parameters <- human_mobility_stage8_params()
  variables <- lapply(parameters, create_variables)

  context <- create_human_mobility_context(parameters, variables, timesteps = 3L)
  expect_false(exists("OD_started_trips", envir = context, inherits = FALSE))
  expect_false(exists("OD_active_overnight_stays", envir = context, inherits = FALSE))
  expect_false(exists("mean_remaining_trip_duration_history", envir = context, inherits = FALSE))

  diagnostic_parameters <- parameters
  diagnostic_parameters[[1L]]$human_mobility_store_diagnostics <- TRUE
  diagnostic_context <- create_human_mobility_context(
    diagnostic_parameters,
    variables,
    timesteps = 3L
  )
  expect_equal(dim(diagnostic_context$OD_started_trips), c(3L, 2L, 2L))
  expect_equal(dim(diagnostic_context$OD_active_overnight_stays), c(3L, 2L, 2L))
  expect_equal(dim(diagnostic_context$mean_remaining_trip_duration_history), c(3L, 2L))

  lean_context <- create_human_mobility_context(
    diagnostic_parameters,
    variables,
    timesteps = 3L,
    render_output = FALSE
  )
  expect_false(exists("OD_started_trips", envir = lean_context, inherits = FALSE))
  expect_false(exists("OD_active_overnight_stays", envir = lean_context, inherits = FALSE))
  expect_false(exists("mean_remaining_trip_duration_history", envir = lean_context, inherits = FALSE))
})

test_that("default mobility outputs exclude large diagnostics and per-human histories", {
  parameters <- human_mobility_stage8_params()
  output <- human_mobility_stage8_run(parameters, timesteps = 2L)

  default_columns <- c("humans_present", "visitors_present", "residents_away", "trips_started")
  forbidden_columns <- c(
    "OD_started_trips",
    "OD_active_overnight_stays",
    "mean_remaining_trip_duration",
    "traveller_experienced_EIR",
    "visitor_FOIM",
    "visitor_FOIM_proportion",
    "human_exposure_lag",
    "human_infectivity_lag"
  )
  forbidden_attributes <- c(
    forbidden_columns,
    "traveller_experienced_EIR",
    "per_human_exposure_history",
    "per_human_infectivity_history"
  )

  for (node_output in output) {
    expect_true(all(default_columns %in% names(node_output)))
    expect_false(any(forbidden_columns %in% names(node_output)))
  }
  expect_false(any(forbidden_attributes %in% names(attributes(output))))
})

test_that("per-human lag buffers remain bounded by configured lag history", {
  exposure_buffer <- HumanExposureLagBuffer$new(
    max_lag = 2,
    default_exposure = c(0, 0),
    default_weighted_exposure = c(0, 0)
  )
  infectivity_buffer <- HumanInfectivityLagBuffer$new(
    max_lag = 2,
    default_infectivity = c(0, 0)
  )

  for (timestep in seq_len(10)) {
    exposure_buffer$save(
      timestep,
      exposure = c(timestep, timestep + 1),
      weighted_exposure = c(timestep + 2, timestep + 3)
    )
    infectivity_buffer$save(timestep, infectivity = c(timestep / 10, timestep / 20))
  }

  exposure_state <- exposure_buffer$save_state()
  infectivity_state <- infectivity_buffer$save_state()

  expect_lte(length(exposure_state$timesteps), 4L)
  expect_lte(nrow(exposure_state$exposure), 4L)
  expect_lte(nrow(exposure_state$weighted_exposure), 4L)
  expect_equal(exposure_state$timesteps, 7:10)

  expect_lte(length(infectivity_state$timesteps), 4L)
  expect_lte(nrow(infectivity_state$infectivity), 4L)
  expect_equal(infectivity_state$timesteps, 7:10)
})

test_that("checkpoint resume preserves mobility outputs, state, and lag buffers", {
  parameters <- human_mobility_stage8_params(list(human_mobility_store_diagnostics = TRUE))

  set.seed(123)
  full <- human_mobility_stage8_run(parameters, timesteps = 6L, return_state = TRUE)
  set.seed(123)
  checkpoint <- human_mobility_stage8_run(parameters, timesteps = 3L, return_state = TRUE)
  set.seed(999)
  resumed <- human_mobility_stage8_run(
    parameters,
    timesteps = 6L,
    return_state = TRUE,
    initial_state = checkpoint$state
  )

  output_columns <- c(
    "humans_present",
    "visitors_present",
    "residents_away",
    "trips_started",
    "EIR_gamb",
    "FOIM_gamb",
    "n_infections"
  )
  for (node in seq_along(full$data)) {
    common_columns <- intersect(output_columns, names(full$data[[node]]))
    expect_equal(
      resumed$data[[node]][common_columns],
      full$data[[node]][4:6, common_columns, drop = FALSE],
      ignore_attr = TRUE,
      tolerance = 1e-10
    )
  }

  for (variable in c("current_node", "travel_destination", "travel_remaining_nights", "is_travelling")) {
    expect_equal(
      human_mobility_stage8_variable_state(resumed$state, variable),
      human_mobility_stage8_variable_state(full$state, variable),
      ignore_attr = TRUE
    )
  }

  expect_equal(
    resumed$state$malariasimulationGD$human_exposure_lag,
    full$state$malariasimulationGD$human_exposure_lag,
    ignore_attr = TRUE,
    tolerance = 1e-10
  )
  expect_equal(
    resumed$state$malariasimulationGD$human_infectivity_lag,
    full$state$malariasimulationGD$human_infectivity_lag,
    ignore_attr = TRUE,
    tolerance = 1e-10
  )
  expect_equal(
    attr(resumed$data, "OD_started_trips"),
    attr(full$data, "OD_started_trips")[4:6, , , drop = FALSE],
    ignore_attr = TRUE
  )
})

test_that("combined old mobility checkpoint without lag buffers is tolerated", {
  parameters <- human_mobility_stage8_params()

  set.seed(123)
  checkpoint <- human_mobility_stage8_run(parameters, timesteps = 2L, return_state = TRUE)
  checkpoint$state$malariasimulationGD$human_exposure_lag <- NULL
  checkpoint$state$malariasimulationGD$human_infectivity_lag <- NULL

  expect_silent(human_mobility_stage8_run(
    parameters,
    timesteps = 4L,
    return_state = TRUE,
    initial_state = checkpoint$state
  ))
})

test_that("native genotype releases, movement, and EIP run with mobility without default large outputs", {
  cube <- human_mobility_stage8_cube()
  mosquito_move_probs <- matrix(c(0.9, 0.1, 0.05, 0.95), nrow = 2, byrow = TRUE)
  mosquito_move_rates <- c(0.1, 0.05)
  parameters <- human_mobility_stage8_params(
    overrides_1 = list(
      cube = cube,
      total_M = 30,
      init_foim = 0.05,
      move_probs = mosquito_move_probs,
      move_rates = mosquito_move_rates
    ),
    overrides_2 = list(
      cube = cube,
      total_M = 30,
      init_foim = 0.05,
      move_probs = mosquito_move_probs,
      move_rates = mosquito_move_rates
    ),
    population_1 = 3L,
    population_2 = 3L
  )
  parameters[[1L]] <- set_releases(parameters[[1L]], list(
    releasesStart = 1L,
    releasesNumber = 1L,
    releaseCount = 12L,
    releaseSex = "both",
    releaseGenotype = "HH",
    releasesInterval = 0L
  ))

  output <- human_mobility_stage8_run(parameters, timesteps = 4L)

  for (node in seq_along(output)) {
    expect_true(all(c("humans_present", "visitors_present", "residents_away", "trips_started") %in% names(output[[node]])))
    expect_false("OD_started_trips" %in% names(output[[node]]))
    expect_false("human_exposure_lag" %in% names(output[[node]]))
    expect_false(is.null(attr(output[[node]], "mosquito_genotype_counts_female")))
    expect_false(is.null(attr(output[[node]], "mosquito_genotype_counts_male")))
    expect_false(is.null(attr(output[[node]], "mosquito_genotype_V")))
  }
  expect_false(is.null(attr(output[[1L]], "mosquito_release_schedule")))
  expect_null(attr(output, "OD_started_trips"))
})

test_that("non-identity mobility shifts exposure and FOIM signals to destination nodes", {
  parameters <- human_mobility_stage8_params(
    population_1 = 2L,
    population_2 = 2L,
    move_probs = matrix(c(0, 1, 0, 1), nrow = 2, byrow = TRUE)
  )
  variables <- lapply(parameters, create_variables)
  variables[[1L]]$infectivity <- individual::DoubleVariable$new(c(1, 1))
  variables[[2L]]$infectivity <- individual::DoubleVariable$new(c(0, 0))
  variables[[1L]]$zeta <- individual::DoubleVariable$new(c(1, 1))
  variables[[2L]]$zeta <- individual::DoubleVariable$new(c(1, 1))

  mobility_context <- create_human_mobility_context(parameters, variables, timesteps = 1L)
  human_mobility_update_context(mobility_context, timestep = 1L)

  expect_equal(variables[[1L]]$current_node$get_values(), c(2L, 2L))
  expect_equal(colSums(mobility_context$active_od), c(0L, 4L))

  exposure_context <- new.env(parent = emptyenv())
  exposure_context$n_nodes <- 2L
  exposure_context$variables <- variables
  exposure_context$weighted_active <- FALSE
  exposure_context$reported_timestep <- NA_real_
  exposure_context$reported <- rep(FALSE, 2L)
  exposure_context$node_exposure <- rep(NA_real_, 2L)
  exposure_context$node_weighted_exposure <- NULL
  exposure_context$buffers <- list(
    HumanExposureLagBuffer$new(max_lag = 2, default_exposure = c(0, 0)),
    HumanExposureLagBuffer$new(max_lag = 2, default_exposure = c(0, 0))
  )
  human_exposure_lag_record_node(exposure_context, node_index = 1L, timestep = 1L, exposure = 10)
  human_exposure_lag_record_node(exposure_context, node_index = 2L, timestep = 1L, exposure = 100)
  expect_equal(exposure_context$buffers[[1L]]$get(1L), c(100, 100))

  infectivity_context <- create_human_infectivity_lag_context(parameters, variables)
  human_infectivity_lag_record_node(infectivity_context, 1L, timestep = 1L)
  human_infectivity_lag_record_node(infectivity_context, 2L, timestep = 1L)
  reservoir <- human_infectivity_lag_get_reservoir(infectivity_context, timestep = 2L)

  expect_equal(reservoir[[1L]], 0)
  expect_gt(reservoir[[2L]], 0)
})
