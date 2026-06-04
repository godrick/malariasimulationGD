human_mobility_stage6_cube <- function() {
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

human_mobility_stage6_params <- function(
  population_1 = 2L,
  population_2 = 2L,
  overrides_1 = list(),
  overrides_2 = list(),
  mobility = TRUE
) {
  make_params <- function(population, overrides) {
    values <- list(
      human_population = population,
      native_mosquito_backend = TRUE,
      individual_mosquitoes = FALSE,
      total_M = 20,
      average_age = 1e12,
      progress_bar = FALSE
    )
    if (isTRUE(mobility)) {
      values$human_mobility_enabled <- TRUE
      values$human_move_probs <- diag(2)
    }
    values[names(overrides)] <- overrides
    get_parameters(values)
  }

  parameters <- list(make_params(population_1, overrides_1), make_params(population_2, overrides_2))
  for (i in seq_along(parameters)) {
    parameters[[i]]$human_mobility_node_index <- as.integer(i)
  }
  parameters
}

human_mobility_stage6_run <- function(parameters, timesteps = 5L) {
  run_metapop_simulation(
    timesteps = timesteps,
    parameters = parameters,
    mixing_tt = 1,
    export_mixing = list(diag(length(parameters))),
    import_mixing = list(diag(length(parameters))),
    p_captured_tt = 1,
    p_captured = list(matrix(0, nrow = length(parameters), ncol = length(parameters))),
    p_success = 0,
    restore_random_state = TRUE
  )
}

human_mobility_stage6_set_node <- function(
  variables,
  infectivity,
  current_node = NULL,
  birth = rep(-1000, length(infectivity)),
  zeta = rep(1, length(infectivity))
) {
  variables$infectivity <- individual::DoubleVariable$new(infectivity)
  variables$birth <- individual::DoubleVariable$new(birth)
  variables$zeta <- individual::DoubleVariable$new(zeta)
  if (!is.null(current_node)) {
    variables$current_node <- individual::IntegerVariable$new(as.integer(current_node))
  }
  variables
}

human_mobility_stage6_solver <- function(step_mock, infectious = 0) {
  list(
    step = step_mock,
    get_states = function() numeric(0),
    get_summary = function(node = 1L) {
      list(
        infectious = infectious,
        female = infectious,
        male = infectious
      )
    }
  )
}

test_that("native mobility process order preserves original solver stepping position", {
  parameters <- human_mobility_stage6_params()
  variables <- lapply(parameters, create_variables)
  mobility_context <- create_human_mobility_context(parameters, variables, timesteps = 1L)
  infectivity_context <- create_human_infectivity_lag_context(parameters, variables)
  step_mock <- mockery::mock()
  solver <- human_mobility_stage6_solver(step_mock)

  processes <- create_processes(
    renderer = NullRender$new(1),
    variables = variables[[1L]],
    events = create_events(parameters[[1L]]),
    parameters = parameters[[1L]],
    models = list(list()),
    solvers = list(solver),
    lagged_eir = list(LaggedValue$new(2, 0)),
    lagged_infectivity = LaggedValue$new(2, 0),
    lagged_transmission_eir = list(LaggedValue$new(2, 0)),
    human_mobility_context = mobility_context,
    human_infectivity_lag_context = infectivity_context,
    enable_rendering = FALSE
  )

  process_names <- names(processes)
  expect_lt(match("mosquito_release_process", process_names), match("human_mobility_process", process_names))
  expect_lt(match("human_mobility_process", process_names), match("human_infectivity_lag_process", process_names))
  expect_lt(match("human_infectivity_lag_process", process_names), match("biting_process", process_names))
  expect_lt(match("biting_process", process_names), match("solver_process", process_names))

  processes$human_mobility_process(1L)
  processes$human_infectivity_lag_process(1L)
  expect_length(mockery::mock_args(step_mock), 0L)

  processes$solver_process(1L)
  expect_length(mockery::mock_args(step_mock), 1L)
})

test_that("native biting handoff does not step solvers and passes mobility-aware scalar FOIM", {
  parameters <- human_mobility_stage6_params()
  variables <- lapply(parameters, create_variables)
  variables[[1L]] <- human_mobility_stage6_set_node(
    variables[[1L]],
    infectivity = c(1, 0),
    current_node = c(2L, 1L)
  )
  variables[[2L]] <- human_mobility_stage6_set_node(
    variables[[2L]],
    infectivity = c(0, 0),
    current_node = c(2L, 2L)
  )
  infectivity_context <- create_human_infectivity_lag_context(parameters, variables)
  human_infectivity_lag_record_node(infectivity_context, 1L, timestep = 1L)
  human_infectivity_lag_record_node(infectivity_context, 2L, timestep = 1L)

  step_mock <- mockery::mock()
  solver <- human_mobility_stage6_solver(step_mock, infectious = 0)
  native_update_mock <- mockery::mock()
  mockery::stub(simulate_bites, "native_mosquito_model_update", native_update_mock)
  mockery::stub(simulate_bites, ".human_blood_meal_rate", function(...) 2)

  simulate_bites(
    renderer = NullRender$new(2),
    solvers = list(solver),
    models = list(list()),
    variables = variables[[2L]],
    events = create_events(parameters[[2L]]),
    age = rep(1000, parameters[[2L]]$human_population),
    parameters = parameters[[2L]],
    timestep = 2L,
    lagged_infectivity = LaggedValue$new(3, 0),
    lagged_eir = list(LaggedValue$new(3, 0)),
    mixing_fn = NULL,
    mixing_index = 2L,
    lagged_transmission_eir = list(LaggedValue$new(3, 0)),
    human_exposure_lag_context = NULL,
    human_infectivity_lag_context = infectivity_context
  )

  expect_length(mockery::mock_args(step_mock), 0L)
  native_args <- mockery::mock_args(native_update_mock)[[1L]]
  foim <- native_args[[4L]]
  expect_length(foim, 1L)
  expect_equal(foim, 2 * (1 / 3), tolerance = 1e-12)
})

test_that("native update rejects non-scalar or invalid FOIM values", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    total_M = 20
  ))

  expect_invalid <- function(foim) {
    models <- parameterise_mosquito_models(parameters, timesteps = 1L)
    expect_error(
      native_mosquito_model_update(
        models[[1L]],
        timestep = 1L,
        mu = parameters$mum[[1L]],
        foim = foim,
        f = parameters$blood_meal_rates[[1L]]
      ),
      "single nonnegative finite"
    )
  }

  expect_invalid(c(0, 1))
  expect_invalid(NA_real_)
  expect_invalid(Inf)
  expect_invalid(-0.1)

  models <- parameterise_mosquito_models(parameters, timesteps = 1L)
  native_mosquito_model_update(
    models[[1L]],
    timestep = 1L,
    mu = parameters$mum[[1L]],
    foim = 0.25,
    f = parameters$blood_meal_rates[[1L]]
  )
  expect_equal(models[[1L]]$shared$pending_inputs$foim[[1L]], 0.25)
  expect_length(models[[1L]]$shared$pending_inputs$foim[[1L]], 1L)
})

test_that("identity explicit mobility leaves deterministic native genotype summaries unchanged", {
  cube <- human_mobility_stage6_cube()
  disabled <- human_mobility_stage6_params(
    overrides_1 = list(cube = cube, init_foim = 0),
    overrides_2 = list(cube = cube, init_foim = 0),
    mobility = FALSE
  )
  enabled <- human_mobility_stage6_params(
    overrides_1 = list(cube = cube, init_foim = 0),
    overrides_2 = list(cube = cube, init_foim = 0),
    mobility = TRUE
  )

  disabled_out <- human_mobility_stage6_run(disabled, timesteps = 5L)
  enabled_out <- human_mobility_stage6_run(enabled, timesteps = 5L)

  for (node in seq_along(disabled_out)) {
    expect_equal(
      attr(enabled_out[[node]], "mosquito_genotype_counts_female"),
      attr(disabled_out[[node]], "mosquito_genotype_counts_female"),
      tolerance = 1e-10
    )
    expect_equal(
      attr(enabled_out[[node]], "mosquito_genotype_counts_male"),
      attr(disabled_out[[node]], "mosquito_genotype_counts_male"),
      tolerance = 1e-10
    )
    expect_equal(
      attr(enabled_out[[node]], "mosquito_genotype_V"),
      attr(disabled_out[[node]], "mosquito_genotype_V"),
      tolerance = 1e-10
    )
  }
})

test_that("native genotype and sex releases still run with explicit mobility active", {
  cube <- human_mobility_stage6_cube()
  parameters <- human_mobility_stage6_params(
    overrides_1 = list(cube = cube),
    overrides_2 = list(cube = cube)
  )
  parameters[[1L]] <- set_releases(parameters[[1L]], list(
    releasesStart = 1L,
    releasesNumber = 1L,
    releaseCount = 20L,
    releaseSex = "both",
    releaseGenotype = "HH",
    releasesInterval = 0L
  ))
  backend <- parameterise_native_metapop_backends(parameters, timesteps = 2L)
  variables <- lapply(parameters, create_variables)
  release_process <- create_mosquito_release_process(
    solvers = backend$solvers[[1L]],
    models = backend$models[[1L]],
    variables = variables[[1L]],
    events = create_events(parameters[[1L]]),
    parameters = parameters[[1L]]
  )

  release_process(1L)
  summary <- backend$solvers[[1L]][[1L]]$get_summary(node = 1L)

  expect_gt(unname(summary$male["HH"]), 0)
  expect_gt(unname(summary$female["HH"]), 0)
})

test_that("native mosquito movement still runs with explicit mobility active", {
  move_probs <- matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE)
  move_rates <- c(0.5, 0)
  parameters <- human_mobility_stage6_params(
    overrides_1 = list(total_M = 0, move_probs = move_probs, move_rates = move_rates),
    overrides_2 = list(total_M = 0, move_probs = move_probs, move_rates = move_rates)
  )
  backend <- parameterise_native_metapop_backends(parameters, timesteps = 2L)
  solver1 <- backend$solvers[[1L]][[1L]]
  idx <- backend$models[[1L]][[1L]]$shared$index
  state <- solver1$get_native_state()
  state[] <- 0
  state[idx$unm_ix[1L, 1L]] <- 100
  solver1$set_native_state(state, t = 0)

  native_mosquito_model_update(
    backend$models[[1L]][[1L]],
    timestep = 0L,
    mu = 0,
    foim = 0,
    f = parameters[[1L]]$blood_meal_rates[[1L]]
  )
  native_mosquito_model_update(
    backend$models[[2L]][[1L]],
    timestep = 0L,
    mu = 0,
    foim = 0,
    f = parameters[[2L]]$blood_meal_rates[[1L]]
  )
  solver1$step()
  backend$solvers[[2L]][[1L]]$step()

  expect_gt(unname(backend$solvers[[2L]][[1L]]$get_summary()$unmated[[1L]]), 0)
})

test_that("native EIP progression still runs with explicit mobility active", {
  cube <- human_mobility_stage6_cube()
  parameters <- human_mobility_stage6_params(
    overrides_1 = list(cube = cube, total_M = 0, init_foim = 0),
    overrides_2 = list(cube = cube, total_M = 0, init_foim = 0)
  )
  backend <- parameterise_native_metapop_backends(parameters, timesteps = 2L)
  solver1 <- backend$solvers[[1L]][[1L]]
  model1 <- backend$models[[1L]][[1L]]
  idx <- model1$shared$index
  state <- solver1$get_native_state()
  state[] <- 0
  state[idx$fem_ix[1L, 1L, 1L, 1L]] <- 20
  state[idx$hS_ix[[1L]]] <- 0
  state[idx$hI_ix[[1L]]] <- 1
  solver1$set_native_state(state, t = 0)

  native_mosquito_model_update(
    model1,
    timestep = 0L,
    mu = 0,
    foim = 1,
    f = 0
  )
  native_mosquito_model_update(
    backend$models[[2L]][[1L]],
    timestep = 0L,
    mu = 0,
    foim = 0,
    f = 0
  )
  solver1$step()
  backend$solvers[[2L]][[1L]]$step()

  state_after <- solver1$get_native_state()
  eip_stages <- if (idx$nStages > 2L) 2L:(idx$nStages - 1L) else integer(0)
  expect_gt(length(eip_stages), 0L)
  expect_gt(sum(state_after[idx$fem_ix[1L, , eip_stages, 1L]]), 0)
})

test_that("native C++ human movement remains disabled in constructor calls", {
  source <- readLines(testthat::test_path("../../R/native_mosquito_backend.R"))
  has_hmove_lines <- grep("has_hmove\\s*=", source, value = TRUE)

  expect_length(has_hmove_lines, 2L)
  expect_true(all(grepl("has_hmove\\s*=\\s*FALSE", has_hmove_lines)))
  expect_false(any(grepl("has_hmove\\s*=\\s*TRUE", source)))
})
