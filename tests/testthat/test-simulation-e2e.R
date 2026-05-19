test_that('Simulation runs for a few timesteps', {
  sim <- run_simulation(100)
  expect_equal(nrow(sim), 100)
})

test_that('run_metapop_simulation fails with incorrect mixing matrix', {
  population <- 4
  timestep <- 5
  parameters <- get_parameters(list(human_population = population))
  paramlist <- list(parameters, parameters)
  # incorrect params
  mixing <- matrix(c(1, 1), nrow = 1, ncol = 2)
  p_captured <- diag(nrow=2)
  expect_error(
    run_metapop_simulation(
      timesteps,
      parameters,
      NULL,
      mixing_tt = 1,
      export_mixing = list(mixing),
      import_mixing = list(mixing),
      p_captured_tt = 1,
      p_captured = list(diag(nrow=2)),
      p_success = 1
    )
  )
})

test_that('run_metapop_simulation integrates two models correctly', {
  population <- 4
  timesteps <- 5
  parameters <- get_parameters(list(human_population = population))
  parametersets <- list(parameters, parameters)
  mixing <- diag(nrow = 2)
  p_captured <- 1 - diag(nrow = 2)

  outputs <- run_metapop_simulation(
    timesteps,
    parametersets,
    NULL,
    mixing_tt = 1,
    export_mixing = list(mixing),
    import_mixing = list(mixing),
    p_captured_tt = 1,
    p_captured = list(p_captured),
    p_success = 1
  )
  expect_equal(length(outputs), 2)
  expect_equal(nrow(outputs[[1]]), 5)
  expect_equal(nrow(outputs[[2]]), 5)
})

test_that('run_metapop_simulation populates incidence defaults for every node', {
  timesteps <- 5
  parameters <- get_parameters(list(
    human_population = 25,
    incidence_rendering_min_ages = 182,
    incidence_rendering_max_ages = 3650,
    clinical_incidence_rendering_min_ages = 182,
    clinical_incidence_rendering_max_ages = 3650
  ))
  parametersets <- list(parameters, parameters)
  mixing <- diag(nrow = 2)

  outputs <- run_metapop_simulation(
    timesteps,
    parametersets,
    NULL,
    mixing_tt = 1,
    export_mixing = list(mixing),
    import_mixing = list(mixing),
    p_captured_tt = 1,
    p_captured = list(matrix(0, nrow = 2, ncol = 2)),
    p_success = 1
  )

  for (out in outputs) {
    expect_true("n_infections" %in% names(out))
    expect_true("n_age_182_3650" %in% names(out))
    expect_true("n_inc_182_3650" %in% names(out))
    expect_true("n_inc_clinical_182_3650" %in% names(out))
    expect_false(anyNA(out$n_infections))
    expect_false(anyNA(out$n_age_182_3650))
    expect_false(anyNA(out$n_inc_182_3650))
    expect_false(anyNA(out$n_inc_clinical_182_3650))
  }
})

test_that("run_metapop_simulation can skip rendering while preserving state", {
  old_warning_state <- isTRUE(native_backend_warning_state$tau_leap_backend)
  native_backend_warning_state$tau_leap_backend <- TRUE
  on.exit({
    native_backend_warning_state$tau_leap_backend <- old_warning_state
  }, add = TRUE)

  timesteps <- 6
  parameters <- get_parameters(list(
    human_population = 20,
    total_M = 25,
    native_mosquito_backend = TRUE,
    individual_mosquitoes = TRUE,
    progress_bar = FALSE
  ))
  parametersets <- list(parameters, parameters)
  mixing <- diag(nrow = 2)
  p_captured <- matrix(0, nrow = 2, ncol = 2)

  set.seed(123)
  full <- run_metapop_simulation(
    timesteps = timesteps,
    parameters = parametersets,
    mixing_tt = 1,
    export_mixing = list(mixing),
    import_mixing = list(mixing),
    p_captured_tt = 1,
    p_captured = list(p_captured),
    p_success = 0,
    return_state = TRUE,
    return_summary = TRUE
  )

  set.seed(123)
  lean <- run_metapop_simulation(
    timesteps = timesteps,
    parameters = parametersets,
    mixing_tt = 1,
    export_mixing = list(mixing),
    import_mixing = list(mixing),
    p_captured_tt = 1,
    p_captured = list(p_captured),
    p_success = 0,
    return_state = TRUE,
    render_output = FALSE,
    return_summary = TRUE
  )

  expect_null(lean$data)
  expect_equal(lean$summary$total_M_by_node, full$summary$total_M_by_node)

  set.seed(321)
  resumed_full <- run_metapop_simulation(
    timesteps = 10,
    parameters = parametersets,
    mixing_tt = 1,
    export_mixing = list(mixing),
    import_mixing = list(mixing),
    p_captured_tt = 1,
    p_captured = list(p_captured),
    p_success = 0,
    initial_state = full$state,
    restore_random_state = TRUE
  )

  set.seed(321)
  resumed_lean <- run_metapop_simulation(
    timesteps = 10,
    parameters = parametersets,
    mixing_tt = 1,
    export_mixing = list(mixing),
    import_mixing = list(mixing),
    p_captured_tt = 1,
    p_captured = list(p_captured),
    p_success = 0,
    initial_state = lean$state,
    restore_random_state = TRUE
  )

  expect_equal(resumed_lean, resumed_full)
})
