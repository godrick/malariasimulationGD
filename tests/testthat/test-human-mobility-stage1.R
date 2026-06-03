human_mobility_stage1_matrix <- function() {
  matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
}

human_mobility_stage1_params <- function(overrides_1 = list(), overrides_2 = list()) {
  base <- list(
    human_population = 6,
    total_M = 20,
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    progress_bar = FALSE
  )
  make_params <- function(overrides) {
    values <- base
    values[names(overrides)] <- overrides
    get_parameters(values)
  }
  list(make_params(overrides_1), make_params(overrides_2))
}

run_human_mobility_stage1 <- function(
  parameters,
  export_mixing = list(diag(2)),
  import_mixing = list(diag(2)),
  p_captured = list(matrix(0, nrow = 2, ncol = 2)),
  p_success = 0
) {
  run_metapop_simulation(
    timesteps = 1,
    parameters = parameters,
    mixing_tt = 1,
    export_mixing = export_mixing,
    import_mixing = import_mixing,
    p_captured_tt = 1,
    p_captured = p_captured,
    p_success = p_success
  )
}

test_that("human mobility disabled follows existing native metapop path", {
  expect_error(
    run_human_mobility_stage1(human_mobility_stage1_params()),
    NA
  )
})

test_that("explicit human mobility accepts one shared row-stochastic matrix", {
  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage1_matrix()
  ))

  expect_error(run_human_mobility_stage1(parameters), NA)
})

test_that("explicit human mobility validates human_move_probs dimensions", {
  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = diag(3)
  ))

  expect_error(
    run_human_mobility_stage1(parameters),
    "n_nodes x n_nodes"
  )
})

test_that("explicit human mobility rejects invalid human_move_probs entries", {
  negative <- human_mobility_stage1_matrix()
  negative[1, 1] <- -0.1
  negative[1, 2] <- 1.1
  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = negative
  ))
  expect_error(run_human_mobility_stage1(parameters), "negative")

  for (bad_value in list(NA_real_, NaN, Inf)) {
    bad_matrix <- human_mobility_stage1_matrix()
    bad_matrix[1, 1] <- bad_value
    parameters <- human_mobility_stage1_params(list(
      human_mobility_enabled = TRUE,
      human_move_probs = bad_matrix
    ))
    expect_error(
      run_human_mobility_stage1(parameters),
      "missing|finite"
    )
  }
})

test_that("explicit human mobility rejects non-row-stochastic human_move_probs", {
  bad_matrix <- human_mobility_stage1_matrix()
  bad_matrix[1, ] <- c(0.7, 0.2)
  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = bad_matrix
  ))

  expect_error(
    run_human_mobility_stage1(parameters),
    "sum to 1"
  )
})

test_that("explicit human mobility verifies multiple non-NULL matrices match", {
  other_matrix <- human_mobility_stage1_matrix()
  other_matrix[1, ] <- c(0.7, 0.3)
  parameters <- human_mobility_stage1_params(
    list(
      human_mobility_enabled = TRUE,
      human_move_probs = human_mobility_stage1_matrix()
    ),
    list(human_move_probs = other_matrix)
  )

  expect_error(
    run_human_mobility_stage1(parameters),
    "identical"
  )
})

test_that("explicit human mobility rejects unsupported mobility parameters", {
  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage1_matrix(),
    human_move_rates = c(1, 1)
  ))
  expect_error(run_human_mobility_stage1(parameters), "human_move_rates")

  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage1_matrix(),
    human_mobility_mode = "native"
  ))
  expect_error(run_human_mobility_stage1(parameters), "human_mobility_mode")
})

test_that("explicit human mobility validates trip-duration settings", {
  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage1_matrix(),
    human_trip_duration_type = "empirical"
  ))
  expect_error(run_human_mobility_stage1(parameters), "human_trip_duration_type")

  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage1_matrix(),
    human_trip_duration_type = "fixed",
    human_trip_duration_mean = 1.5
  ))
  expect_error(run_human_mobility_stage1(parameters), "positive integer")

  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage1_matrix(),
    human_trip_duration_type = "geometric",
    human_trip_duration_mean = 0.5
  ))
  expect_error(run_human_mobility_stage1(parameters), ">= 1")
})

test_that("explicit human mobility rejects non-native metapop backend", {
  parameters <- human_mobility_stage1_params(
    list(
      human_mobility_enabled = TRUE,
      human_move_probs = human_mobility_stage1_matrix(),
      native_mosquito_backend = FALSE
    ),
    list(native_mosquito_backend = FALSE)
  )

  expect_error(
    run_human_mobility_stage1(parameters),
    "native metapop mosquito backend"
  )
})

test_that("explicit human mobility rejects transmission mixing and border capture", {
  parameters <- human_mobility_stage1_params(list(
    human_mobility_enabled = TRUE,
    human_move_probs = human_mobility_stage1_matrix()
  ))
  non_identity <- list(matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE))

  expect_error(
    run_human_mobility_stage1(parameters, export_mixing = non_identity),
    "export_mixing"
  )
  expect_error(
    run_human_mobility_stage1(parameters, import_mixing = non_identity),
    "import_mixing"
  )

  p_captured <- matrix(0, nrow = 2, ncol = 2)
  p_captured[1, 2] <- 0.1
  expect_error(
    run_human_mobility_stage1(parameters, p_captured = list(p_captured)),
    "p_captured"
  )
})
