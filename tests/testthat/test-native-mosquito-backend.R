make_native_test_cube <- function() {
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

test_that("native deterministic mosquito backend parameterises and steps", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    total_M = 100
  ))

  models <- parameterise_mosquito_models(parameters, timesteps = 2)
  solvers <- parameterise_solvers(models, parameters)

  expect_true(inherits(models[[1]], "NativeMosquitoModel"))
  expect_true(inherits(solvers[[1]], "NativeMosquitoSolver"))

  native_mosquito_model_update(
    models[[1]],
    timestep = 0,
    mu = parameters$mum[[1]],
    foim = parameters$init_foim,
    f = parameters$blood_meal_rates[[1]]
  )
  solvers[[1]]$step()

  summary <- solvers[[1]]$get_summary()
  expect_true(all(summary$totals >= 0))
  expect_false(anyNA(summary$totals))
})

test_that("native deterministic backend stays at equilibrium after set_equilibrium", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE
  ))
  parameters <- set_equilibrium(parameters, 5)

  models <- parameterise_mosquito_models(parameters, timesteps = 30)
  solvers <- parameterise_solvers(models, parameters)
  model <- models[[1]]
  solver <- solvers[[1]]

  initial <- solver$get_summary()$totals

  for (t in 0:29) {
    native_mosquito_model_update(
      model,
      timestep = t,
      mu = parameters$mum[[1]],
      foim = parameters$init_foim,
      f = parameters$blood_meal_rates[[1]]
    )
    solver$step()
  }

  final <- solver$get_summary()$totals
  expect_equal(final, initial, tolerance = 1e-5)
})

test_that("native deterministic backend stays at equilibrium with multistage overrides", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    del = 2 / (1 / 3),
    dl = 3 / (1 / 7),
    dpl = 2 / 1,
    me = 0.05,
    ml = 0.15,
    mup = 0.05,
    mum = 0.132,
    beta = 16
  ))
  parameters$native_mosquito_nE <- 2L
  parameters$native_mosquito_nL <- 3L
  parameters$native_mosquito_nP <- 2L
  parameters$native_mosquito_nu <- 1 / (4 / 24)
  parameters <- set_equilibrium(parameters, 8.015726)

  models <- parameterise_mosquito_models(parameters, timesteps = 30)
  solvers <- parameterise_solvers(models, parameters)
  model <- models[[1]]
  solver <- solvers[[1]]

  initial <- solver$get_summary()$totals

  for (t in 0:29) {
    native_mosquito_model_update(
      model,
      timestep = t,
      mu = parameters$mum[[1]],
      foim = parameters$init_foim,
      f = parameters$blood_meal_rates[[1]]
    )
    solver$step()
  }

  final <- solver$get_summary()$totals
  expect_equal(final, initial, tolerance = 1e-5)
})

test_that("native stochastic mosquito backend parameterises and steps", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = TRUE,
    total_M = 100
  ))

  models <- expect_warning(
    parameterise_mosquito_models(parameters, timesteps = 2),
    "count-based tau-leap mosquito engine"
  )
  solvers <- parameterise_solvers(models, parameters)

  expect_true(inherits(models[[1]], "NativeMosquitoModel"))
  expect_true(inherits(solvers[[1]], "NativeMosquitoSolver"))

  set.seed(1)
  native_mosquito_model_update(
    models[[1]],
    timestep = 0,
    mu = parameters$mum[[1]],
    foim = parameters$init_foim,
    f = parameters$blood_meal_rates[[1]]
  )
  solvers[[1]]$step()

  summary <- solvers[[1]]$get_summary()
  expect_true(all(summary$totals >= 0))
  expect_false(anyNA(summary$totals))
})

test_that("native backend rejects incompatible restored state lengths", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    total_M = 100
  ))

  models <- parameterise_mosquito_models(parameters, timesteps = 2)
  solvers <- parameterise_solvers(models, parameters)
  solver <- solvers[[1]]
  state <- solver$get_native_state()

  expect_error(
    solver$set_native_state(state[-1], t = 0),
    "Native mosquito state length mismatch during restore/set"
  )

  expect_error(
    solver$restore_state(
      0,
      list(
        t = 0,
        state = state[-1],
        pending_inputs = native_empty_pending_inputs(1L),
        last_completed_timestep = NULL
      )
    ),
    "Native mosquito state length mismatch during restore/set"
  )
})

test_that("native backend respects cube$c when computing mosquito infection", {
  cube <- make_native_test_cube()
  cube$c <- setNames(c(1, 0), cube$genotypesID)

  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    cube = cube,
    total_M = 0,
    init_foim = 0
  ))

  models <- parameterise_mosquito_models(parameters, timesteps = 1)
  solvers <- parameterise_solvers(models, parameters)
  model <- models[[1]]
  solver <- solvers[[1]]
  idx <- model$shared$index
  state <- solver$get_native_state()
  state[] <- 0
  state[idx$fem_ix[1, 1, 1, 1]] <- 20
  state[idx$fem_ix[2, 2, 1, 1]] <- 20
  state[idx$hS_ix[[1]]] <- 0
  state[idx$hI_ix[[1]]] <- 1
  solver$set_native_state(state, t = 0)

  native_mosquito_model_update(
    model,
    timestep = 0,
    mu = 0,
    foim = 1,
    f = 0
  )
  solver$step()

  state_after <- solver$get_native_state()
  eip_stages <- if (idx$nStages > 2L) 2L:(idx$nStages - 1L) else integer(0)
  ww_exposed <- sum(state_after[idx$fem_ix[1, , eip_stages, 1]])
  hh_exposed <- sum(state_after[idx$fem_ix[2, , eip_stages, 1]])

  expect_gt(ww_exposed, 0)
  expect_equal(hh_exposed, 0, tolerance = 1e-8)
})

test_that("shell human transmission weighting falls back to cube$b", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE
  ))
  cube <- make_native_test_cube()
  cube$b <- setNames(c(0.55, 0), cube$genotypesID)
  parameters$cube <- cube

  fake_solver <- list(
    get_summary = function(node = 1L) {
      list(infectious = setNames(c(10, 20), cube$genotypesID))
    }
  )

  expect_equal(
    calculate_infectious(1, list(fake_solver), variables = NULL, parameters = parameters),
    30
  )

  expect_equal(
    calculate_transmission_infectious(1, list(fake_solver), variables = NULL, parameters = parameters),
    10 * (0.55 / parameters$b0)
  )

  parameters$vector_infectivity_g <- setNames(c(1, 1), cube$genotypesID)
  parameters <- validate_vector_infectivity_g_parameters(parameters)

  expect_equal(
    calculate_infectious(1, list(fake_solver), variables = NULL, parameters = parameters),
    30
  )

  expect_equal(
    calculate_transmission_infectious(1, list(fake_solver), variables = NULL, parameters = parameters),
    30
  )
})

test_that("native backend resolves per-species beta without collapsing to a scalar", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    species = c("sp1", "sp2"),
    species_proportions = c(0.5, 0.5),
    beta = c(sp1 = 11, sp2 = 29),
    blood_meal_rates = c(0.25, 0.5),
    foraging_time = c(0.5, 0.5),
    Q0 = c(0.9, 0.9),
    phi_bednets = c(0.8, 0.8),
    phi_indoors = c(0.85, 0.85),
    mum = c(0.12, 0.18)
  ))

  models <- parameterise_mosquito_models(parameters, timesteps = 2)

  native_mosquito_model_update(
    models[[1]],
    timestep = 0,
    mu = 0.12,
    foim = 0,
    f = 0.25
  )
  native_mosquito_model_update(
    models[[2]],
    timestep = 0,
    mu = 0.18,
    foim = 0,
    f = 0.5
  )

  expect_equal(models[[1]]$species_beta, 11)
  expect_equal(models[[2]]$species_beta, 29)
  expect_equal(
    models[[1]]$shared$pending_inputs$beta[[1]],
    eggs_laid(11, 0.12, 0.25)
  )
  expect_equal(
    models[[2]]$shared$pending_inputs$beta[[1]],
    eggs_laid(29, 0.18, 0.5)
  )
})

test_that("native metapop mosquito backend shares a moving backend across nodes", {
  move_probs <- matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE)
  move_rates <- c(5, 0)

  p1 <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    total_M = 100,
    move_probs = move_probs,
    move_rates = move_rates
  ))
  p2 <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    total_M = 0,
    move_probs = move_probs,
    move_rates = move_rates
  ))

  backend <- parameterise_native_metapop_backends(list(p1, p2), timesteps = 2)

  native_mosquito_model_update(
    backend$models[[1]][[1]],
    timestep = 0,
    mu = p1$mum[[1]],
    foim = p1$init_foim,
    f = p1$blood_meal_rates[[1]]
  )
  native_mosquito_model_update(
    backend$models[[2]][[1]],
    timestep = 0,
    mu = p2$mum[[1]],
    foim = p2$init_foim,
    f = p2$blood_meal_rates[[1]]
  )

  backend$solvers[[1]][[1]]$step()
  backend$solvers[[2]][[1]]$step()

  node2 <- backend$solvers[[2]][[1]]$get_summary()
  expect_gt(sum(node2$totals[c("Sm", "Pm", "Im")]), 0)
})

test_that("native deterministic mosquito movement includes unmated females", {
  move_probs <- matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE)
  move_rates <- c(0.5, 0)

  p1 <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    total_M = 0,
    move_probs = move_probs,
    move_rates = move_rates
  ))
  p2 <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    total_M = 0,
    move_probs = move_probs,
    move_rates = move_rates
  ))

  backend <- parameterise_native_metapop_backends(list(p1, p2), timesteps = 2)
  solver1 <- backend$solvers[[1]][[1]]
  model1 <- backend$models[[1]][[1]]
  idx <- model1$shared$index
  state <- solver1$get_native_state()
  state[] <- 0
  state[idx$unm_ix[1, 1]] <- 100
  solver1$set_native_state(state, t = 0)

  native_mosquito_model_update(
    backend$models[[1]][[1]],
    timestep = 0,
    mu = 0,
    foim = 0,
    f = p1$blood_meal_rates[[1]]
  )
  native_mosquito_model_update(
    backend$models[[2]][[1]],
    timestep = 0,
    mu = 0,
    foim = 0,
    f = p2$blood_meal_rates[[1]]
  )

  solver1$step()
  backend$solvers[[2]][[1]]$step()

  node1 <- backend$solvers[[1]][[1]]$get_summary()
  node2 <- backend$solvers[[2]][[1]]$get_summary()
  expect_lt(unname(node1$unmated[[1]]), 100)
  expect_gt(unname(node2$unmated[[1]]), 0)
})

test_that("native stochastic mosquito movement includes unmated females", {
  move_probs <- matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE)
  move_rates <- c(1, 0)

  p1 <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = TRUE,
    total_M = 0,
    mosquito_tau_step = 1,
    move_probs = move_probs,
    move_rates = move_rates
  ))
  p2 <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = TRUE,
    total_M = 0,
    mosquito_tau_step = 1,
    move_probs = move_probs,
    move_rates = move_rates
  ))

  backend <- expect_warning(
    parameterise_native_metapop_backends(list(p1, p2), timesteps = 2),
    "count-based tau-leap mosquito engine"
  )
  solver1 <- backend$solvers[[1]][[1]]
  model1 <- backend$models[[1]][[1]]
  idx <- model1$shared$index
  state <- solver1$get_native_state()
  state[] <- 0
  state[idx$unm_ix[1, 1]] <- 1000
  solver1$set_native_state(state, t = 0)

  set.seed(1)
  native_mosquito_model_update(
    backend$models[[1]][[1]],
    timestep = 0,
    mu = 0,
    foim = 0,
    f = p1$blood_meal_rates[[1]]
  )
  native_mosquito_model_update(
    backend$models[[2]][[1]],
    timestep = 0,
    mu = 0,
    foim = 0,
    f = p2$blood_meal_rates[[1]]
  )

  solver1$step()
  backend$solvers[[2]][[1]]$step()

  node2 <- backend$solvers[[2]][[1]]$get_summary()
  expect_gt(unname(node2$unmated[[1]]), 0)
})

test_that("native deterministic metapop releases propagate into genotype outputs", {
  cube <- make_native_test_cube()
  base <- get_parameters(list(
    human_population = 50,
    total_M = 20,
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    cube = cube
  ))
  p1 <- set_releases(base, list(
    releasesStart = 5L,
    releasesNumber = 1L,
    releaseCount = 40L,
    releaseSex = "M",
    releaseGenotype = "HH",
    releasesInterval = 0L
  ))
  p2 <- base
  identity_mix <- list(diag(2))

  res <- run_metapop_simulation(
    timesteps = 20,
    parameters = list(p1, p2),
    mixing_tt = 1,
    export_mixing = identity_mix,
    import_mixing = identity_mix,
    p_captured_tt = 1,
    p_captured = list(matrix(0, 2, 2)),
    p_success = 0
  )

  release_node_female <- attr(res[[1]], "mosquito_genotype_counts_female")
  control_node_female <- attr(res[[2]], "mosquito_genotype_counts_female")

  expect_false(is.null(release_node_female))
  expect_false(is.null(control_node_female))
  expect_gt(max(release_node_female[, "HH"]), 0)
  expect_equal(max(control_node_female[, "HH"]), 0)
})

test_that("legacy mosquito event wiring is validated against the selected backend", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = TRUE
  ))
  events <- create_events(parameters)

  expect_false(legacy_individual_mosquito_events_enabled(parameters, events))

  bad_events <- events
  bad_events$mosquito_infection <- individual::TargetedEvent$new(parameters$mosquito_limit)
  bad_events$mosquito_death <- individual::TargetedEvent$new(parameters$mosquito_limit)

  expect_error(
    legacy_individual_mosquito_events_enabled(parameters, bad_events),
    "inconsistent with the selected mosquito backend"
  )
})

test_that("native carrying-capacity lookup matches direct daily calculations", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    model_seasonality = TRUE,
    total_M = 100
  ))
  parameters <- set_carrying_capacity(
    parameters,
    timesteps = c(2, 5),
    carrying_capacity_scalers = matrix(c(1.2, 0.8), ncol = 1)
  )

  lookup <- native_build_carrying_capacity_lookup(list(parameters), species_i = 1, timesteps = 10)

  for (tt in 0:10) {
    expect_equal(
      native_carrying_capacity_at(lookup, tt),
      vapply(
        list(parameters),
        function(p) native_effective_carrying_capacity(p, 1, tt),
        numeric(1)
      )
    )
  }

  for (tt in c(120L, 365L, 5841L)) {
    expect_equal(
      native_carrying_capacity_at(lookup, tt),
      vapply(
        list(parameters),
        function(p) native_effective_carrying_capacity(p, 1, tt),
        numeric(1)
      )
    )
  }
})

test_that("native backend applies restart time offsets to runtime timesteps", {
  parameters <- get_parameters(list(
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    total_M = 100
  ))
  parameters$vector_control_time_offset <- 20L

  models <- parameterise_mosquito_models(parameters, timesteps = 2)
  native_mosquito_model_update(
    models[[1]],
    timestep = 1L,
    mu = parameters$mum[[1]],
    foim = parameters$init_foim,
    f = parameters$blood_meal_rates[[1]]
  )

  expect_equal(models[[1]]$shared$time_offset, 20L)
  expect_equal(models[[1]]$shared$pending_inputs$timestep, 21L)
})
