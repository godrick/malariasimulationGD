test_that('create_variables allows empty species', {
  params <- get_parameters(list(
    individual_mosquitoes = TRUE
  ))
  params <- set_species(
    params,
    species=list(
      gamb_params,
      arab_params,
      fun_params
    ),
    proportions=c(1,0,0)
  )
  variables <- create_variables(params)
  expect_equal(variables$species$get_size_of('gamb'), params$mosquito_limit)
})

test_that('create_variables allows multiple species', {
  params <- get_parameters(list(
    individual_mosquitoes = TRUE
  ))
  params <- set_species(
    params,
    species=list(
      gamb_params,
      arab_params
    ),
    proportions=c(.9, .1)
  )
  variables <- create_variables(params)
  expect_equal(
    variables$species$get_size_of('arab'),
    params$total_M * .1
  )
})

test_that('create_variables allows multiple species w different total_M', {
  params <- get_parameters(list(
    individual_mosquitoes = TRUE
  ))
  params <- set_species(
    params,
    species=list(
      gamb_params,
      arab_params
    ),
    proportions=c(.9, .1)
  )
  params <- parameterise_total_M(params, 1000)
  variables <- create_variables(params)
  expect_equal(
    variables$species$get_size_of('arab'),
    params$total_M * .1
  )
})

test_that('stationary human initializer shifts saved timestep variables back to day 0', {
  params <- get_parameters()
  fake_library <- list(
    reference_timestep = c(40L, 50L, 60L),
    variables = list(
      state = c("S", "A", "U"),
      birth = c(30L, 10L, -10L),
      last_boosted_ica = c(39, 30, -1),
      icm = c(0.1, 0.2, 0.3),
      ica = c(1, 2, 3),
      zeta = c(1, 1.5, 2),
      zeta_group = c("1", "2", "3"),
      infectivity = c(0, 0.05, 0.0062),
      progression_rates = c(0, 1 / 200, 1 / 100),
      drug = c(0L, 0L, 0L),
      drug_time = c(-1L, -1L, -1L),
      last_pev_timestep = c(-1L, -1L, -1L),
      last_eff_pev_timestep = c(-1L, -1L, -1L),
      pev_profile = c(-1L, -1L, -1L),
      tbv_vaccinated = c(-1, -1, -1),
      net_time = c(1L, -1L, 1L),
      spray_time = c(5L, -1L, 7L),
      last_boosted_ib = c(39, 30, -1),
      last_boosted_iva = c(39, 30, -1),
      last_boosted_id = c(39, 30, -1),
      ivm = c(0.01, 0.02, 0.03),
      ib = c(1, 2, 3),
      iva = c(1, 2, 3),
      id = c(1, 2, 3)
    )
  )

  resample_stub <- stationary_human_initializer_resample
  mockery::stub(
    resample_stub,
    "sample.int",
    function(n, size, replace = FALSE, prob = NULL) c(3L, 1L)
  )
  sampled <- resample_stub(fake_library, size = 2L, parameters = params)

  expect_equal(sampled$state, c("U", "S"))
  expect_equal(sampled$birth, c(-70L, -10L))
  expect_equal(sampled$last_boosted_ica, c(-1, -1))
  expect_equal(sampled$last_boosted_ib, c(-1, -1))
  expect_equal(sampled$net_time, c(1L, 1L))
  expect_equal(sampled$spray_time, c(7L, 5L))
})

test_that('stationary human initializer pools snapshots from a resumed stochastic path', {
  clear_stationary_human_initializer_cache()

  params <- get_parameters(list(
    human_initialization = "stochastic_resample",
    human_initialization_burnin_timesteps = 30L,
    human_initialization_n_snapshots = 3L,
    human_initialization_snapshot_spacing = 10L,
    progress_bar = FALSE
  ))
  params$init_EIR <- 1

  call_log <- new.env(parent = emptyenv())
  call_log$calls <- list()

  library_stub <- get_stationary_human_initialization_library
  mockery::stub(
    library_stub,
    "run_resumable_simulation",
    function(timesteps,
             parameters = NULL,
             correlations = NULL,
             initial_state = NULL,
             restore_random_state = FALSE) {
      call_log$calls <- c(call_log$calls, list(list(
        timesteps = as.integer(timesteps),
        initial_state = initial_state,
        restore_random_state = restore_random_state
      )))
      list(state = list(timesteps = as.integer(timesteps)))
    }
  )
  mockery::stub(
    library_stub,
    "stationary_human_initializer_extract_library",
    function(state, parameters) {
      list(
        timesteps = as.integer(state$timesteps),
        reference_timestep = rep.int(as.integer(state$timesteps), 2L),
        variables = list(
          state = c("S", "A"),
          birth = c(as.integer(state$timesteps) - 5L, -10L)
        )
      )
    }
  )

  library <- library_stub(params)

  expect_equal(vapply(call_log$calls, `[[`, integer(1), "timesteps"), c(10L, 20L, 30L))
  expect_null(call_log$calls[[1]]$initial_state)
  expect_false(call_log$calls[[1]]$restore_random_state)
  expect_equal(call_log$calls[[2]]$initial_state$timesteps, 10L)
  expect_true(call_log$calls[[2]]$restore_random_state)
  expect_equal(call_log$calls[[3]]$initial_state$timesteps, 20L)
  expect_true(call_log$calls[[3]]$restore_random_state)
  expect_equal(library$reference_timestep, c(10L, 10L, 20L, 20L, 30L, 30L))
  expect_equal(library$variables$birth, c(5L, -10L, 15L, -10L, 25L, -10L))
})

test_that('stationary human initializer can select node-conditioned libraries', {
  params <- get_parameters(list(
    human_initialization = "stochastic_resample",
    human_initialization_node_index = 2L
  ))

  node_library <- new_stationary_human_initializer_node_conditioned_library(
    list(
      "1" = list(
        timesteps = 10L,
        reference_timestep = 10L,
        variables = list(state = "S", birth = -10L)
      ),
      "2" = list(
        timesteps = 20L,
        reference_timestep = 20L,
        variables = list(state = "A", birth = -25L)
      )
    )
  )

  params$human_initialization_library <- node_library
  selected <- get_stationary_human_initialization_library(params)

  expect_equal(selected$timesteps, 20L)
  expect_equal(selected$variables$state, "A")
  expect_equal(selected$variables$birth, -25L)
})

test_that('stationary human initializer carries node-specific vector regime metadata', {
  params <- get_parameters(list(
    human_initialization = "stochastic_resample",
    human_initialization_node_index = 2L
  ))

  node_library <- new_stationary_human_initializer_node_conditioned_library(
    list(
      "1" = list(
        timesteps = 10L,
        reference_timestep = 10L,
        variables = list(state = "S", birth = -10L)
      ),
      "2" = list(
        timesteps = 20L,
        reference_timestep = 20L,
        variables = list(state = "A", birth = -25L)
      )
    ),
    metadata = list(
      stationary_total_M_by_node = c("1" = 101, "2" = 202),
      stationary_init_foim_by_node = c("1" = 0.01, "2" = 0.02)
    )
  )

  params$human_initialization_library <- node_library
  selected <- get_stationary_human_initialization_library(params)

  expect_equal(selected$metadata$stationary_total_M, 202)
  expect_equal(selected$metadata$stationary_init_foim, 0.02)
})

test_that('stationary checkpoint library can select snapshots by index or timestep', {
  make_checkpoint <- function(timesteps) {
    list(
      state = list(
        timesteps = as.integer(timesteps),
        individual = list(dummy = TRUE),
        malariasimulationGD = list(dummy = TRUE)
      ),
      metadata = list(snapshot_index = as.integer(timesteps / 10))
    )
  }

  library <- new_stationary_checkpoint_library(
    checkpoints = list(
      make_checkpoint(10L),
      make_checkpoint(20L),
      make_checkpoint(30L)
    ),
    metadata = list(label = "test")
  )

  selected_default <- stationary_checkpoint_library_select(library)
  selected_by_index <- stationary_checkpoint_library_select(library, snapshot_index = 2L)
  selected_by_timestep <- stationary_checkpoint_library_select(library, timesteps = 10L)

  expect_equal(selected_default$state$timesteps, 30L)
  expect_equal(selected_by_index$state$timesteps, 20L)
  expect_equal(selected_by_timestep$state$timesteps, 10L)
})

test_that('stationary initializer snapshot timing is computed generically', {
  expect_equal(
    stationary_initializer_snapshot_timesteps(
      burnin_timesteps = 30L,
      n_snapshots = 3L,
      snapshot_spacing = 10L
    ),
    c(10L, 20L, 30L)
  )
})

test_that('stationary human initializer extracts node-conditioned metapop human libraries', {
  params <- get_parameters()
  categories <- stationary_human_initializer_categorical_variables(params)
  human_names <- stationary_human_initializer_human_variable_names(params)

  make_saved_categorical <- function(values, allowed) {
    out <- setNames(vector("list", length(allowed)), allowed)
    for (lvl in allowed) {
      out[[lvl]] <- which(values == lvl)
    }
    out
  }

  make_saved_block <- function(state_values, birth_values) {
    block <- list(
      state = make_saved_categorical(state_values, categories$state),
      birth = as.integer(birth_values),
      last_boosted_ica = c(-1, -1),
      icm = c(0.1, 0.2),
      ica = c(1, 2),
      zeta = c(1, 1.5),
      zeta_group = make_saved_categorical(c("1", "2"), categories$zeta_group),
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
    block[human_names]
  }

  saved_variables <- c(
    make_saved_block(c("S", "A"), c(-10L, -25L)),
    make_saved_block(c("U", "S"), c(-30L, -40L))
  )

  state <- list(
    timesteps = 30L,
    individual = list(variables = saved_variables)
  )

  library <- stationary_human_initializer_extract_metapop_library(
    state = state,
    parameters = list(params, params)
  )

  expect_true(stationary_human_initializer_is_node_conditioned_library(library))
  expect_equal(names(library$node_libraries), c("1", "2"))
  expect_equal(library$node_libraries[[1]]$variables$state, c("S", "A"))
  expect_equal(library$node_libraries[[2]]$variables$state, c("U", "S"))
  expect_equal(library$node_libraries[[1]]$variables$birth, c(-10L, -25L))
  expect_equal(library$node_libraries[[2]]$variables$birth, c(-30L, -40L))
})

test_that('stationary human initializer reads legacy human libraries without shifting fields', {
  params <- get_parameters()
  categories <- stationary_human_initializer_categorical_variables(params)
  legacy_names <- stationary_human_initializer_human_variable_names(
    params,
    include_human_slot_contact_multiplier = FALSE
  )

  make_saved_categorical <- function(values, allowed) {
    out <- setNames(vector("list", length(allowed)), allowed)
    for (lvl in allowed) {
      out[[lvl]] <- which(values == lvl)
    }
    out
  }

  legacy_block <- list(
    state = make_saved_categorical(c("S", "A"), categories$state),
    birth = c(-10L, -25L),
    last_boosted_ica = c(-1, -1),
    icm = c(0.1, 0.2),
    ica = c(1, 2),
    zeta = c(1, 1.5),
    zeta_group = make_saved_categorical(c("1", "2"), categories$zeta_group),
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
    last_boosted_iva = c(-2, -2),
    last_boosted_id = c(-3, -3),
    ivm = c(0.01, 0.02),
    ib = c(11, 12),
    iva = c(21, 22),
    id = c(31, 32)
  )
  saved_variables <- c(
    legacy_block[legacy_names],
    list(species = make_saved_categorical(c("gamb", "gamb"), "gamb"))
  )

  library <- stationary_human_initializer_extract_library(
    state = list(timesteps = 30L, individual = list(variables = saved_variables)),
    parameters = params
  )

  expect_equal(library$variables$state, c("S", "A"))
  expect_equal(library$variables$last_boosted_ib, c(-1, -1))
  expect_equal(library$variables$ib, c(11, 12))
  expect_equal(library$variables$id, c(31, 32))
  expect_equal(library$variables$human_slot_contact_multiplier, c(1, 1))
})

test_that('stationary human initializer keeps per-snapshot vector state metadata', {
  params <- get_parameters()
  categories <- stationary_human_initializer_categorical_variables(params)
  human_names <- stationary_human_initializer_human_variable_names(params)

  make_saved_categorical <- function(values, allowed) {
    out <- setNames(vector("list", length(allowed)), allowed)
    for (lvl in allowed) {
      out[[lvl]] <- which(values == lvl)
    }
    out
  }

  make_saved_block <- function(state_values, birth_values) {
    block <- list(
      state = make_saved_categorical(state_values, categories$state),
      birth = as.integer(birth_values),
      last_boosted_ica = c(-1, -1),
      icm = c(0.1, 0.2),
      ica = c(1, 2),
      zeta = c(1, 1.5),
      zeta_group = make_saved_categorical(c("1", "2"), categories$zeta_group),
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
    block[human_names]
  }

  state <- list(
    timesteps = 30L,
    individual = list(variables = c(
      make_saved_block(c("S", "A"), c(-10L, -25L)),
      make_saved_block(c("U", "S"), c(-30L, -40L))
    )),
    malariasimulationGD = list(
      "rng",
      list(list(dummy = TRUE), list(dummy = TRUE)),
      list(list(NULL), list(NULL)),
      list(list(state = 11), list(state = 22)),
      list(list(data.frame(timestep = 1, value = 1)), list(data.frame(timestep = 1, value = 2))),
      list(list(data.frame(timestep = 1, value = 3)), list(data.frame(timestep = 1, value = 4))),
      list(data.frame(timestep = 1, value = 5), data.frame(timestep = 1, value = 6))
    )
  )

  library <- stationary_human_initializer_extract_metapop_library(
    state = state,
    parameters = list(params, params)
  )

  expect_equal(library$node_libraries[[1]]$metadata$vector_state$timesteps, 30L)
  expect_equal(library$node_libraries[[1]]$metadata$vector_state$correlations$dummy, TRUE)
  expect_equal(library$node_libraries[[2]]$metadata$vector_state$solvers$state, 22)
})

test_that('stationary human initializer extracts full metapop snapshots with shared runtime state', {
  params <- get_parameters()
  categories <- stationary_human_initializer_categorical_variables(params)
  human_names <- stationary_human_initializer_human_variable_names(params)

  make_saved_categorical <- function(values, allowed) {
    out <- setNames(vector("list", length(allowed)), allowed)
    for (lvl in allowed) {
      out[[lvl]] <- which(values == lvl)
    }
    out
  }

  make_saved_block <- function(state_values, birth_values) {
    block <- list(
      state = make_saved_categorical(state_values, categories$state),
      birth = as.integer(birth_values),
      last_boosted_ica = c(-1, -1),
      icm = c(0.1, 0.2),
      ica = c(1, 2),
      zeta = c(1, 1.5),
      zeta_group = make_saved_categorical(c("1", "2"), categories$zeta_group),
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
    block[human_names]
  }

  state <- list(
    timesteps = 30L,
    individual = list(variables = c(
      make_saved_block(c("S", "A"), c(-10L, -25L)),
      make_saved_block(c("U", "S"), c(-30L, -40L))
    )),
    malariasimulationGD = list(
      "rng",
      list(list(dummy = TRUE), list(dummy = TRUE)),
      list(list(NULL), list(NULL)),
      list(list(state = 11), list(state = 22)),
      list(list(data.frame(timestep = 1, value = 1)), list(data.frame(timestep = 1, value = 2))),
      list(list(data.frame(timestep = 1, value = 3)), list(data.frame(timestep = 1, value = 4))),
      list(data.frame(timestep = 1, value = 5), data.frame(timestep = 1, value = 6))
    )
  )

  snapshot <- stationary_human_initializer_extract_metapop_snapshot(
    state = state,
    parameters = list(params, params)
  )

  expect_equal(snapshot$node_libraries[[1]]$variables$state, c("S", "A"))
  expect_equal(snapshot$node_libraries[[2]]$variables$state, c("U", "S"))
  expect_equal(snapshot$metadata$vector_state$timesteps, 30L)
  expect_equal(snapshot$metadata$vector_state$correlations[[1]]$dummy, TRUE)
  expect_equal(snapshot$metadata$vector_state$solvers[[2]]$state, 22)
})

test_that('stationary human initializer rebases retained runtime state to fresh-start day 0', {
  vector_state <- list(
    timesteps = 30L,
    correlations = list(mvnorm = matrix(0, nrow = 2, ncol = 1)),
    vector_models = list(NULL),
    solvers = list(list(
      t = 30L,
      state = 11:12,
      pending_inputs = list(timestep = 30L),
      last_completed_timestep = 30L
    )),
    lagged_eir = list(data.frame(timestep = 18:30, value = seq_len(13L))),
    lagged_transmission_eir = list(data.frame(timestep = 18:30, value = seq_len(13L) + 100)),
    lagged_infectivity = data.frame(timestep = 17:30, value = seq_len(14L) + 200)
  )

  rebased <- stationary_human_initializer_rebase_vector_state(vector_state)

  expect_equal(rebased$timesteps, 0L)
  expect_equal(rebased$correlations$mvnorm, matrix(0, nrow = 2, ncol = 1))
  expect_equal(rebased$solvers[[1]]$t, 30L)
  expect_equal(rebased$solvers[[1]]$pending_inputs$timestep, 30L)
  expect_equal(rebased$solvers[[1]]$last_completed_timestep, 30L)
  expect_equal(rebased$lagged_eir[[1]]$timestep, -12:0)
  expect_equal(rebased$lagged_transmission_eir[[1]]$timestep, -12:0)
  expect_equal(rebased$lagged_infectivity$timestep, -13:0)
})

test_that('create_variables can initialise humans from a stationary resample library', {
  clear_stationary_human_initializer_cache()

  params <- get_parameters(list(
    human_population = 3L,
    human_initialization = "stochastic_resample",
    human_initialization_burnin_timesteps = 1L,
    progress_bar = FALSE
  ))
  params$init_EIR <- 1
  params$init_foim <- 0

  fake_library <- list(
    timesteps = 20L,
    variables = list(
      state = c("S", "A", "U"),
      birth = c(10L, -5L, -100L),
      last_boosted_ica = c(-1, 19, 5),
      icm = c(0.1, 0.2, 0.3),
      ica = c(1, 2, 3),
      zeta = c(1, 1.5, 2),
      zeta_group = c("1", "2", "3"),
      infectivity = c(0, 0.05, 0.0062),
      progression_rates = c(0, 1 / 200, 1 / 100),
      drug = c(0L, 0L, 0L),
      drug_time = c(-1L, -1L, -1L),
      last_pev_timestep = c(-1L, -1L, -1L),
      last_eff_pev_timestep = c(-1L, -1L, -1L),
      pev_profile = c(-1L, -1L, -1L),
      tbv_vaccinated = c(-1, -1, -1),
      net_time = c(1L, -1L, 1L),
      spray_time = c(-1L, -1L, -1L),
      last_boosted_ib = c(-1, 19, 5),
      last_boosted_iva = c(-1, 19, 5),
      last_boosted_id = c(-1, 19, 5),
      ivm = c(0.01, 0.02, 0.03),
      ib = c(1, 2, 3),
      iva = c(1, 2, 3),
      id = c(1, 2, 3)
    )
  )

  create_variables_stub <- create_variables
  mockery::stub(
    create_variables_stub,
    "stationary_human_initializer_resample",
    function(library, size, parameters) {
      list(
        state = c("S", "A", "U"),
        birth = c(-10L, -25L, -120L),
        last_boosted_ica = c(-1, -1, -15),
        icm = c(0.1, 0.2, 0.3),
        ica = c(1, 2, 3),
        zeta = c(1, 1.5, 2),
        zeta_group = c("1", "2", "3"),
        infectivity = c(0, 0.05, 0.0062),
        progression_rates = c(0, 1 / 200, 1 / 100),
        drug = c(0L, 0L, 0L),
        drug_time = c(-1L, -1L, -1L),
        last_pev_timestep = c(-1L, -1L, -1L),
        last_eff_pev_timestep = c(-1L, -1L, -1L),
        pev_profile = c(-1L, -1L, -1L),
        tbv_vaccinated = c(-1, -1, -1),
        net_time = c(-19L, -1L, -19L),
        spray_time = c(-1L, -1L, -1L),
        last_boosted_ib = c(-1, -1, -15),
        last_boosted_iva = c(-1, -1, -15),
        last_boosted_id = c(-1, -1, -15),
        ivm = c(0.01, 0.02, 0.03),
        ib = c(1, 2, 3),
        iva = c(1, 2, 3),
        id = c(1, 2, 3)
      )
    }
  )
  mockery::stub(
    create_variables_stub,
    "get_stationary_human_initialization_library",
    function(parameters) fake_library
  )

  variables <- create_variables_stub(params)

  expect_equal(variables$state$get_values(), c("S", "A", "U"))
  expect_equal(variables$birth$get_values(), c(-10L, -25L, -120L))
  expect_equal(variables$last_boosted_ica$get_values(), c(-1, -1, -15))
  expect_equal(variables$ib$get_values(), c(1, 2, 3))
  expect_equal(variables$progression_rates$get_values(), c(0, 1 / 200, 1 / 100))
  expect_equal(variables$net_time$get_values(), c(-19L, -1L, -19L))
})

test_that('stochastic resample initialization preserves sampled infectivity values', {
  clear_stationary_human_initializer_cache()

  params <- get_parameters(list(
    human_population = 2L,
    human_initialization = "stochastic_resample",
    human_initialization_burnin_timesteps = 1L,
    progress_bar = FALSE
  ))
  params$init_EIR <- 1
  params$init_foim <- 50

  fake_library <- list(
    timesteps = 20L,
    variables = list(
      state = c("A", "U"),
      birth = c(-25L, -120L),
      last_boosted_ica = c(-1, -15),
      icm = c(0.2, 0.3),
      ica = c(2, 3),
      zeta = c(1.5, 2),
      zeta_group = c("2", "3"),
      infectivity = c(0.05, 0.0062),
      progression_rates = c(1 / 200, 1 / 100),
      drug = c(0L, 0L),
      drug_time = c(-1L, -1L),
      last_pev_timestep = c(-1L, -1L),
      last_eff_pev_timestep = c(-1L, -1L),
      pev_profile = c(-1L, -1L),
      tbv_vaccinated = c(-1, -1),
      net_time = c(-1L, -1L),
      spray_time = c(-1L, -1L),
      last_boosted_ib = c(-1, -15),
      last_boosted_iva = c(-1, -15),
      last_boosted_id = c(-1, -15),
      ivm = c(0.02, 0.03),
      ib = c(2, 3),
      iva = c(2, 3),
      id = c(2, 3)
    )
  )

  create_variables_stub <- create_variables
  mockery::stub(
    create_variables_stub,
    "get_stationary_human_initialization_library",
    function(parameters) fake_library
  )

  variables <- create_variables_stub(params)

  expect_equal(sort(variables$infectivity$get_values()), c(0.0062, 0.05))
})

test_that('stochastic resample initialization records stationary vector context', {
  clear_stationary_human_initializer_cache()

  params <- get_parameters(list(
    human_population = 2L,
    human_initialization = "stochastic_resample",
    human_initialization_burnin_timesteps = 1L,
    progress_bar = FALSE
  ))
  params$init_EIR <- 1

  fake_library <- list(
    timesteps = 20L,
    metadata = list(stationary_total_M = 123.4),
    variables = list(
      state = c("A", "U"),
      birth = c(-25L, -120L),
      last_boosted_ica = c(-1, -15),
      icm = c(0.2, 0.3),
      ica = c(2, 3),
      zeta = c(1.5, 2),
      zeta_group = c("2", "3"),
      infectivity = c(0.05, 0.0062),
      progression_rates = c(1 / 200, 1 / 100),
      drug = c(0L, 0L),
      drug_time = c(-1L, -1L),
      last_pev_timestep = c(-1L, -1L),
      last_eff_pev_timestep = c(-1L, -1L),
      pev_profile = c(-1L, -1L),
      tbv_vaccinated = c(-1, -1),
      net_time = c(-1L, -1L),
      spray_time = c(-1L, -1L),
      last_boosted_ib = c(-1, -15),
      last_boosted_iva = c(-1, -15),
      last_boosted_id = c(-1, -15),
      ivm = c(0.02, 0.03),
      ib = c(2, 3),
      iva = c(2, 3),
      id = c(2, 3)
    )
  )

  create_variables_stub <- create_variables
  mockery::stub(
    create_variables_stub,
    "get_stationary_human_initialization_library",
    function(parameters) fake_library
  )

  variables <- create_variables_stub(params)
  context <- attr(variables, "stationary_initialization_context")

  expect_equal(context$total_M, 123.4)
  expect_true(is.numeric(context$init_foim))
  expect_true(is.finite(context$init_foim))
  expect_gt(context$init_foim, 0)
})

test_that('stochastic resample picks a single retained snapshot vector state', {
  params <- get_parameters()

  snapshot_1 <- list(
    timesteps = 10L,
    reference_timestep = c(10L, 10L),
    variables = list(
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
    ),
    metadata = list(vector_state = list(timesteps = 10L))
  )
  snapshot_2 <- snapshot_1
  snapshot_2$timesteps <- 20L
  snapshot_2$reference_timestep <- c(20L, 20L)
  snapshot_2$metadata <- list(vector_state = list(timesteps = 20L))
  snapshot_2$variables$state <- c("U", "A")

  library <- list(
    metadata = list(snapshot_libraries = list(snapshot_1, snapshot_2)),
    variables = snapshot_1$variables
  )

  resample_stub <- stationary_human_initializer_resample
  call_count <- 0L
  mockery::stub(
    resample_stub,
    "sample.int",
    function(n, size, replace = FALSE, prob = NULL) {
      call_count <<- call_count + 1L
      if (call_count == 1L) {
        2L
      } else {
        c(1L, 2L)
      }
    }
  )
  sampled <- resample_stub(library, size = 2L, parameters = params)

  expect_equal(sampled$state, c("U", "A"))
  expect_equal(attr(sampled, "stationary_vector_state")$timesteps, 20L)
})

test_that('stationary human initializer applies sampled vector context to parameters', {
  params <- get_parameters()
  updated <- stationary_human_initializer_apply_context(
    parameters = params,
    context = list(total_M = 4321, init_foim = 0.0123, vector_control_time_offset = 20L)
  )

  expect_equal(updated$total_M, 4321)
  expect_equal(updated$init_foim, 0.0123)
  expect_equal(updated$vector_control_time_offset, 20L)
})

test_that('stationary human initializer can preserve baseline parameters during snapshot replay', {
  params <- get_parameters()
  baseline_total_M <- params$total_M
  baseline_init_foim <- params$init_foim

  context <- stationary_human_initializer_sample_context(
    sampled = list(
      birth = c(-25L, -120L),
      zeta = c(1.5, 2),
      net_time = c(-1L, -1L),
      spray_time = c(-1L, -1L),
      infectivity = c(0.05, 0.0062)
    ),
    parameters = params,
    library = list(metadata = list(
      preserve_baseline_parameters = TRUE,
      snapshot_timestep = 20L,
      stationary_total_M = 4321,
      stationary_init_foim = 0.0123
    ))
  )

  expect_true(context$preserve_baseline_parameters)
  expect_null(context$total_M)
  expect_null(context$init_foim)
  expect_equal(context$vector_control_time_offset, 20L)

  updated <- stationary_human_initializer_apply_context(params, context)
  expect_equal(updated$total_M, baseline_total_M)
  expect_equal(updated$init_foim, baseline_init_foim)
})

test_that('stochastic resample can preserve an exact retained human microstate', {
  params <- get_parameters()

  library <- list(
    timesteps = 20L,
    reference_timestep = c(20L, 20L),
    metadata = list(preserve_human_microstate = TRUE),
    variables = list(
      state = c("A", "U"),
      birth = c(-25L, -120L),
      last_boosted_ica = c(-1, 5),
      icm = c(0.2, 0.3),
      ica = c(2, 3),
      zeta = c(1.5, 2),
      zeta_group = c("2", "3"),
      infectivity = c(0.05, 0.0062),
      progression_rates = c(1 / 200, 1 / 100),
      drug = c(0L, 0L),
      drug_time = c(-1L, -1L),
      last_pev_timestep = c(-1L, -1L),
      last_eff_pev_timestep = c(-1L, -1L),
      pev_profile = c(-1L, -1L),
      tbv_vaccinated = c(-1, -1),
      net_time = c(-1L, -1L),
      spray_time = c(-1L, -1L),
      last_boosted_ib = c(-1, 5),
      last_boosted_iva = c(-1, 5),
      last_boosted_id = c(-1, 5),
      ivm = c(0.02, 0.03),
      ib = c(2, 3),
      iva = c(2, 3),
      id = c(2, 3)
    )
  )

  sampled <- stationary_human_initializer_resample(library, size = 2L, parameters = params)

  expect_equal(sampled$state, c("A", "U"))
  expect_equal(sampled$birth, c(-45L, -140L))
  expect_equal(sampled$last_boosted_ica, c(-1, -15))
  expect_equal(sampled$ib, c(2, 3))
})

test_that('create_variables can initialise humans from a node-conditioned stationary library', {
  clear_stationary_human_initializer_cache()

  params <- get_parameters(list(
    human_population = 2L,
    human_initialization = "stochastic_resample",
    human_initialization_node_index = 2L,
    progress_bar = FALSE
  ))
  params$init_EIR <- 1
  params$init_foim <- 0

  node_library <- new_stationary_human_initializer_node_conditioned_library(
    list(
      "1" = list(
        timesteps = 10L,
        reference_timestep = c(10L, 10L),
        variables = list(
          state = c("S", "S"),
          birth = c(-5L, -6L),
          last_boosted_ica = c(-1, -1),
          icm = c(0.1, 0.1),
          ica = c(1, 1),
          zeta = c(1, 1),
          zeta_group = c("1", "1"),
          infectivity = c(0, 0),
          progression_rates = c(0, 0),
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
          ivm = c(0.01, 0.01),
          ib = c(1, 1),
          iva = c(1, 1),
          id = c(1, 1)
        )
      ),
      "2" = list(
        timesteps = 20L,
        reference_timestep = c(20L, 20L),
        variables = list(
          state = c("A", "U"),
          birth = c(-25L, -120L),
          last_boosted_ica = c(-1, -15),
          icm = c(0.2, 0.3),
          ica = c(2, 3),
          zeta = c(1.5, 2),
          zeta_group = c("2", "3"),
          infectivity = c(0.05, 0.0062),
          progression_rates = c(1 / 200, 1 / 100),
          drug = c(0L, 0L),
          drug_time = c(-1L, -1L),
          last_pev_timestep = c(-1L, -1L),
          last_eff_pev_timestep = c(-1L, -1L),
          pev_profile = c(-1L, -1L),
          tbv_vaccinated = c(-1, -1),
          net_time = c(-1L, -1L),
          spray_time = c(-1L, -1L),
          last_boosted_ib = c(-1, -15),
          last_boosted_iva = c(-1, -15),
          last_boosted_id = c(-1, -15),
          ivm = c(0.02, 0.03),
          ib = c(2, 3),
          iva = c(2, 3),
          id = c(2, 3)
        )
      )
    )
  )

  params$human_initialization_library <- node_library

  variables <- create_variables(params)

  observed <- data.frame(
    state = variables$state$get_values(),
    birth = variables$birth$get_values(),
    last_boosted_ica = variables$last_boosted_ica$get_values(),
    ib = variables$ib$get_values(),
    progression_rates = variables$progression_rates$get_values()
  )
  observed <- observed[order(observed$birth), , drop = FALSE]

  expect_equal(observed$state, c("U", "A"))
  expect_equal(observed$birth, c(-140L, -45L))
  expect_equal(observed$last_boosted_ica, c(-35, -1))
  expect_equal(observed$ib, c(3, 2))
  expect_equal(observed$progression_rates, c(1 / 100, 1 / 200))
})

test_that('stationary human initializer prepares one shared metapop snapshot for all nodes', {
  make_node_library <- function(timesteps, states, births) {
    list(
      timesteps = as.integer(timesteps),
      reference_timestep = rep.int(as.integer(timesteps), length(states)),
      variables = list(
        state = states,
        birth = as.integer(births),
        last_boosted_ica = rep(-1, length(states)),
        icm = c(0.1, 0.2)[seq_along(states)],
        ica = seq_along(states),
        zeta = c(1, 1.5)[seq_along(states)],
        zeta_group = c("1", "2")[seq_along(states)],
        infectivity = c(0, 0.05)[seq_along(states)],
        progression_rates = c(0, 1 / 200)[seq_along(states)],
        drug = rep(0L, length(states)),
        drug_time = rep(-1L, length(states)),
        last_pev_timestep = rep(-1L, length(states)),
        last_eff_pev_timestep = rep(-1L, length(states)),
        pev_profile = rep(-1L, length(states)),
        tbv_vaccinated = rep(-1, length(states)),
        net_time = rep(-1L, length(states)),
        spray_time = rep(-1L, length(states)),
        last_boosted_ib = rep(-1, length(states)),
        last_boosted_iva = rep(-1, length(states)),
        last_boosted_id = rep(-1, length(states)),
        ivm = c(0.01, 0.02)[seq_along(states)],
        ib = seq_along(states),
        iva = seq_along(states),
        id = seq_along(states)
      )
    )
  }

  snapshot_1 <- list(
    node_libraries = list(
      "1" = make_node_library(10L, c("S", "A"), c(-10L, -25L)),
      "2" = make_node_library(10L, c("U", "S"), c(-30L, -40L))
    ),
    metadata = list(
      snapshot_index = 1L,
      snapshot_timestep = 10L,
      stationary_total_M_by_node = c("1" = 101, "2" = 202),
      stationary_init_foim_by_node = c("1" = 0.01, "2" = 0.02),
      vector_state = list(timesteps = 10L)
    )
  )
  snapshot_2 <- list(
    node_libraries = list(
      "1" = make_node_library(20L, c("A", "U"), c(-12L, -28L)),
      "2" = make_node_library(20L, c("D", "A"), c(-35L, -45L))
    ),
    metadata = list(
      snapshot_index = 2L,
      snapshot_timestep = 20L,
      stationary_total_M_by_node = c("1" = 303, "2" = 404),
      stationary_init_foim_by_node = c("1" = 0.03, "2" = 0.04),
      vector_state = list(timesteps = 20L)
    )
  )

  metapop_library <- stationary_human_initializer_set_accepted_metapop_snapshots(
    new_stationary_human_initializer_metapop_library(
      snapshots = list(snapshot_1, snapshot_2)
    ),
    accepted_snapshot_indices = 2L
  )
  params <- lapply(
    1:2,
    function(nd) {
      set_stationary_human_initialization_library(
        get_parameters(list(human_initialization = "stochastic_resample")),
        library = metapop_library,
        node_index = nd
      )
    }
  )

  prepared <- stationary_human_initializer_prepare_metapop_parameters(params)

  expect_equal(prepared$context$vector_state$timesteps, 20L)
  expect_equal(prepared$context$selected_snapshot_index, 2L)
  expect_equal(
    prepared$parameters[[1]]$human_initialization_library$variables$state,
    c("A", "U")
  )
  expect_equal(
    prepared$parameters[[2]]$human_initialization_library$variables$state,
    c("D", "A")
  )
  expect_true(prepared$parameters[[1]]$human_initialization_library$metadata$preserve_human_microstate)
  expect_true(prepared$parameters[[1]]$human_initialization_library$metadata$preserve_baseline_parameters)
  expect_equal(prepared$parameters[[1]]$human_initialization_library$metadata$stationary_total_M, 303)
  expect_equal(prepared$parameters[[2]]$human_initialization_library$metadata$stationary_init_foim, 0.04)
  expect_null(prepared$parameters[[1]]$human_initialization_node_index)
  expect_null(prepared$parameters[[2]]$human_initialization_node_index)
})

test_that('metapop stationary libraries can restrict runtime sampling to accepted snapshots', {
  make_node_library <- function(timesteps, states) {
    list(
      timesteps = as.integer(timesteps),
      reference_timestep = rep.int(as.integer(timesteps), length(states)),
      variables = list(
        state = states,
        birth = c(-10L, -20L),
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
    )
  }

  snapshots <- list(
    list(
      node_libraries = list("1" = make_node_library(10L, c("S", "A"))),
      metadata = list(snapshot_timestep = 10L)
    ),
    list(
      node_libraries = list("1" = make_node_library(20L, c("U", "A"))),
      metadata = list(snapshot_timestep = 20L)
    ),
    list(
      node_libraries = list("1" = make_node_library(30L, c("D", "U"))),
      metadata = list(snapshot_timestep = 30L)
    )
  )

  library <- stationary_human_initializer_set_accepted_metapop_snapshots(
    new_stationary_human_initializer_metapop_library(snapshots = snapshots),
    accepted_snapshot_indices = c(2L, 3L)
  )

  select_stub <- stationary_human_initializer_select_metapop_snapshot
  mockery::stub(
    select_stub,
    "sample.int",
    function(n, size, replace = FALSE, prob = NULL) 1L
  )

  selected <- select_stub(library)
  expect_equal(selected$metadata$selected_snapshot_index, 2L)
  expect_equal(selected$node_libraries[[1]]$variables$state, c("U", "A"))
})

test_that('metapop stationary libraries fail closed when no accepted snapshots remain', {
  make_node_library <- function(timesteps, states) {
    list(
      timesteps = as.integer(timesteps),
      reference_timestep = rep.int(as.integer(timesteps), length(states)),
      variables = list(
        state = states,
        birth = c(-10L, -20L),
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
    )
  }

  library <- stationary_human_initializer_set_accepted_metapop_snapshots(
    new_stationary_human_initializer_metapop_library(
      snapshots = list(
        list(
          node_libraries = list("1" = make_node_library(10L, c("S", "A"))),
          metadata = list(snapshot_timestep = 10L)
        ),
        list(
          node_libraries = list("1" = make_node_library(20L, c("U", "A"))),
          metadata = list(snapshot_timestep = 20L)
        )
      )
    ),
    accepted_snapshot_indices = integer(0)
  )

  expect_error(
    stationary_human_initializer_select_metapop_snapshot(library),
    "has no accepted snapshots"
  )
})
