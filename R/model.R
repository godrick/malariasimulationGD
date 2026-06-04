#' @title Run the simulation
#'
#' @description
#' Run the simulation for some time given some parameters. This currently
#' returns a dataframe with the number of individuals in each state at each
#' timestep.
#'
#' The resulting dataframe contains the following columns:
#'
#'  * timestep: the timestep for the row
#'  * infectivity: the infectivity from humans towards mosquitoes
#'  * FOIM: the force of infection towards mosquitoes (per species)
#'  * mu: the death rate of adult mosquitoes (per species)
#'  * EIR: the Entomological Inoculation Rate (per timestep, per species, over
#'  the whole population)
#'  * infectivity_weighted_I_* and vector_infectivity_mean_*: optional
#'  genotype-weighted human transmission signal and mean relative mosquito-to-
#'  human transmission multiplier among infectious females when
#'  `vector_infectivity_g` or `cube$b` is enabled
#'  * n_bitten: number of humans bitten by an infectious mosquito
#'  * n_treated: number of humans treated for clinical or severe malaria this timestep
#'  * n_infections: number of humans who get an asymptomatic, clinical or severe malaria this timestep
#'  * natural_deaths: number of humans who die from aging
#'  * humans_present: number of humans physically present in a node after the
#' explicit overnight mobility update
#'  * visitors_present: number of present humans whose home node differs from
#' the rendered node after the explicit overnight mobility update
#'  * residents_away: number of residents physically away from their home node
#' after the explicit overnight mobility update
#'  * trips_started: number of residents who started an off-home trip this
#' timestep under explicit overnight mobility
#'
#' When explicit human mobility is enabled in `run_metapop_simulation()`, these
#' four mobility columns are the only default mobility outputs. If
#' `human_mobility_store_diagnostics = TRUE`, the returned metapopulation output
#' list also carries optional attributes `OD_started_trips`,
#' `OD_active_overnight_stays`, and `mean_remaining_trip_duration`.
#'
#'  * S_count: number of humans who are Susceptible
#'  * A_count: number of humans who are Asymptomatic
#'  * D_count: number of humans who have the clinical malaria
#'  * U_count: number of subpatent infections in humans
#'  * Tr_count: number of detectable infections being treated in humans
#'  * ica_mean: the mean acquired immunity to clinical infection over the population of humans
#'  * icm_mean: the mean maternal immunity to clinical infection over the population of humans
#'  * ib_mean: the mean blood immunity to all infection over the population of humans
#'  * id_mean: the mean immunity from detection through microscopy over the population of humans
#'  * n: number of humans between an inclusive age range at this timestep. This
#' defaults to n_730_3650. Other age ranges can be set with
#' prevalence_rendering_min_ages and prevalence_rendering_max_ages parameters.
#'  * n_detect_lm (or pcr): number of humans with an infection detectable by microscopy (or pcr) between an inclusive age range at this timestep. This
#' defaults to n_detect_730_3650. Other age ranges can be set with
#' prevalence_rendering_min_ages and prevalence_rendering_max_ages parameters.
#'  * p_detect_lm (or pcr): the sum of probabilities of detection by microscopy (or pcr) between an
#' inclusive age range at this timestep. This
#' defaults to p_detect_730_3650. Other age ranges can be set with
#' prevalence_rendering_min_ages and prevalence_rendering_max_ages parameters.
#'  * n_inc: number of new infections for humans between an inclusive age range at this timestep.
#' incidence columns can be set with
#' incidence_rendering_min_ages and incidence_rendering_max_ages parameters.
#'  * p_inc: sum of probabilities of infection for humans between an inclusive age range at this timestep.
#' incidence columns can be set with
#' incidence_rendering_min_ages and incidence_rendering_max_ages parameters.
#'  * n_inc_clinical: number of new clinical infections for humans between an inclusive age range at this timestep.
#' clinical incidence columns can be set with
#' clinical_incidence_rendering_min_ages and clinical_incidence_rendering_max_ages parameters.
#'  * p_inc_clinical: sub of probabilities of clinical infection for humans between an inclusive age range at this timestep.
#' clinical incidence columns can be set with
#' clinical_incidence_rendering_min_ages and clinical_incidence_rendering_max_ages parameters.
#'  * n_inc_severe: number of new severe infections for humans between an inclusive age range at this timestep.
#' severe incidence columns can be set with
#' severe_incidence_rendering_min_ages and severe_incidence_rendering_max_ages parameters.
#'  * p_inc_severe: the sum of probabilities of severe infection for humans between an inclusive age range at this timestep.
#' severe incidence columns can be set with
#' severe_incidence_rendering_min_ages and severe_incidence_rendering_max_ages parameters.
#'  * E_count: number of mosquitoes in the early larval stage (per species)
#'  * L_count: number of mosquitoes in the late larval stage (per species)
#'  * P_count: number of mosquitoes in the pupal stage (per species)
#'  * Sm_count: number of adult female mosquitoes who are Susceptible (per
#'  species)
#'  * Pm_count: number of adult female mosquitoes who are incubating (per
#'  species)
#'  * Im_count: number of adult female mosquitoes who are infectious (per
#'  species)
#'  * rate_D_A: rate that humans transition from clinical disease to
#' asymptomatic
#'  * rate_A_U: rate that humans transition from asymptomatic to
#' subpatent
#'  * rate_U_S: rate that humans transition from subpatent to
#' susceptible
#'  * net_usage: the number people protected by a bed net
#'  * mosquito_deaths: number of adult female mosquitoes who die this timestep
#'  * n_drug_efficacy_failures: number of clinically treated individuals whose treatment failed due to drug efficacy
#'  * n_early_treatment_failure: number of clinically treated individuals who experienced early treatment failure
#'  * n_successfully_treated: number of clinically treated individuals who are treated successfully (includes individuals who experience slow parasite clearance)
#'  * n_slow_parasite_clearance: number of clinically treated individuals who experienced slow parasite clearance
#'
#' @param timesteps the number of timesteps to run the simulation for (in days)
#' @param parameters a named list of parameters to use
#' @param correlations correlation parameters
#' @return dataframe of results
#' @export
run_simulation <- function(
    timesteps,
    parameters = NULL,
    correlations = NULL
) {
  run_resumable_simulation(timesteps, parameters, correlations)$data
}

#' @title Run the simulation in a resumable way
#'
#' @description this function accepts an initial simulation state as an argument, and returns the
#' final state after running all of its timesteps. This allows one run to be resumed, possibly
#' having changed some of the parameters.
#' @param timesteps the timestep at which to stop the simulation
#' @param parameters a named list of parameters to use
#' @param correlations correlation parameters
#' @param initial_state the state from which the simulation is resumed
#' @param restore_random_state if TRUE, restore the random number generator's state from the checkpoint.
#' @return a list with `data` (dataframe of results) and `state` (final simulation state).
#' When `parameters$cube` is provided in hybrid mosquito mode, an additional
#' `mosquito_genotypes` entry is returned containing per-timestep female and male
#' adult genotype counts and the egg viability fraction `V`, and the same arrays
#' are attached as attributes on `data`. When genotype-specific mosquito-to-
#' human transmission is enabled via `vector_infectivity_g` or `cube$b`, and
#' `debug_genotypes = TRUE`, `data` also carries a
#' `mosquito_infectious_genotype_counts` attribute (timestep x genotype).
#' @export
run_resumable_simulation <- function(
    timesteps,
    parameters = NULL,
    correlations = NULL,
    initial_state = NULL,
    restore_random_state = FALSE
) {
  random_seed(ceiling(runif(1) * .Machine$integer.max))
  if (is.null(parameters)) {
    parameters <- get_parameters()
  }
  if (is.null(parameters$native_mosquito_backend)) {
    parameters$native_mosquito_backend <- TRUE
  }
  parameters <- validate_vector_infectivity_g_parameters(parameters)
  if (!is.null(parameters$releases) && is.null(parameters$releases_schedule)) {
    parameters <- set_releases(parameters, parameters$releases)
  }
  if (is.null(correlations)) {
    correlations <- get_correlation_parameters(parameters)
  }
  variables <- create_variables(parameters)
  stationary_initialization_context <- attr(variables, "stationary_initialization_context")
  if (is.null(initial_state)) {
    parameters <- stationary_human_initializer_apply_context(
      parameters,
      stationary_initialization_context
    )
  }
  if ((parameters$individual_mosquitoes || native_mosquito_backend_enabled(parameters)) &&
      !is.null(parameters$cube)) {
    cube_info <- cube_genotype_info(parameters$cube)
    parameters$mosquito_genotype_history <- new.env(parent = emptyenv())
    parameters$mosquito_genotype_history$female <- matrix(
      0,
      nrow = timesteps,
      ncol = cube_info$G,
      dimnames = list(NULL, cube_info$genotypesID)
    )
    parameters$mosquito_genotype_history$male <- matrix(
      0,
      nrow = timesteps,
      ncol = cube_info$G,
      dimnames = list(NULL, cube_info$genotypesID)
    )
    parameters$mosquito_genotype_history$V <- matrix(
      1,
      nrow = timesteps,
      ncol = length(parameters$species),
      dimnames = list(NULL, parameters$species)
    )
    parameters$mosquito_genotype_history$total_adults <- numeric(timesteps)
    parameters$mosquito_aquatic_genotype_history <- new.env(parent = emptyenv())
    parameters$mosquito_aquatic_genotype_history$E <- matrix(
      0,
      nrow = timesteps,
      ncol = cube_info$G,
      dimnames = list(NULL, cube_info$genotypesID)
    )
    parameters$mosquito_aquatic_genotype_history$L <- matrix(
      0,
      nrow = timesteps,
      ncol = cube_info$G,
      dimnames = list(NULL, cube_info$genotypesID)
    )
    parameters$mosquito_aquatic_genotype_history$P <- matrix(
      0,
      nrow = timesteps,
      ncol = cube_info$G,
      dimnames = list(NULL, cube_info$genotypesID)
    )
    if (isTRUE(parameters$debug_genotypes) &&
        (!is.null(parameters$vector_infectivity_g_by_species) ||
         (!is.null(parameters$cube) && !is.null(parameters$cube$b)))) {
      parameters$mosquito_infectious_genotype_history <- matrix(
        0,
        nrow = timesteps,
        ncol = cube_info$G,
        dimnames = list(NULL, cube_info$genotypesID)
      )
    }
  }
  events <- create_events(parameters)
  initialise_events(events, variables, parameters)
  renderer <- individual::Render$new(timesteps)
  populate_incidence_rendering_columns(renderer, parameters)
  attach_event_listeners(
    events,
    variables,
    parameters,
    correlations,
    renderer
  )
  vector_models <- parameterise_mosquito_models(parameters, timesteps)
  solvers <- parameterise_solvers(vector_models, parameters)

  lagged_eir <- create_lagged_eir(variables, solvers, parameters)
  lagged_transmission_eir <- create_lagged_transmission_eir(variables, solvers, parameters)
  lagged_infectivity <- create_lagged_infectivity(variables, parameters)
  if (is.null(initial_state)) {
    stationary_human_initializer_restore_vector_state(
      stationary_initialization_context,
      correlations = correlations,
      vector_models = vector_models,
      solvers = solvers,
      lagged_eir = lagged_eir,
      lagged_transmission_eir = lagged_transmission_eir,
      lagged_infectivity = lagged_infectivity
    )
  }

  stateful_objects <- list(
    random_state = RandomState$new(restore_random_state),
    correlations = correlations,
    vector_models = vector_models,
    solvers = solvers,
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_transmission_eir,
    lagged_infectivity = lagged_infectivity)

  if (!is.null(initial_state)) {
    individual::restore_object_state(
      initial_state$timesteps,
      stateful_objects,
      initial_state$malariasimulationGD)
  }

  individual_state <- individual::simulation_loop(
    processes = create_processes(
      renderer = renderer,
      variables = variables,
      events = events,
      parameters = parameters,
      models = vector_models,
      solvers = solvers,
      correlations = correlations,
      lagged_eir = lagged_eir,
      lagged_infectivity = lagged_infectivity,
      timesteps = timesteps,
      lagged_transmission_eir = lagged_transmission_eir
    ),
    variables = variables,
    events = events,
    timesteps = timesteps,
    state = initial_state$individual,
    restore_random_state = restore_random_state
  )

  final_state <- list(
    timesteps = timesteps,
    individual = individual_state,
    malariasimulationGD = individual::save_object_state(stateful_objects)
  )

  data <- renderer$to_dataframe()
  genotype_outputs <- NULL
  aquatic_genotype_outputs <- NULL
  infectious_genotype_counts_output <- parameters$mosquito_infectious_genotype_history
  release_schedule_output <- parameters$releases_schedule
  if (!is.null(parameters$mosquito_genotype_history)) {
    genotype_outputs <- list(
      female = parameters$mosquito_genotype_history$female,
      male = parameters$mosquito_genotype_history$male,
      V = parameters$mosquito_genotype_history$V,
      total_adults = parameters$mosquito_genotype_history$total_adults
    )
  }
  if (!is.null(parameters$mosquito_aquatic_genotype_history)) {
    aquatic_genotype_outputs <- list(
      E = parameters$mosquito_aquatic_genotype_history$E,
      L = parameters$mosquito_aquatic_genotype_history$L,
      P = parameters$mosquito_aquatic_genotype_history$P
    )
  }
  if (!is.null(initial_state)) {
    # Drop the timesteps we didn't simulate from the data.
    # It would just be full of NA.
    data <- data[-(1:initial_state$timesteps),]
    if (!is.null(release_schedule_output) && nrow(release_schedule_output) > 0) {
      release_schedule_output <- release_schedule_output[
        release_schedule_output$timestep > initial_state$timesteps,
        ,
        drop = FALSE
      ]
    }
    if (!is.null(genotype_outputs)) {
      keep <- -(1:initial_state$timesteps)
      genotype_outputs$female <- genotype_outputs$female[keep, , drop = FALSE]
      genotype_outputs$male <- genotype_outputs$male[keep, , drop = FALSE]
      genotype_outputs$V <- genotype_outputs$V[keep, , drop = FALSE]
      genotype_outputs$total_adults <- genotype_outputs$total_adults[keep]
    }
    if (!is.null(aquatic_genotype_outputs)) {
      keep <- -(1:initial_state$timesteps)
      aquatic_genotype_outputs$E <- aquatic_genotype_outputs$E[keep, , drop = FALSE]
      aquatic_genotype_outputs$L <- aquatic_genotype_outputs$L[keep, , drop = FALSE]
      aquatic_genotype_outputs$P <- aquatic_genotype_outputs$P[keep, , drop = FALSE]
    }
    if (!is.null(infectious_genotype_counts_output)) {
      infectious_genotype_counts_output <- infectious_genotype_counts_output[
        -(1:initial_state$timesteps),
        ,
        drop = FALSE
      ]
    }
  }
  if (!is.null(genotype_outputs)) {
    attr(data, "mosquito_genotype_counts_female") <- genotype_outputs$female
    attr(data, "mosquito_genotype_counts_male") <- genotype_outputs$male
    attr(data, "mosquito_genotype_V") <- genotype_outputs$V
    attr(data, "mosquito_genotype_total_adults") <- genotype_outputs$total_adults
  }
  if (!is.null(aquatic_genotype_outputs)) {
    attr(data, "mosquito_aquatic_genotype_E") <- aquatic_genotype_outputs$E
    attr(data, "mosquito_aquatic_genotype_L") <- aquatic_genotype_outputs$L
    attr(data, "mosquito_aquatic_genotype_P") <- aquatic_genotype_outputs$P
  }
  if (!is.null(infectious_genotype_counts_output)) {
    attr(data, "mosquito_infectious_genotype_counts") <- infectious_genotype_counts_output
  }
  if (!is.null(release_schedule_output)) {
    if (nrow(release_schedule_output) > 0 && "timestep" %in% names(data)) {
      for (sp in unique(release_schedule_output$species)) {
        col <- paste0("n_released_", sp)
        data[[col]] <- 0L
        sp_rows <- release_schedule_output$species == sp
        counts_by_t <- tapply(
          release_schedule_output$count[sp_rows],
          release_schedule_output$timestep[sp_rows],
          sum
        )
        if (length(counts_by_t) > 0) {
          t_idx <- match(as.integer(names(counts_by_t)), data$timestep)
          valid <- !is.na(t_idx)
          data[[col]][t_idx[valid]] <- as.integer(counts_by_t[valid])
        }
      }
    }
    attr(data, "mosquito_release_schedule") <- release_schedule_output
  }

  res <- list(data = data, state = final_state)
  if (!is.null(genotype_outputs)) {
    res$mosquito_genotypes <- genotype_outputs
  }
  if (!is.null(aquatic_genotype_outputs)) {
    res$mosquito_aquatic_genotypes <- aquatic_genotype_outputs
  }
  res
}

#' @title Run a metapopulation model
#'
#' @param timesteps the number of timesteps to run the simulation for (in days)
#' @param parameters a list of model parameter lists for each population
#' @param correlations a list of correlation parameters for each population
#' (default: NULL)
#' @param mixing_tt a vector of time steps for each mixing matrix
#' @param export_mixing a list of matrices of coefficients for exportation of infectivity.
#' Rows = origin sites, columns = destinations. Each matrix element
#' describes the mixing pattern from destination to origin. Each matrix element must
#' be between 0 and 1. Each matrix is activated at the corresponding timestep in mixing_tt
#' @param import_mixing a list of matrices of coefficients for importation of
#' infectivity.
#' @param p_captured_tt a vector of time steps for each p_captured matrix
#' @param p_captured a list of matrices representing the probability that
#' travel between sites is intervened by a test and treat border check.
#' Dimensions are the same as for `export_mixing`
#' @param p_success the probability that an individual who has tested positive
#' (through an RDT) successfully clears their infection through treatment
#' @param initial_state optional resumable metapopulation state to continue from
#' @param restore_random_state if TRUE, restore the random number generator
#'   state from `initial_state`
#' @param return_state if TRUE, also return the final resumable state
#' @param render_output if FALSE, skip output rendering/data-frame materialization
#'   and run in a state-only mode suitable for checkpoint burn-in
#' @param return_summary if TRUE, return a lightweight final-state summary
#' @return a list of dataframe of model outputs as in `run_simulation`. When
#' explicit human mobility diagnostics are enabled, the output list also has
#' attributes `OD_started_trips`, `OD_active_overnight_stays`, and
#' `mean_remaining_trip_duration`.
#' @export
run_metapop_simulation <- function(
  timesteps,
  parameters,
  correlations = NULL,
  mixing_tt,
  export_mixing,
  import_mixing,
  p_captured_tt,
  p_captured,
  p_success,
  initial_state = NULL,
  restore_random_state = FALSE,
  return_state = FALSE,
  render_output = TRUE,
  return_summary = FALSE
  ) {
  random_seed(ceiling(runif(1) * .Machine$integer.max))

  for (mixing in list(export_mixing, import_mixing)) {
    if (!is.list(mixing)) {
      stop('mixing arguments must be a list of mixing matrices')
    }

    if (length(mixing_tt) != length(mixing)) {
      stop('mixing_tt must be the same length as mixing matrices')
    }

    for (i in seq_along(mixing)) {
      if (nrow(mixing[[i]]) != ncol(mixing[[i]])) {
        stop(sprintf('mixing matrix %d must be square', i))
      }
      if (nrow(mixing[[i]]) != length(parameters)) {
        stop(sprintf("mixing matrix %d's rows must match length of parameters", i))
      }
      if (!all(vlapply(seq_along(parameters), function(x) approx_sum(mixing[[i]][x,], 1)))) {
        warning(sprintf("all of mixing matrix %d's rows must sum to 1", i))
      }
      if (!all(vlapply(seq_along(parameters), function(x) approx_sum(mixing[[i]][,x], 1)))) {
        warning(sprintf('mixing matrix %d is asymmetrical', i))
      }
    }
    if (length(mixing_tt) != length(mixing)) {
      stop('mixing_tt must be the same size as mixing')
    }
  }

  for (i in seq_along(p_captured)) {
    if (nrow(p_captured[[i]]) != ncol(p_captured[[i]])) {
      stop(sprintf('p_captured matrix %d must be square', i))
    }
    if (!all(diag(p_captured[[i]]) == 0)) {
      warning(sprintf('p_captured matrix %d has a non-zero diagonal', i))
    }
  }

  if (!is.numeric(mixing_tt)) {
    stop('mixing_tt must be numeric')
  }

  if (length(p_captured_tt) != length(p_captured)) {
    stop('p_captured_tt must be the same length as p_captured')
  }

  if (is.null(correlations)) {
    correlations <- lapply(parameters, get_correlation_parameters)
  }
  parameters <- lapply(parameters, function(p) {
    if (is.null(p$native_mosquito_backend)) {
      p$native_mosquito_backend <- TRUE
    }
    p
  })
  validate_explicit_human_mobility_metapop(
    parameters = parameters,
    export_mixing = export_mixing,
    import_mixing = import_mixing,
    p_captured = p_captured,
    p_success = p_success
  )
  human_mobility_move_probs <- human_mobility_resolve_move_probs(parameters)
  if (!is.null(human_mobility_move_probs)) {
    parameters <- lapply(seq_along(parameters), function(i) {
      p <- parameters[[i]]
      p$human_mobility_enabled <- TRUE
      p$human_move_probs <- human_mobility_move_probs
      p$human_mobility_node_index <- as.integer(i)
      p$human_mobility_n_nodes <- as.integer(length(parameters))
      p
    })
  }
  metapop_stationary_initialization_context <- NULL
  if (is.null(initial_state)) {
    metapop_initialization <- stationary_human_initializer_prepare_metapop_parameters(parameters)
    parameters <- metapop_initialization$parameters
    metapop_stationary_initialization_context <- metapop_initialization$context
  }
  parameters <- lapply(parameters, validate_vector_infectivity_g_parameters)
  parameters <- lapply(parameters, function(p) {
    p$mosquito_genotype_history <- NULL
    p$mosquito_aquatic_genotype_history <- NULL
    p$mosquito_infectious_genotype_history <- NULL
    if (!is.null(p$releases) && is.null(p$releases_schedule)) {
      p <- set_releases(p, p$releases)
    }
    if (isTRUE(render_output) && !is.null(p$cube)) {
      cube_info <- cube_genotype_info(p$cube)
      p$mosquito_genotype_history <- new.env(parent = emptyenv())
      p$mosquito_genotype_history$female <- matrix(
        0, nrow = timesteps, ncol = cube_info$G,
        dimnames = list(NULL, cube_info$genotypesID)
      )
      p$mosquito_genotype_history$male <- matrix(
        0, nrow = timesteps, ncol = cube_info$G,
        dimnames = list(NULL, cube_info$genotypesID)
      )
      p$mosquito_genotype_history$V <- matrix(
        1, nrow = timesteps, ncol = length(p$species),
        dimnames = list(NULL, p$species)
      )
      p$mosquito_genotype_history$total_adults <- numeric(timesteps)
      p$mosquito_aquatic_genotype_history <- new.env(parent = emptyenv())
      p$mosquito_aquatic_genotype_history$E <- matrix(
        0, nrow = timesteps, ncol = cube_info$G,
        dimnames = list(NULL, cube_info$genotypesID)
      )
      p$mosquito_aquatic_genotype_history$L <- matrix(
        0, nrow = timesteps, ncol = cube_info$G,
        dimnames = list(NULL, cube_info$genotypesID)
      )
      p$mosquito_aquatic_genotype_history$P <- matrix(
        0, nrow = timesteps, ncol = cube_info$G,
        dimnames = list(NULL, cube_info$genotypesID)
      )
    }
    p
  })
  variables <- lapply(parameters, create_variables)
  stationary_initialization_context <- lapply(
    variables,
    function(v) attr(v, "stationary_initialization_context")
  )
  if (is.null(initial_state)) {
    parameters <- Map(
      function(p, v) {
        stationary_human_initializer_apply_context(
          p,
          attr(v, "stationary_initialization_context")
        )
      },
      parameters,
      stationary_initialization_context
    )
  }
  events <- lapply(parameters, create_events)
  renderer <- lapply(
    parameters,
    function(.) {
      if (isTRUE(render_output)) {
        individual::Render$new(timesteps)
      } else {
        NullRender$new(timesteps)
      }
    }
  )
  if (isTRUE(render_output)) {
    populate_metapopulation_incidence_rendering_columns(renderer, parameters)
  }
  human_mobility_context <- create_human_mobility_context(
    parameters = parameters,
    variables = variables,
    timesteps = timesteps,
    render_output = render_output
  )
  for (i in seq_along(parameters)) {
    # NOTE: forceAndCall is necessary here to make sure i refers to the current
    # iteration
    forceAndCall(
      3,
      initialise_events,
      events[[i]],
      variables[[i]],
      parameters[[i]]
    )
    forceAndCall(
      5,
      attach_event_listeners,
      events[[i]],
      variables[[i]],
      parameters[[i]],
      correlations[[i]],
      renderer[[i]]
    )
  }
  native_backend <- native_metapop_backend_enabled(parameters)
  if (native_backend) {
    backend <- parameterise_native_metapop_backends(parameters, timesteps)
    vector_models <- backend$models
    solvers <- backend$solvers
  } else {
    vector_models <- lapply(parameters, parameterise_mosquito_models, timesteps = timesteps)
    solvers <- lapply(
      seq_along(parameters),
      function(i) parameterise_solvers(vector_models[[i]], parameters[[i]])
    )
  }
  lagged_eir <- lapply(
    seq_along(parameters),
    function(i) create_lagged_eir(variables[[i]], solvers[[i]], parameters[[i]])
  )
  lagged_transmission_eir <- lapply(
    seq_along(parameters),
    function(i) create_lagged_transmission_eir(variables[[i]], solvers[[i]], parameters[[i]])
  )
  lagged_infectivity <- lapply(
    seq_along(parameters),
    function(i) create_lagged_infectivity(variables[[i]], parameters[[i]])
  )
  human_exposure_lag_context <- create_human_exposure_lag_context(
    parameters = parameters,
    variables = variables,
    solvers = solvers,
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_transmission_eir
  )
  human_infectivity_lag_context <- create_human_infectivity_lag_context(
    parameters = parameters,
    variables = variables
  )
  if (is.null(initial_state)) {
    if (!is.null(metapop_stationary_initialization_context$vector_state)) {
      stationary_human_initializer_restore_vector_state(
        metapop_stationary_initialization_context,
        correlations = correlations,
        vector_models = vector_models,
        solvers = solvers,
        lagged_eir = lagged_eir,
        lagged_transmission_eir = lagged_transmission_eir,
        lagged_infectivity = lagged_infectivity
      )
    } else {
      for (i in seq_along(parameters)) {
        stationary_human_initializer_restore_vector_state(
          stationary_initialization_context[[i]],
          correlations = correlations[[i]],
          vector_models = vector_models[[i]],
          solvers = solvers[[i]],
          lagged_eir = lagged_eir[[i]],
          lagged_transmission_eir = lagged_transmission_eir[[i]],
          lagged_infectivity = lagged_infectivity[[i]]
        )
      }
    }
  }

  stateful_objects <- list(
    random_state = RandomState$new(restore_random_state),
    correlations = correlations,
    vector_models = vector_models,
    solvers = solvers,
    lagged_eir = lagged_eir,
    lagged_transmission_eir = lagged_transmission_eir,
    lagged_infectivity = lagged_infectivity
  )

  if (!is.null(initial_state)) {
    individual::restore_object_state(
      initial_state$timesteps,
      stateful_objects,
      initial_state$malariasimulationGD
    )
    human_exposure_lag_restore_state(
      human_exposure_lag_context,
      initial_state$timesteps,
      initial_state$malariasimulationGD$human_exposure_lag
    )
    human_infectivity_lag_restore_state(
      human_infectivity_lag_context,
      initial_state$timesteps,
      initial_state$malariasimulationGD$human_infectivity_lag
    )
  }

  mixing_fn <- time_cached(
    create_transmission_mixer(
      variables = variables,
      parameters = parameters,
      lagged_eir = lagged_eir,
      lagged_infectivity = lagged_infectivity,
      mixing_tt = mixing_tt,
      export_mixing = export_mixing,
      import_mixing = import_mixing,
      p_captured_tt = p_captured_tt,
      p_captured = p_captured,
      p_success = p_success,
      lagged_transmission_eir = lagged_transmission_eir
    )
  )
    
  processes <- lapply(
    seq_along(parameters),
    function(i) {
      create_processes(
        renderer = renderer[[i]],
        variables = variables[[i]],
        events = events[[i]],
        parameters = parameters[[i]],
        models = vector_models[[i]],
        solvers = solvers[[i]],
        correlations = correlations[[i]],
        lagged_eir = lagged_eir[[i]],
        lagged_infectivity = lagged_infectivity[[i]],
        timesteps = timesteps,
        mixing_fn = mixing_fn,
        mixing_index = i,
        lagged_transmission_eir = lagged_transmission_eir[[i]],
        enable_rendering = render_output,
        human_mobility_context = human_mobility_context,
        human_exposure_lag_context = human_exposure_lag_context,
        human_infectivity_lag_context = human_infectivity_lag_context
      )
    }
  )

  individual_state <- individual::simulation_loop(
    processes = if (native_backend) native_interleave_metapop_processes(processes) else unlist(processes),
    variables = unlist(variables),
    events = unlist(events),
    timesteps = timesteps,
    state = if (!is.null(initial_state)) initial_state$individual else NULL,
    restore_random_state = restore_random_state
  )

  malariasimulationGD_state <- individual::save_object_state(stateful_objects)
  human_exposure_lag_state <- human_exposure_lag_save_state(human_exposure_lag_context)
  if (!is.null(human_exposure_lag_state)) {
    malariasimulationGD_state$human_exposure_lag <- human_exposure_lag_state
  }
  human_infectivity_lag_state <- human_infectivity_lag_save_state(human_infectivity_lag_context)
  if (!is.null(human_infectivity_lag_state)) {
    malariasimulationGD_state$human_infectivity_lag <- human_infectivity_lag_state
  }

  final_state <- list(
    timesteps = timesteps,
    individual = individual_state,
    malariasimulationGD = malariasimulationGD_state
  )
  final_summary <- NULL
  if (isTRUE(return_summary)) {
    final_summary <- list(
      total_M_by_node = final_total_M_by_node(
        solvers = solvers,
        parameters = parameters,
        variables = variables
      )
    )
  }

  trim_resumed_output <- function(data, initial_timesteps) {
    if (is.null(initial_timesteps) || initial_timesteps <= 0L) {
      return(data)
    }

    initial_timesteps <- as.integer(initial_timesteps)
    keep <- seq.int(initial_timesteps + 1L, nrow(data))
    trimmed <- data[keep, , drop = FALSE]
    if ("timestep" %in% names(trimmed)) {
      trimmed$timestep <- seq_len(nrow(trimmed))
    }

    attr_names <- names(attributes(data))
    preserve_attrs <- setdiff(attr_names, c("names", "row.names", "class"))
    for (attr_name in preserve_attrs) {
      attr_val <- attr(data, attr_name)
      if (is.matrix(attr_val) && nrow(attr_val) == nrow(data)) {
        attr(trimmed, attr_name) <- attr_val[keep, , drop = FALSE]
      } else if (is.numeric(attr_val) && length(attr_val) == nrow(data)) {
        attr(trimmed, attr_name) <- attr_val[keep]
      } else {
        attr(trimmed, attr_name) <- attr_val
      }
    }

    release_schedule <- attr(trimmed, "mosquito_release_schedule")
    if (!is.null(release_schedule) && "timestep" %in% names(release_schedule)) {
      keep_release <- release_schedule$timestep > initial_timesteps
      release_schedule <- release_schedule[keep_release, , drop = FALSE]
      release_schedule$timestep <- release_schedule$timestep - initial_timesteps
      attr(trimmed, "mosquito_release_schedule") <- release_schedule
    }

    trimmed
  }

  outputs <- if (isTRUE(render_output)) {
    lapply(
      seq_along(renderer),
      function(i) {
        data <- renderer[[i]]$to_dataframe()
        if (!is.null(parameters[[i]]$mosquito_genotype_history)) {
          attr(data, "mosquito_genotype_counts_female") <- parameters[[i]]$mosquito_genotype_history$female
          attr(data, "mosquito_genotype_counts_male") <- parameters[[i]]$mosquito_genotype_history$male
          attr(data, "mosquito_genotype_V") <- parameters[[i]]$mosquito_genotype_history$V
          attr(data, "mosquito_genotype_total_adults") <- parameters[[i]]$mosquito_genotype_history$total_adults
        }
        if (!is.null(parameters[[i]]$mosquito_aquatic_genotype_history)) {
          attr(data, "mosquito_aquatic_genotype_E") <- parameters[[i]]$mosquito_aquatic_genotype_history$E
          attr(data, "mosquito_aquatic_genotype_L") <- parameters[[i]]$mosquito_aquatic_genotype_history$L
          attr(data, "mosquito_aquatic_genotype_P") <- parameters[[i]]$mosquito_aquatic_genotype_history$P
        }
        if (!is.null(parameters[[i]]$releases_schedule)) {
          attr(data, "mosquito_release_schedule") <- parameters[[i]]$releases_schedule
        }
        if (!is.null(initial_state)) {
          data <- trim_resumed_output(data, initial_state$timesteps)
        }
        data
      }
    )
  } else {
    NULL
  }
  outputs <- attach_human_mobility_diagnostics(
    outputs,
    human_mobility_context,
    initial_timesteps = if (!is.null(initial_state)) initial_state$timesteps else NULL
  )

  if (!isTRUE(return_state)) {
    if (isTRUE(return_summary)) {
      return(list(data = outputs, summary = final_summary))
    }
    return(outputs)
  }

  result <- list(
    data = outputs,
    state = final_state
  )
  if (isTRUE(return_summary)) {
    result$summary <- final_summary
  }
  result
}

#' @title Run the simulation with repetitions
#'
#' @param timesteps the number of timesteps to run the simulation for
#' @param repetitions n times to run the simulation
#' @param overrides a named list of parameters to use instead of defaults
#' @param parallel execute runs in parallel
#' @export
run_simulation_with_repetitions <- function(
    timesteps,
    repetitions,
    overrides = list(),
    parallel = FALSE
) {
  if (parallel) {
    fapply <- parallel::mclapply
  } else {
    fapply <- lapply
  }
  dfs <- fapply(
    seq(repetitions),
    function(repetition) {
      df <- run_simulation(timesteps, overrides)
      df$repetition <- repetition
      df
    }
  )
  do.call("rbind", dfs)
}
