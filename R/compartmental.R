ODE_INDICES <- c(E = 1, L = 2, P = 3)
ADULT_ODE_INDICES <- c(Sm = 4, Pm = 5, Im = 6)

parameterise_mosquito_models <- function(parameters, timesteps = 0) {
  if (native_mosquito_backend_enabled(parameters)) {
    native_warn_tau_leap_backend(parameters)
    return(
      lapply(
        seq_along(parameters$species),
        function(i) create_native_mosquito_model(parameters, i, timesteps)
      )
    )
  }

  lapply(
    seq_along(parameters$species),
    function(i) {
      p <- parameters$species_proportions[[i]]
      m <- p * parameters$total_M
      # Baseline carrying capacity
      k0 <- calculate_carrying_capacity(parameters, m, i)
      # Create the carrying capacity object
      k_timeseries <- create_timeseries(size = length(parameters$carrying_capacity_timesteps), k0)
      if(parameters$carrying_capacity){
        for(j in 1:length(parameters$carrying_capacity_timesteps)){
          timeseries_push(
            k_timeseries,
            parameters$carrying_capacity_scalers[j,i] * k0,
            parameters$carrying_capacity_timesteps[j]
          )
        }
      }
      growth_model <- create_aquatic_mosquito_model(
        species_beta_value(parameters, i),
        parameters$del,
        parameters$me,
        k_timeseries,
        parameters$gamma,
        parameters$dl,
        parameters$ml,
        parameters$dpl,
        parameters$mup,
        m,
        parameters$model_seasonality,
        parameters$g0,
        parameters$g,
        parameters$h,
        calculate_R_bar(parameters),
        parameters$mum[[i]],
        parameters$blood_meal_rates[[i]],
        parameters$rainfall_floor
      )

      if (parameters$individual_mosquitoes && !is.null(parameters$cube)) {
        cube_info <- cube_genotype_info(parameters$cube)
        egg_p <- rep.int(0, cube_info$G)
        egg_p[[cube_info$wild_type_index]] <- 1
        aquatic_mosquito_model_set_egg_proportions(growth_model, egg_p)
      }
      
      if (!parameters$individual_mosquitoes) {
        susceptible <- initial_mosquito_counts(
          parameters,
          i,
          parameters$init_foim,
          m
        )[ADULT_ODE_INDICES['Sm']]
        return(
          AdultMosquitoModel$new(create_adult_mosquito_model(
            growth_model,
            parameters$mum[[i]],
            parameters$dem,
            susceptible * parameters$init_foim,
            parameters$init_foim
          ))
        )
      }
      genotype_state <- NULL
      if (!is.null(parameters$cube)) {
        cube_info <- cube_genotype_info(parameters$cube)
        initial_female_adults <- floor(sum(
          initial_mosquito_counts(
            parameters,
            i,
            parameters$init_foim,
            m
          )[ADULT_ODE_INDICES]
        ))
        male_counts <- rep.int(0, cube_info$G)
        male_counts[[cube_info$wild_type_index]] <- initial_female_adults
        genotype_state <- list(
          male_counts = male_counts,
          last_V = 1
        )
      }
      AquaticMosquitoModel$new(
        growth_model,
        cube = parameters$cube,
        genotype_state = genotype_state
      )
    }
  )
}

parameterise_solvers <- function(models, parameters) {
  if (native_mosquito_backend_enabled(parameters)) {
    return(lapply(models, create_native_mosquito_solver, parameters = parameters))
  }

  lapply(
    seq_along(models),
    function(i) {
      m <- parameters$species_proportions[[i]] * parameters$total_M
      init <- initial_mosquito_counts(parameters, i, parameters$init_foim, m)
      if (!parameters$individual_mosquitoes) {
        return(
          Solver$new(create_adult_solver(
            models[[i]]$.model,
            init,
            parameters$r_tol,
            parameters$a_tol,
            parameters$ode_max_steps
          ))
        )
      }
      init_aquatic <- init[ODE_INDICES]
      if (!is.null(parameters$cube)) {
        cube_info <- cube_genotype_info(parameters$cube)
        if (cube_info$G > 1L) {
          init_aquatic <- rep.int(0, 3L * cube_info$G)
          e_idx <- aquatic_genotype_stage_indices(cube_info$G, "E")
          l_idx <- aquatic_genotype_stage_indices(cube_info$G, "L")
          p_idx <- aquatic_genotype_stage_indices(cube_info$G, "P")
          wt <- cube_info$wild_type_index
          init_aquatic[[e_idx[[wt]]]] <- init[[ODE_INDICES[["E"]]]]
          init_aquatic[[l_idx[[wt]]]] <- init[[ODE_INDICES[["L"]]]]
          init_aquatic[[p_idx[[wt]]]] <- init[[ODE_INDICES[["P"]]]]
        }
      }
      Solver$new(create_aquatic_solver(
        models[[i]]$.model,
        init_aquatic,
        parameters$r_tol,
        parameters$a_tol,
        parameters$ode_max_steps
      ))
    }
  )
}

create_compartmental_rendering_process <- function(renderer, solvers, parameters) {
  if (native_mosquito_backend_enabled(parameters)) {
    return(function(timestep) {
      aquatic_history <- NULL
      cube_info <- NULL
      aquatic_E <- aquatic_L <- aquatic_P <- NULL
      genotype_history <- NULL
      female_total <- male_total <- NULL
      V_by_species <- NULL
      if (!is.null(parameters$cube) && !is.null(parameters$mosquito_aquatic_genotype_history)) {
        cube_info <- cube_genotype_info(parameters$cube)
        aquatic_history <- parameters$mosquito_aquatic_genotype_history
        aquatic_E <- rep.int(0, cube_info$G)
        aquatic_L <- rep.int(0, cube_info$G)
        aquatic_P <- rep.int(0, cube_info$G)
      }
      if (!is.null(parameters$cube) && !is.null(parameters$mosquito_genotype_history)) {
        cube_info <- cube_genotype_info(parameters$cube)
        genotype_history <- parameters$mosquito_genotype_history
        female_total <- rep.int(0, cube_info$G)
        male_total <- rep.int(0, cube_info$G)
        V_by_species <- rep(1, length(parameters$species))
      }

      for (s_i in seq_along(solvers)) {
        summary <- solvers[[s_i]]$get_summary()
        for (stage_name in c("E", "L", "P", "Sm", "Pm", "Im")) {
          renderer$render(
            paste0(stage_name, "_", parameters$species[[s_i]], "_count"),
            summary$totals[[stage_name]],
            timestep
          )
        }
        if (!is.null(aquatic_history)) {
          aquatic_E <- aquatic_E + summary$aquatic$E
          aquatic_L <- aquatic_L + summary$aquatic$L
          aquatic_P <- aquatic_P + summary$aquatic$P
        }
        if (!is.null(genotype_history)) {
          female_total <- female_total + summary$female
          male_total <- male_total + summary$male
          V_by_species[[s_i]] <- calc_pg_V_from_cube(parameters$cube, summary$female, summary$male)$V
        }
      }

      if (!is.null(aquatic_history)) {
        aquatic_history$E[timestep, ] <- aquatic_E
        aquatic_history$L[timestep, ] <- aquatic_L
        aquatic_history$P[timestep, ] <- aquatic_P
      }
      if (!is.null(genotype_history)) {
        genotype_history$female[timestep, ] <- female_total
        genotype_history$male[timestep, ] <- male_total
        genotype_history$V[timestep, ] <- V_by_species
        genotype_history$total_adults[timestep] <- sum(female_total + male_total)
      }
    })
  }

  if (parameters$individual_mosquitoes) {
    indices <- ODE_INDICES
  } else {
    indices <- c(ODE_INDICES, ADULT_ODE_INDICES)
  }
  
  function(timestep) {
    aquatic_history <- NULL
    cube_info <- NULL
    aquatic_E <- aquatic_L <- aquatic_P <- NULL
    if (parameters$individual_mosquitoes && !is.null(parameters$cube) &&
        !is.null(parameters$mosquito_aquatic_genotype_history)) {
      cube_info <- cube_genotype_info(parameters$cube)
      aquatic_history <- parameters$mosquito_aquatic_genotype_history
      aquatic_E <- rep.int(0, cube_info$G)
      aquatic_L <- rep.int(0, cube_info$G)
      aquatic_P <- rep.int(0, cube_info$G)
    }
    for (s_i in seq_along(solvers)) {
      if (parameters$species_proportions[[s_i]] > 0) {
        row <- solvers[[s_i]]$get_states()
      } else {
        row <- rep(0, if (is.null(cube_info)) length(indices) else 3L * cube_info$G)
      }
      if (parameters$individual_mosquitoes) {
        stage_totals <- if (is.null(parameters$cube)) {
          setNames(row[ODE_INDICES], names(ODE_INDICES))
        } else {
          aquatic_stage_totals(row, parameters$cube)
        }
        for (stage_name in names(ODE_INDICES)) {
          renderer$render(
            paste0(stage_name, '_', parameters$species[[s_i]], '_count'),
            stage_totals[[stage_name]],
            timestep
          )
        }
        if (!is.null(aquatic_history)) {
          aquatic_E <- aquatic_E + aquatic_stage_values_by_genotype(row, parameters$cube, "E")
          aquatic_L <- aquatic_L + aquatic_stage_values_by_genotype(row, parameters$cube, "L")
          aquatic_P <- aquatic_P + aquatic_stage_values_by_genotype(row, parameters$cube, "P")
        }
      } else {
        for (i in seq_along(indices)) {
          renderer$render(
            paste0(names(indices)[[i]], '_', parameters$species[[s_i]], '_count'),
            row[[i]],
            timestep
          )
        }
      }
    }
    if (!is.null(aquatic_history)) {
      aquatic_history$E[timestep, ] <- aquatic_E
      aquatic_history$L[timestep, ] <- aquatic_L
      aquatic_history$P[timestep, ] <- aquatic_P
    }
  }
}

#' @title Step mosquito solver
#' @description calculates total_M per species and updates the vector ode
#'
#' @param solvers for each species
#' @noRd
create_solver_stepping_process <- function(solvers, parameters) {
  function(timestep) {
    for (i in seq_along(solvers)) {
      if (parameters$species_proportions[[i]] > 0) {
        solvers[[i]]$step()
      }
    }
  }
}

Solver <- R6::R6Class(
  'Solver',
  private = list(
    .solver = NULL
  ),
  public = list(
    initialize = function(solver) {
      private$.solver <- solver
    },
    step = function() {
      solver_step(private$.solver)
    },
    get_states = function() {
      solver_get_states(private$.solver)
    },

    # This is the same as `get_states`, just exposed under the interface that
    # is expected of stateful objects.
    save_state = function() {
      solver_get_states(private$.solver)
    },
    restore_state = function(t, state) {
      solver_set_states(private$.solver, t, state)
    }
  )
)

AquaticMosquitoModel <- R6::R6Class(
  'AquaticMosquitoModel',
  public = list(
    .model = NULL,
    cube = NULL,
    genotype_state = NULL,
    initialize = function(model, cube = NULL, genotype_state = NULL) {
      self$.model <- model
      self$cube <- cube
      self$genotype_state <- genotype_state
    },

    # The aquatic ODE state is stored in the solver. We only persist the
    # auxiliary genotype tracking state when present.
    save_state = function() {
      if (is.null(self$genotype_state)) {
        return(NULL)
      }
      list(genotype_state = self$genotype_state)
    },
    restore_state = function(t, state) {
      if (is.null(state)) {
        return(invisible(NULL))
      }
      if (is.list(state) && !is.null(state$genotype_state)) {
        self$genotype_state <- state$genotype_state
      }
    }
  )
)

AdultMosquitoModel <- R6::R6Class(
  'AdultMosquitoModel',
  public = list(
    .model = NULL,
    initialize = function(model) {
      self$.model <- model
    },
    save_state = function() {
      adult_mosquito_model_save_state(self$.model)
    },
    restore_state = function(t, state) {
      adult_mosquito_model_restore_state(self$.model, state)
    }
  )
)
