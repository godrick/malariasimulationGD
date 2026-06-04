#' @title Biting process
#' @description
#' This is the biting process. It results in human and mosquito infection and
#' mosquito death.
#' @param renderer the model renderer object
#' @param solvers mosquito ode solvers
#' @param models mosquito ode models
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param parameters model pararmeters
#' @param lagged_infectivity a list of LaggedValue objects with historical sums
#' of infectivity, one for every metapopulation
#' @param lagged_eir a LaggedValue class with historical EIRs
#' @param mixing_fn a function to retrieve the mixed EIR and infectivity based
#' on the other populations
#' @param mixing_index an index for this population's position in the
#' lagged_infectivity list (default: 1)
#' @param infection_outcome competing hazards object for infection rates
#' @param human_exposure_lag_context optional per-human exposure lag context
#' @param human_infectivity_lag_context optional per-human infectivity lag context
#' @param timestep the current timestep
#' @noRd
create_biting_process <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  parameters,
  lagged_infectivity,
  lagged_eir,
  mixing_fn = NULL,
  mixing_index = 1,
  infection_outcome,
  lagged_transmission_eir = lagged_eir,
  human_exposure_lag_context = NULL,
  human_infectivity_lag_context = NULL
  ) {
  function(timestep) {
    # Calculate combined EIR
    age <- get_age(variables$birth$get_values(), timestep)
    bitten <- simulate_bites(
      renderer,
      solvers,
      models,
      variables,
      events,
      age,
      parameters,
      timestep,
      lagged_infectivity,
      lagged_eir,
      mixing_fn,
      mixing_index,
      lagged_transmission_eir,
      human_exposure_lag_context,
      human_infectivity_lag_context
    )

    if (human_mobility_enabled(parameters)) {
      infection_input <- human_exposure_lag_get_infection_input(
        human_exposure_lag_context,
        mixing_index,
        timestep,
        parameters
      )
      simulate_infection(
        variables,
        events,
        bitten$bitten_humans,
        bitten$n_bites_per_person,
        age,
        parameters,
        timestep,
        renderer,
        infection_outcome,
        infection_input$transmission_multiplier,
        infection_exposure = infection_input$infection_exposure
      )
      return(invisible(NULL))
    }

    simulate_infection(
      variables,
      events,
      bitten$bitten_humans,
      bitten$n_bites_per_person,
      age,
      parameters,
      timestep,
      renderer,
      infection_outcome,
      bitten$transmission_multiplier
    )
  }
}

#' @importFrom stats rpois
simulate_bites <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  age,
  parameters,
  timestep,
  lagged_infectivity,
  lagged_eir,
  mixing_fn = NULL,
  mixing_index = 1,
  lagged_transmission_eir = lagged_eir,
  human_exposure_lag_context = NULL,
  human_infectivity_lag_context = NULL
  ) {
  bitten_humans <- individual::Bitset$new(parameters$human_population)
  n_bites_per_person <- numeric(0)
  vector_infectivity_active <- !is.null(parameters$vector_infectivity_g_by_species) ||
    (!is.null(parameters$cube) && !is.null(parameters$cube$b))
  infectivity_weighting_active <- FALSE
  native_backend <- native_mosquito_backend_enabled(parameters)
  
  human_infectivity <- human_infectivity_lag_current_values(
    timestep,
    variables,
    parameters
  )
  renderer$render('infectivity', mean(human_infectivity), timestep)
  
  # Calculate pi (the relative biting rate for each human)
  psi <- unique_biting_rate(age, parameters)
  zeta <- variables$zeta$get_values()
  human_slot_contact_multiplier <- human_slot_contact_multiplier_values(variables)
  slot_contact_rate_multiplier <- human_slot_contact_rate_multiplier(
    zeta,
    psi,
    human_slot_contact_multiplier
  )
  .pi <- human_pi(zeta, psi, human_slot_contact_multiplier)
  
  # Get some indices for later
  genotype_tracking <- !is.null(parameters$cube)
  if (native_backend) {
    if (vector_infectivity_active && !genotype_tracking) {
      stop("Internal error: genotype-specific vector_infectivity_g requires genotype tracking state")
    }
    if (genotype_tracking) {
      cube_info <- cube_genotype_info(parameters$cube)
      female_geno_totals <- rep.int(0, cube_info$G)
      male_geno_totals <- rep.int(0, cube_info$G)
      V_by_species <- rep(1, length(parameters$species))
      infectious_geno_totals <- rep.int(0, cube_info$G)
    }
  } else if (parameters$individual_mosquitoes) {
    infectious_index <- variables$mosquito_state$get_index_of('Im')
    susceptible_index <- variables$mosquito_state$get_index_of('Sm')
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
    omega_by_species <- vector("list", length(parameters$species))
    genotype_tracking <- !is.null(parameters$cube) && !is.null(variables$geno_id)
    if (vector_infectivity_active && !genotype_tracking) {
      stop("Internal error: genotype-specific vector_infectivity_g requires genotype tracking state")
    }
    if (genotype_tracking) {
      cube_info <- cube_genotype_info(parameters$cube)
      female_geno_totals <- rep.int(0, cube_info$G)
      male_geno_totals <- rep.int(0, cube_info$G)
      mu_by_species <- rep(0, length(parameters$species))
      female_totals_by_species <- rep(0L, length(parameters$species))
      V_by_species <- rep(1, length(parameters$species))
      infectious_geno_totals <- rep.int(0, cube_info$G)
    }
  }
  if (!native_backend && parameters$individual_mosquitoes && genotype_tracking && genotype_debug_enabled(parameters, timestep)) {
    for (s_i in seq_along(parameters$species)) {
      genotype_debug_log_counts(
        parameters,
        timestep,
        "BITE_START",
        parameters$species[[s_i]],
        genotype_debug_species_counts(variables, models, parameters, s_i),
        extra = "visible state at start of biting (after release process; emergence females may still be queued)"
      )
    }
  }
  
  EIR <- 0
  transmission_signal <- 0
  current_exposure_by_species <- numeric(length(parameters$species))
  current_weighted_exposure_by_species <- numeric(length(parameters$species))

  for (s_i in seq_along(parameters$species)) {
    species_name <- parameters$species[[s_i]]
    solver_states <- solvers[[s_i]]$get_states()
    p_bitten <- prob_bitten(timestep, variables, s_i, parameters)
    Q0 <- parameters$Q0[[s_i]]
    W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
    Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
    f <- blood_meal_rate(s_i, Z, parameters)
    a <- apply_runtime_contact_multiplier(
      .human_blood_meal_rate(f, s_i, W, parameters) *
        slot_contact_rate_multiplier,
      s_i,
      parameters
    )
    lambda <- effective_biting_rates(a, .pi, p_bitten)

    if (native_backend) {
      summary <- solvers[[s_i]]$get_summary()
      transmission_weights <- human_transmission_weights_for_species(
        parameters,
        species_name
      )
      infectivity_weighting_active <- infectivity_weighting_active || !is.null(transmission_weights)
      n_infectious <- sum(summary$infectious)
      n_transmission_infectious <- n_infectious
      if (!is.null(transmission_weights)) {
        n_transmission_infectious <- sum(summary$infectious * transmission_weights)
        infectious_geno_totals <- infectious_geno_totals + summary$infectious
        renderer$render(
          paste0("infectivity_weighted_I_", species_name),
          n_transmission_infectious,
          timestep
        )
        renderer$render(
          paste0("vector_infectivity_mean_", species_name),
          if (n_infectious > 0) n_transmission_infectious / n_infectious else NA_real_,
          timestep
        )
      }
      if (genotype_tracking) {
        female_geno_totals <- female_geno_totals + summary$female
        male_geno_totals <- male_geno_totals + summary$male
        V_by_species[[s_i]] <- calc_pg_V_from_cube(models[[s_i]]$cube, summary$female, summary$male)$V
      }
    } else if (parameters$individual_mosquitoes) {
      species_index <- variables$species$get_index_of(
        parameters$species[[s_i]]
      )$and(adult_index)
      transmission_weights <- human_transmission_weights_for_species(
        parameters,
        species_name
      )
      infectivity_weighting_active <- infectivity_weighting_active || !is.null(transmission_weights)
      if (is.null(transmission_weights)) {
        n_infectious <- calculate_infectious_individual(
          s_i,
          variables,
          infectious_index,
          adult_index,
          species_index,
          parameters
        )
        n_transmission_infectious <- n_infectious
      } else {
        infectious_counts_g <- calculate_infectious_individual_genotype_counts(
          variables,
          infectious_index,
          species_index,
          cube_info$G
        )
        n_infectious <- sum(infectious_counts_g)
        n_transmission_infectious <- sum(infectious_counts_g * transmission_weights)
        infectious_geno_totals <- infectious_geno_totals + infectious_counts_g
        renderer$render(
          paste0("infectivity_weighted_I_", species_name),
          n_transmission_infectious,
          timestep
        )
        renderer$render(
          paste0("vector_infectivity_mean_", species_name),
          if (n_infectious > 0) n_transmission_infectious / n_infectious else NA_real_,
          timestep
        )
      }
    } else {
      n_infectious <- calculate_infectious_compartmental(solver_states)
      n_transmission_infectious <- n_infectious
    }
    
    # Store both unweighted bite exposure and genotype-weighted human
    # transmission signal for later. The latter is used only to scale infection
    # probability after bites are generated, matching the Imperial-style
    # genotype-specific b0 semantics.
    species_exposure <- n_infectious * a
    species_weighted_exposure <- n_transmission_infectious * a
    lagged_eir[[s_i]]$save(
      species_exposure,
      timestep
    )
    lagged_transmission_eir[[s_i]]$save(
      species_weighted_exposure,
      timestep
    )
    current_exposure_by_species[[s_i]] <- species_exposure
    current_weighted_exposure_by_species[[s_i]] <- species_weighted_exposure

    # lagged EIR
    if (is.null(mixing_fn)) {
      species_eir <- lagged_eir[[s_i]]$get(timestep - parameters$de)
      species_transmission_eir <- lagged_transmission_eir[[s_i]]$get(timestep - parameters$de)
    } else {
      mixed_transmission <- mixing_fn(timestep=timestep)
      species_eir <- mixed_transmission$eir[mixing_index, s_i]
      if (is.null(mixed_transmission$transmission_eir)) {
        species_transmission_eir <- species_eir
      } else {
        species_transmission_eir <- mixed_transmission$transmission_eir[mixing_index, s_i]
      }
    }

    renderer$render(paste0('EIR_', species_name), species_eir, timestep)
    EIR <- EIR + species_eir
    transmission_signal <- transmission_signal + species_transmission_eir
    if(parameters$parasite == "falciparum"){
      # p.f model factors eir by psi
      expected_bites <- species_eir * mean(psi)
    } else if (parameters$parasite == "vivax"){
      # p.v model standardises biting rate het to eir
      expected_bites <- species_eir
    }

    if (!is.finite(expected_bites)) {
      stop(
        sprintf(
          paste(
            "Non-finite expected_bites in simulate_bites",
            "(timestep=%s, species=%s, species_eir=%s, mean_psi=%s)."
          ),
          as.character(timestep),
          species_name,
          format(species_eir, digits = 16),
          format(mean(psi), digits = 16)
        ),
        call. = FALSE
      )
    }

    if (expected_bites > 0) {
      n_bites <- rpois(1, expected_bites)
      if (n_bites > 0) {
        bitten <- fast_weighted_sample(n_bites, lambda)
        bitten_humans$insert(bitten)
        renderer$render('n_bitten', bitten_humans$size(), timestep)
        if(parameters$parasite == "vivax"){
          # p.v must pass through the number of bites per person
          n_bites_per_person <- tabulate(bitten, nbins = length(lambda))
        }
      }
    }

    lagged_infectivity$save(sum(human_infectivity * .pi), timestep)

    if (human_mobility_enabled(parameters)) {
      infectivity <- human_infectivity_lag_get_node_reservoir(
        human_infectivity_lag_context,
        mixing_index,
        timestep
      )
    } else if (is.null(mixing_fn)) {
      infectivity <- lagged_infectivity$get(timestep - parameters$delay_gam)
    } else {
      infectivity <- mixing_fn(timestep=timestep)$inf[[mixing_index]]
    }

    foim <- calculate_foim(a, infectivity)
    renderer$render(paste0('FOIM_', species_name), foim, timestep)
    mu <- death_rate(f, W, Z, s_i, parameters)
    renderer$render(paste0('mu_', species_name), mu, timestep)
    
    if (native_backend) {
      native_mosquito_model_update(
        models[[s_i]],
        timestep,
        mu,
        foim,
        f
      )
    } else if (parameters$individual_mosquitoes) {
      # update the ODE with stats for ovoposition calculations
      effective_total_M <- species_index$size()
      if (genotype_tracking && !is.null(models[[s_i]]$genotype_state)) {
        female_counts <- tabulate(
          variables$geno_id$get_values(species_index),
          nbins = cube_info$G
        )
        male_counts <- models[[s_i]]$genotype_state$male_counts
        pgv <- calc_pg_V_from_cube(models[[s_i]]$cube, female_counts, male_counts)
        V_by_species[[s_i]] <- pgv$V
        effective_total_M <- effective_total_M * pgv$V
        aquatic_mosquito_model_set_egg_proportions(models[[s_i]]$.model, as.numeric(pgv$p))
        female_geno_totals <- female_geno_totals + female_counts
        male_geno_totals <- male_geno_totals + male_counts
        mu_by_species[[s_i]] <- mu
        omega_g <- cube_omega_vector(models[[s_i]]$cube, cube_info$G, cube_info$genotypesID)
        omega_by_species[[s_i]] <- if (all(omega_g == 1)) NULL else omega_g
        female_totals_by_species[[s_i]] <- species_index$size()
        models[[s_i]]$genotype_state$last_V <- pgv$V
      }

      aquatic_mosquito_model_update(
        models[[s_i]]$.model,
        effective_total_M,
        f,
        mu
      )
      
      # update the individual mosquitoes
      susceptible_species_index <- susceptible_index$copy()$and(species_index)
      omega_for_species <- NULL
      if (exists("omega_by_species", inherits = FALSE) &&
          length(omega_by_species) >= s_i) {
        omega_for_species <- omega_by_species[[s_i]]
      }
      
      biting_effects_individual(
        variables,
        foim,
        events,
        s_i,
        susceptible_species_index,
        species_index,
        mu,
        parameters,
        timestep,
        omega_by_genotype = omega_for_species
      )
    } else {
      adult_mosquito_model_update(
        models[[s_i]]$.model,
        mu,
        foim,
        solver_states[[ADULT_ODE_INDICES['Sm']]],
        f
      )
    }
  }

  if ((native_backend || parameters$individual_mosquitoes) && infectivity_weighting_active &&
      !is.null(parameters$mosquito_infectious_genotype_history)) {
    parameters$mosquito_infectious_genotype_history[timestep, ] <- infectious_geno_totals
  }

  if (!native_backend && parameters$individual_mosquitoes && genotype_tracking) {
    history <- parameters$mosquito_genotype_history
    if (!is.null(history)) {
      history$female[timestep, ] <- female_geno_totals
      history$male[timestep, ] <- male_geno_totals
      history$V[timestep, ] <- V_by_species
      history$total_adults[timestep] <- sum(female_geno_totals + male_geno_totals)
      if (genotype_debug_enabled(parameters, timestep)) {
        for (s_i in seq_along(parameters$species)) {
          genotype_debug_log(
            parameters,
            timestep,
            "HWRITE",
            parameters$species[[s_i]],
            sprintf(
              "row=%d writes visible counts before male death update: F{%s} M{%s}",
              timestep,
              genotype_debug_fmt_counts(female_geno_totals, names(parameters$mosquito_genotype_history$female[timestep, ])),
              genotype_debug_fmt_counts(male_geno_totals, names(parameters$mosquito_genotype_history$male[timestep, ]))
            )
          )
        }
      }
    }
    for (s_i in seq_along(parameters$species)) {
      if (is.null(models[[s_i]]$genotype_state)) {
        next
      }
      male_counts <- models[[s_i]]$genotype_state$male_counts
      omega_by_genotype <- NULL
      if (exists("omega_by_species", inherits = FALSE) &&
          length(omega_by_species) >= s_i) {
        omega_by_genotype <- omega_by_species[[s_i]]
      }
      if (length(male_counts) == 1L &&
          (is.null(omega_by_genotype) || omega_by_genotype[[1]] == 1)) {
        # Keep the trivial cube case RNG-free and aligned with the implicit 1:1 sex ratio.
        male_counts[[1]] <- female_totals_by_species[[s_i]]
      } else if (is.null(omega_by_genotype)) {
        death_prob <- min(mu_by_species[[s_i]], 1)
        male_counts <- stats::rbinom(length(male_counts), size = male_counts, prob = 1 - death_prob)
      } else {
        death_prob <- pmin(mu_by_species[[s_i]] * as.numeric(omega_by_genotype), 1)
        male_counts <- stats::rbinom(length(male_counts), size = male_counts, prob = 1 - death_prob)
      }
      models[[s_i]]$genotype_state$male_counts <- male_counts
      if (genotype_debug_enabled(parameters, timestep)) {
        genotype_debug_log_counts(
          parameters,
          timestep,
          "AFTER_MALE_DEATH",
          parameters$species[[s_i]],
          genotype_debug_species_counts(variables, models, parameters, s_i),
          extra = "female deaths are scheduled via mosquito_death event (delay=0)"
        )
      }
    }
  }

  human_exposure_lag_record_node(
    human_exposure_lag_context,
    mixing_index,
    timestep,
    current_exposure_by_species,
    current_weighted_exposure_by_species
  )

  transmission_multiplier <- if (EIR > 0) transmission_signal / EIR else 1

  list(
    bitten_humans = bitten_humans,
    n_bites_per_person = n_bites_per_person,
    transmission_multiplier = transmission_multiplier
  )
}


# =================
# Utility functions
# =================

calculate_eir <- function(species, solvers, variables, parameters, timestep) {
  a <- human_blood_meal_rate(species, variables, parameters, timestep)
  infectious <- calculate_infectious(species, solvers, variables, parameters)
  infectious * a
}

calculate_transmission_eir <- function(species, solvers, variables, parameters, timestep) {
  a <- human_blood_meal_rate(species, variables, parameters, timestep)
  infectious <- calculate_transmission_infectious(species, solvers, variables, parameters)
  infectious * a
}

effective_biting_rates <- function(a, .pi, p_bitten) {
  a * .pi * p_bitten$prob_bitten / sum(.pi * p_bitten$prob_bitten_survives)
}

legacy_contact_multiplier_for_human_slots <- function(parameters) {
  contact_multiplier <- parameters$contact_multiplier
  if (is.null(contact_multiplier)) {
    return(1)
  }

  value <- NULL
  if (is.list(contact_multiplier)) {
    if (length(contact_multiplier) == 1L) {
      value <- contact_multiplier[[1L]]
    }
  } else if (is.numeric(contact_multiplier)) {
    if (length(contact_multiplier) == 1L) {
      value <- contact_multiplier[[1L]]
    }
  }

  value <- as.numeric(value)
  if (length(value) != 1L || !is.finite(value) || value <= 0) {
    stop(
      "`parameters$contact_multiplier` must be NULL or a single positive finite number.",
      call. = FALSE
    )
  }

  value
}

# Scalar legacy contact multipliers are folded into human_slot_contact_multiplier
# at human-variable initialization. Runtime scaling is retained as a no-op
# compatibility wrapper for older internal call sites.
runtime_contact_multiplier <- function(parameters, species) {
  species_names <- as.character(parameters$species)
  species <- as.integer(species)
  if (length(species) != 1L || is.na(species) ||
      species < 1L || species > length(species_names)) {
    stop("species must index one configured mosquito species.", call. = FALSE)
  }
  legacy_contact_multiplier_for_human_slots(parameters)
  1
}

apply_runtime_contact_multiplier <- function(a, species, parameters) {
  a * runtime_contact_multiplier(parameters, species)
}

calculate_infectious <- function(species, solvers, variables, parameters) {
  if (native_mosquito_backend_enabled(parameters)) {
    summary <- solvers[[species]]$get_summary()
    return(sum(summary$infectious))
  }
  if (parameters$individual_mosquitoes) {
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
    species_name <- parameters$species[[species]]
    species_index <- variables$species$get_index_of(
      species_name
    )$and(adult_index)
    return(
      calculate_infectious_individual(
        species,
        variables,
        variables$mosquito_state$get_index_of('Im'),
        adult_index,
        species_index,
        parameters
      )
    )
  }
  calculate_infectious_compartmental(solvers[[species]]$get_states())
}

calculate_transmission_infectious <- function(species, solvers, variables, parameters) {
  species_name <- parameters$species[[species]]
  transmission_weights <- human_transmission_weights_for_species(parameters, species_name)
  if (is.null(transmission_weights)) {
    return(calculate_infectious(species, solvers, variables, parameters))
  }

  if (native_mosquito_backend_enabled(parameters)) {
    summary <- solvers[[species]]$get_summary()
    return(sum(summary$infectious * transmission_weights))
  }
  if (parameters$individual_mosquitoes) {
    if (is.null(parameters$cube) || is.null(variables$geno_id)) {
      stop("Internal error: genotype-specific vector_infectivity_g requires genotype tracking state")
    }
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
    species_index <- variables$species$get_index_of(
      species_name
    )$and(adult_index)
    cube_info <- cube_genotype_info(parameters$cube)
    counts_g <- calculate_infectious_individual_genotype_counts(
      variables,
      variables$mosquito_state$get_index_of('Im'),
      species_index,
      cube_info$G
    )
    return(sum(counts_g * transmission_weights))
  }
  calculate_infectious(species, solvers, variables, parameters)
}

calculate_infectious_individual_genotype_counts <- function(
  variables,
  infectious_index,
  species_index,
  G
  ) {
  infectious_species_index <- infectious_index$copy()$and(species_index)
  if (infectious_species_index$size() == 0) {
    return(rep.int(0, G))
  }
  tabulate(variables$geno_id$get_values(infectious_species_index), nbins = G)
}

calculate_infectious_individual <- function(
  species,
  variables,
  infectious_index,
  adult_index,
  species_index,
  parameters
  ) {
  infectious_index$copy()$and(species_index)$size()
}

calculate_infectious_compartmental <- function(solver_states) {
  max(solver_states[[ADULT_ODE_INDICES['Im']]], 0)
}

intervention_coefficient <- function(p_bitten) {
  p_bitten$prob_bitten / sum(p_bitten$prob_bitten_survives)
}

validate_human_slot_contact_multiplier <- function(values, size, label = "human_slot_contact_multiplier") {
  if (!is.numeric(values)) {
    stop(sprintf("`%s` must be numeric.", label), call. = FALSE)
  }
  if (length(values) == 1L) {
    values <- rep(as.numeric(values), size)
  }
  if (length(values) != size) {
    stop(
      sprintf(
        "`%s` must have length 1 or match the human population size (%d).",
        label,
        size
      ),
      call. = FALSE
    )
  }
  values <- as.numeric(values)
  if (any(!is.finite(values)) || any(values <= 0)) {
    stop(sprintf("`%s` must contain positive finite values.", label), call. = FALSE)
  }
  values
}

resolve_human_slot_contact_multiplier <- function(parameters, size) {
  values <- parameters$human_slot_contact_multiplier
  if (is.null(values)) {
    values <- 1
  }
  values <- validate_human_slot_contact_multiplier(values, size)
  values * legacy_contact_multiplier_for_human_slots(parameters)
}

human_slot_contact_multiplier_values <- function(variables) {
  if (is.null(variables$human_slot_contact_multiplier)) {
    return(NULL)
  }
  variables$human_slot_contact_multiplier$get_values()
}

human_biting_weights <- function(zeta, psi, human_slot_contact_multiplier = NULL) {
  weights <- zeta * psi
  if (!is.null(human_slot_contact_multiplier)) {
    human_slot_contact_multiplier <- validate_human_slot_contact_multiplier(
      human_slot_contact_multiplier,
      length(weights)
    )
    weights <- weights * human_slot_contact_multiplier
  }
  weights
}

human_pi <- function(zeta, psi, human_slot_contact_multiplier = NULL) {
  # human_pi() is conditional on a bite happening in this node, so these
  # weights must be normalized. The same slot multipliers enter the
  # unconditional node contact rate through human_slot_contact_rate_multiplier().
  weights <- human_biting_weights(zeta, psi, human_slot_contact_multiplier)
  weights / sum(weights)
}

human_slot_contact_rate_multiplier <- function(zeta, psi, human_slot_contact_multiplier = NULL) {
  if (is.null(human_slot_contact_multiplier)) {
    return(1)
  }
  base_weights <- zeta * psi
  human_slot_contact_multiplier <- validate_human_slot_contact_multiplier(
    human_slot_contact_multiplier,
    length(base_weights)
  )
  base_total <- sum(base_weights)
  if (!is.finite(base_total) || base_total <= 0) {
    stop("Base human biting weights must have positive finite total.", call. = FALSE)
  }
  sum(base_weights * human_slot_contact_multiplier) / base_total
}

blood_meal_rate <- function(v, z, parameters) {
  gonotrophic_cycle <- get_gonotrophic_cycle(v, parameters)
  interrupted_foraging_time <- parameters$foraging_time[[v]] / (1 - z)
  1 / (interrupted_foraging_time + gonotrophic_cycle)
}

human_blood_meal_rate <- function(species, variables, parameters, timestep) {
  age <- get_age(variables$birth$get_values(), timestep)
  psi <- unique_biting_rate(age, parameters)
  zeta <- variables$zeta$get_values()
  p_bitten <- prob_bitten(timestep, variables, species, parameters)
  human_slot_contact_multiplier <- human_slot_contact_multiplier_values(variables)
  slot_contact_rate_multiplier <- human_slot_contact_rate_multiplier(
    zeta,
    psi,
    human_slot_contact_multiplier
  )
  .pi <- human_pi(zeta, psi, human_slot_contact_multiplier)
  Q0 <- parameters$Q0[[species]]
  W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
  Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
  f <- blood_meal_rate(species, Z, parameters)
  apply_runtime_contact_multiplier(
    .human_blood_meal_rate(f, species, W, parameters) *
      slot_contact_rate_multiplier,
    species,
    parameters
  )
}

.human_blood_meal_rate <- function(f, v, W, parameters) {
  Q <- 1 - (1 - parameters$Q0[[v]]) / W
  Q * f
}

average_p_repelled <- function(p_repelled, .pi, Q0) {
  Q0 * sum(.pi * p_repelled)
}

average_p_successful <- function(prob_bitten_survives, .pi, Q0) {
  (1 - Q0) + Q0 * sum(.pi *  prob_bitten_survives)
}

# Unique biting rate (psi) for a human of a given age
unique_biting_rate <- function(age, parameters) {
  1 - parameters$rho * exp(- age / parameters$a0)
}

#' @title Calculate the force of infection towards mosquitoes
#'
#' @param a human blood meal rate
#' @param infectivity_sum the sum of each individual's infectivity 
#' @noRd
calculate_foim <- function(a, infectivity_sum) {
  a * infectivity_sum
}
