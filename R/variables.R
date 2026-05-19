#' @title Define model variables 
#' @description
#' create_variables creates the human and mosquito variables for
#' the model. Variables are used to track real data for each individual over
#' time, they are read and updated by processes
#'
#' The human variables are defined as:
#'
#' * state - the state of a human individual S|D|A|U|Tr|NonExistent
#' * birth - an integer representing the timestep when this individual was born
#' * last_boosted_* - the last timestep at which this individual's immunity was
#' boosted for tracking grace periods in the boost of immunity
#' * IAM - Maternal anti-parasite immunity (p.v only)
#' * ICM - Maternal immunity to clinical disease
#' * IVM - Maternal immunity to severe disease (p.f only)
#' * IB  - Pre-erythrocytic immunity (p.f only)
#' * IAA  - Acquired anti-parasite immunity (p.v only)
#' * ICA  - Acquired immunity to clinical disease
#' * IVA  - Acquired immunity to severe disease (p.f only)
#' * ID - Acquired immunity to detectability (p.f only)
#' * hypnozoites - Hypnozoite batch number (p.v only)
#' * zeta - Heterogeneity of human individuals
#' * zeta_group - Discretised heterogeneity of human individuals
#' * human_slot_contact_multiplier - Optional persistent human-slot contact
#' multiplier. Constant values scale the node's total contact rate;
#' heterogeneous values also redistribute contacts across human slots. This is
#' suitable for household/compound exposure heterogeneity represented on human
#' slots. Values are measured relative to the reference value 1 and survive
#' rebirth. Scalar legacy `contact_multiplier` values are folded into this
#' variable during initialization.
#' * last_pev_timestep - The timestep of the last pev vaccination (-1 if there
#' * last_eff_pev_timestep - The timestep of the last efficacious pev
#' vaccination, including final primary dose and booster doses (-1 if there have not been any)
#' * pev_profile - The index of the efficacy profile of any pev vaccinations.
#' Not set until the final primary dose.
#' This is only set on the final primary dose and subsequent booster doses
#' (-1 otherwise)
#' * tbv_vaccinated - The timstep of the last tbv vaccination (-1 if there
#' haven't been any
#' * net_time - The timestep when a net was last put up (-1 if never)
#' * spray_time - The timestep when the house was last sprayed (-1 if never)
#' * infectivity - The onward infectiousness to mosquitos
#' * drug - The last prescribed drug
#' * drug_time - The timestep of the last drug
#' * ls_drug - The last prescribed drug
#' * ls_drug_time - The timestep of the last drug
#'
#' Antimalarial resistance variables are:
#' * dt - the delay for humans to move from state Tr to state S
#'
#' Mosquito variables are: 
#' * mosquito_state - the state of the mosquito, a category Sm|Pm|Im|NonExistent
#' * species - the species of mosquito, this is a category gamb|fun|arab
#'
#' @param parameters, model parameters created by `get_parameters`
#' @noRd
#' @importFrom stats rexp rnorm
create_variables <- function(parameters) {
  size <- get_human_population(parameters, 0)
  states <- c('S', 'D', 'A', 'U', 'Tr')
  mode <- human_initialization_mode(parameters)
  stationary_initialization_context <- NULL

  if (mode == "stochastic_resample" && !is.null(parameters$init_EIR)) {
    library <- get_stationary_human_initialization_library(parameters)
    sampled <- stationary_human_initializer_resample(
      library,
      size = size,
      parameters = parameters
    )

    initial_age <- -sampled$birth
    state <- individual::CategoricalVariable$new(states, sampled$state)
    birth <- individual::IntegerVariable$new(as.integer(sampled$birth))
    last_boosted_ica <- individual::DoubleVariable$new(sampled$last_boosted_ica)
    icm <- individual::DoubleVariable$new(sampled$icm)
    ica <- individual::DoubleVariable$new(sampled$ica)
    zeta <- individual::DoubleVariable$new(sampled$zeta)
    zeta_group <- individual::CategoricalVariable$new(
      to_char_vector(seq(parameters$n_heterogeneity_groups)),
      sampled$zeta_group
    )
    if (is.null(sampled$human_slot_contact_multiplier)) {
      sampled$human_slot_contact_multiplier <- resolve_human_slot_contact_multiplier(
        parameters,
        size
      )
    }
    human_slot_contact_multiplier <- individual::DoubleVariable$new(
      validate_human_slot_contact_multiplier(
        sampled$human_slot_contact_multiplier,
        size
      )
    )
    infectivity_values <- sampled$infectivity
    progression_rate_values <- sampled$progression_rates
    drug <- individual::IntegerVariable$new(as.integer(sampled$drug))
    drug_time <- individual::IntegerVariable$new(as.integer(sampled$drug_time))
    last_pev_timestep <- individual::IntegerVariable$new(as.integer(sampled$last_pev_timestep))
    last_eff_pev_timestep <- individual::IntegerVariable$new(as.integer(sampled$last_eff_pev_timestep))
    pev_profile <- individual::IntegerVariable$new(as.integer(sampled$pev_profile))
    tbv_vaccinated <- individual::DoubleVariable$new(sampled$tbv_vaccinated)
    net_time <- individual::IntegerVariable$new(as.integer(sampled$net_time))
    spray_time <- individual::IntegerVariable$new(as.integer(sampled$spray_time))

    if (parameters$parasite == "falciparum") {
      last_boosted_ib <- individual::DoubleVariable$new(sampled$last_boosted_ib)
      last_boosted_iva <- individual::DoubleVariable$new(sampled$last_boosted_iva)
      last_boosted_id <- individual::DoubleVariable$new(sampled$last_boosted_id)
      ivm <- individual::DoubleVariable$new(sampled$ivm)
      ib <- individual::DoubleVariable$new(sampled$ib)
      iva <- individual::DoubleVariable$new(sampled$iva)
      id <- individual::DoubleVariable$new(sampled$id)
    } else if (parameters$parasite == "vivax") {
      last_boosted_iaa <- individual::DoubleVariable$new(sampled$last_boosted_iaa)
      iaa <- individual::DoubleVariable$new(sampled$iaa)
      iam <- individual::DoubleVariable$new(sampled$iam)
      hypnozoites <- individual::IntegerVariable$new(as.integer(sampled$hypnozoites))
      if (any(parameters$drug_hypnozoite_efficacy > 0)) {
        ls_drug <- individual::IntegerVariable$new(as.integer(sampled$ls_drug))
        ls_drug_time <- individual::IntegerVariable$new(as.integer(sampled$ls_drug_time))
      }
    }

    stationary_initialization_context <- stationary_human_initializer_sample_context(
      sampled = sampled,
      parameters = parameters,
      library = library
    )
  } else {
    initial_age <- calculate_initial_ages(parameters)

    if (parameters$enable_heterogeneity) {
      quads <- statmod::gauss.quad.prob(
        parameters$n_heterogeneity_groups,
        dist='normal'
      )
      groups <- sample.int(
        parameters$n_heterogeneity_groups,
        size,
        replace = TRUE,
        prob = quads$weights
      )
      zeta_norm <- quads$nodes[groups]
      zeta <- individual::DoubleVariable$new(
        calculate_zeta(zeta_norm, parameters)
      )
      zeta_group <- individual::CategoricalVariable$new(
        to_char_vector(seq(parameters$n_heterogeneity_groups)),
        to_char_vector(groups)
      )
      if (!is.null(parameters$init_EIR)) {
        eq <- calculate_eq(quads$nodes, parameters)
      } else {
        eq <- NULL
      }
    } else {
      zeta <- individual::DoubleVariable$new(rep(1, size))
      groups <- rep(1, size)
      zeta_group <- individual::CategoricalVariable$new(
        to_char_vector(seq(parameters$n_heterogeneity_groups)),
        to_char_vector(groups)
      )
      if (!is.null(parameters$init_EIR)) {
        if(parameters$parasite == "falciparum"){
          eq <- list(
            malariaEquilibrium::human_equilibrium_no_het(
              parameters$init_EIR,
              equilibrium_treatment_coverage(parameters),
              parameters$eq_params,
              EQUILIBRIUM_AGES
            )
          )
        } else if (parameters$parasite == "vivax"){
          eq <- malariaEquilibriumVivax::vivax_equilibrium(
            EIR = parameters$init_EIR,
            ft = equilibrium_treatment_coverage(parameters),
            p = translate_vivax_parameters(parameters),
            age = EQUILIBRIUM_AGES
          )$states
        }
      } else {
        eq <- NULL
      }
    }

    human_slot_contact_multiplier <- individual::DoubleVariable$new(
      resolve_human_slot_contact_multiplier(parameters, size)
    )

    if(parameters$parasite == "falciparum"){
      initial_states <- initial_state(parameters, initial_age, groups, eq, states)
      hypnozoite_v <- NULL

    } else if (parameters$parasite == "vivax"){

      eq_v_output <- initial_state_vivax(parameters, initial_age, groups, eq, states)

      # Human states
      initial_states <- eq_v_output$human_states
      hypnozoite_v <- eq_v_output$hypnozoites_v

      ## Initial hypnozoites
      hypnozoites <- individual::IntegerVariable$new(hypnozoite_v)

    }

    state <- individual::CategoricalVariable$new(states, initial_states)
    birth <- individual::IntegerVariable$new(-initial_age)

    # Maternal immunity to clinical disease
    icm <- individual::DoubleVariable$new(
      initial_immunity(
        parameters$init_icm,
        initial_age,
        groups,
        eq,
        parameters,
        'ICM',
        hypnozoite_v
      )
    )

    # Acquired immunity to clinical disease
    last_boosted_ica <- individual::DoubleVariable$new(rep(-1, size))
    ica <- individual::DoubleVariable$new(
      initial_immunity(
        parameters$init_ica,
        initial_age,
        groups,
        eq,
        parameters,
        'ICA',
        hypnozoite_v
      )
    )

    if(parameters$parasite == "falciparum"){
      # Pre-erythoctic immunity
      last_boosted_ib <- individual::DoubleVariable$new(rep(-1, size))
      ib  <- individual::DoubleVariable$new(
        initial_immunity(
          parameters$init_ib,
          initial_age,
          groups,
          eq,
          parameters,
          'IB'
        )
      )

      # Maternal immunity to severe disease
      ivm <- individual::DoubleVariable$new(
        initial_immunity(
          parameters$init_ivm,
          initial_age,
          groups,
          eq,
          parameters,
          'IVM'
        )
      )

      # Acquired immunity to severe disease
      last_boosted_iva <- individual::DoubleVariable$new(rep(-1, size))
      iva <- individual::DoubleVariable$new(
        initial_immunity(
          parameters$init_iva,
          initial_age,
          groups,
          eq,
          parameters,
          'IVA'
        )
      )

      # Acquired immunity to lm detectability
      last_boosted_id <- individual::DoubleVariable$new(rep(-1, size))
      id <- individual::DoubleVariable$new(
        initial_immunity(
          parameters$init_id,
          initial_age,
          groups,
          eq,
          parameters,
          'ID'
        )
      )

    } else if (parameters$parasite == "vivax"){
      # Maternal anti-parasite immunity
      iam <- individual::DoubleVariable$new(
        initial_immunity(
          parameters$init_iam,
          initial_age,
          groups,
          eq,
          parameters,
          'IAM',
          hypnozoite_v
        )
      )

      # Acquired anti-parasite immunity
      last_boosted_iaa <- individual::DoubleVariable$new(rep(-1, size))
      iaa <- individual::DoubleVariable$new(
        initial_immunity(
          parameters$init_iaa,
          initial_age,
          groups,
          eq,
          parameters,
          'IAA',
          hypnozoite_v
        )
      )
    }

    # Initialise infectiousness of humans -> mosquitoes
    # NOTE: not yet supporting initialisation of infectiousness of Treated individuals
    infectivity_values <- rep(0, get_human_population(parameters, 0))

    # Calculate the indices of individuals in each infectious state
    diseased <- state$get_index_of('D')$to_vector()
    asymptomatic <- state$get_index_of('A')$to_vector()
    subpatent <- state$get_index_of('U')$to_vector()
    treated <- state$get_index_of('Tr')$to_vector()

    # Set the initial infectivity values for each individual
    infectivity_values[diseased] <- parameters$cd
    if(parameters$parasite == "falciparum"){
      # p.f has immunity-determined asymptomatic infectivity
      infectivity_values[asymptomatic] <- asymptomatic_infectivity(
        initial_age[asymptomatic],
        id$get_values(asymptomatic),
        parameters
      )
    } else if (parameters$parasite == "vivax"){
      # p.v has constant asymptomatic infectivity
      infectivity_values[asymptomatic] <- parameters$ca
    }
    infectivity_values[subpatent] <- parameters$cu

    # Set disease progression rates for each individual
    progression_rate_values <- rep(0, get_human_population(parameters, 0))
    progression_rate_values[diseased] <- 1/parameters$dd
    progression_rate_values[asymptomatic] <- 1/parameters$da
    if(parameters$parasite == "falciparum"){
      # p.f subpatent recovery rate is constant
      progression_rate_values[subpatent] <- 1/parameters$du
    } else if (parameters$parasite == "vivax"){
      # p.v subpatent recovery rate is immunity-dependent
      progression_rate_values[subpatent] <- 1/anti_parasite_immunity(
        parameters$dpcr_min, parameters$dpcr_max, parameters$apcr50, parameters$kpcr,
        iaa$get_values(subpatent),
        iam$get_values(subpatent)
      )
    }
    progression_rate_values[treated] <- 1/parameters$dt

    drug <- individual::IntegerVariable$new(rep(0, size))
    drug_time <- individual::IntegerVariable$new(rep(-1, size))

    if(any(parameters$drug_hypnozoite_efficacy > 0)){
      ls_drug <- individual::IntegerVariable$new(rep(0, size))
      ls_drug_time <- individual::IntegerVariable$new(rep(-1, size))
    }

    last_pev_timestep <- individual::IntegerVariable$new(rep(-1, size))
    last_eff_pev_timestep <- individual::IntegerVariable$new(rep(-1, size))
    pev_profile <- individual::IntegerVariable$new(rep(-1, size))

    tbv_vaccinated <- individual::DoubleVariable$new(rep(-1, size))

    # Init vector controls. This lets a background bednet regime start in-place
    # rather than arriving as an avoidable day-1 intervention shock.
    net_time_values <- rep(-1L, size)
    if (isTRUE(parameters$bednets) && !is.null(parameters$initial_bednet_coverage)) {
      coverage <- as.numeric(parameters$initial_bednet_coverage)
      if (!is.finite(coverage) || length(coverage) != 1L || coverage < 0 || coverage > 1) {
        stop("initial_bednet_coverage must be a single finite number in [0, 1].")
      }
      initial_net_time <- parameters$initial_bednet_time
      if (!is.numeric(initial_net_time) || length(initial_net_time) != 1L ||
          !is.finite(initial_net_time)) {
        stop("initial_bednet_time must be a single finite numeric timestep.")
      }
      n_init <- as.integer(round(coverage * size))
      if (n_init > 0L) {
        target <- sample.int(size, size = n_init, replace = FALSE)
        net_time_values[target] <- as.integer(round(initial_net_time))
      }
    }
    net_time <- individual::IntegerVariable$new(net_time_values)
    spray_time <- individual::IntegerVariable$new(rep(-1, size))
  }

  progression_rates <- individual::DoubleVariable$new(progression_rate_values)

  if (mode != "stochastic_resample") {
    infectivity_values <- match_initial_human_infectivity(
      infectivity_values = infectivity_values,
      birth = birth,
      zeta = zeta,
      human_slot_contact_multiplier = human_slot_contact_multiplier,
      net_time = net_time,
      spray_time = spray_time,
      parameters = parameters
    )
  }
  infectivity <- individual::DoubleVariable$new(infectivity_values)

  variables <- list(
    state = state,
    birth = birth,
    last_boosted_ica = last_boosted_ica,
    icm = icm,
    ica = ica,
    zeta = zeta,
    zeta_group = zeta_group,
    infectivity = infectivity,
    progression_rates = progression_rates,
    drug = drug,
    drug_time = drug_time,
    last_pev_timestep = last_pev_timestep,
    last_eff_pev_timestep = last_eff_pev_timestep,
    pev_profile = pev_profile,
    tbv_vaccinated = tbv_vaccinated,
    net_time = net_time,
    spray_time = spray_time
  )
  
  if(parameters$parasite == "falciparum"){
    variables <- c(variables,
                   last_boosted_ib = last_boosted_ib,
                   last_boosted_iva = last_boosted_iva,
                   last_boosted_id = last_boosted_id,
                   ivm = ivm,
                   ib = ib,
                   iva = iva,
                   id = id
    )
  } else if (parameters$parasite == "vivax"){
    variables <- c(variables,
                   last_boosted_iaa = last_boosted_iaa,
                   iaa = iaa,
                   iam = iam,
                   hypnozoites = hypnozoites
    )
  }

  variables <- c(
    variables,
    human_slot_contact_multiplier = human_slot_contact_multiplier
  )
  
  # Add variables for individual mosquitoes
  if (legacy_individual_mosquito_backend_enabled(parameters)) {
    species_values <- NULL
    state_values <- NULL
    n_initialised <- 0
    for (i in seq_along(parameters$species)) {
      mosquito_counts <- floor(
        initial_mosquito_counts(
          parameters,
          i,
          parameters$init_foim,
          parameters$total_M * parameters$species_proportions[[i]]
        )
      )

      species_M <- sum(mosquito_counts[ADULT_ODE_INDICES])

      if (species_M > 0) {
        if (length(species_values) > parameters$mosquito_limit) {
          stop('Mosquito limit not high enough')
        }

        species_values <- c(
          species_values,
          rep(parameters$species[[i]], species_M)
        )
        state_values <- c(
          state_values,
          rep(
            c('Sm', 'Pm', 'Im'),
            times = mosquito_counts[ADULT_ODE_INDICES]
          )
        )
      }
    }

    # fill excess mosquitoes
    excess <- parameters$mosquito_limit - length(species_values)
    species_values <- c(
      species_values,
      rep(parameters$species[[1]], excess)
    )
    state_values <- c(state_values, rep('NonExistent', excess))

    # initialise variables
    species <- individual::CategoricalVariable$new(
      parameters$species,
      species_values
    )
    mosquito_state <- individual::CategoricalVariable$new(
      c('Sm', 'Pm', 'Im', 'NonExistent'),
      state_values
    )
    wt_genotype <- 1L
    if (!is.null(parameters$cube)) {
      wt_genotype <- cube_wild_type_index(parameters$cube)
    }
    geno_id <- individual::IntegerVariable$new(
      rep.int(wt_genotype, parameters$mosquito_limit)
    )
    variables <- c(
      variables,
      species = species,
      mosquito_state = mosquito_state,
      geno_id = geno_id
    )
  }

  if (!is.null(stationary_initialization_context)) {
    attr(variables, "stationary_initialization_context") <- stationary_initialization_context
  }

  variables
}


create_export_variable <- function(metapop_params) {
  individual::DoubleVariable$new(rep(0, length(metapop_params$x)))
}

# =========
# Utilities
# =========

match_initial_human_infectivity <- function(
  infectivity_values,
  birth,
  zeta,
  human_slot_contact_multiplier = NULL,
  net_time,
  spray_time,
  parameters
) {
  if (is.null(parameters$init_foim) || !is.finite(parameters$init_foim) ||
      parameters$init_foim <= 0) {
    return(infectivity_values)
  }

  age <- get_age(birth$get_values(), 0)
  psi <- unique_biting_rate(age, parameters)
  human_slot_contact_values <- if (is.null(human_slot_contact_multiplier)) {
    NULL
  } else {
    human_slot_contact_multiplier$get_values()
  }
  .pi <- human_pi(zeta$get_values(), psi, human_slot_contact_values)
  current_infectivity <- sum(.pi * infectivity_values)
  if (!is.finite(current_infectivity) || current_infectivity <= 0) {
    return(infectivity_values)
  }

  vc_variables <- list(
    birth = birth,
    zeta = zeta,
    human_slot_contact_multiplier = human_slot_contact_multiplier,
    net_time = net_time,
    spray_time = spray_time
  )
  a0 <- human_blood_meal_rate(1, vc_variables, parameters, 0)
  current_foim <- a0 * current_infectivity
  if (!is.finite(current_foim) || current_foim <= 0) {
    return(infectivity_values)
  }

  scale <- parameters$init_foim / current_foim
  if (!is.finite(scale) || scale <= 0 || abs(scale - 1) < 1e-8) {
    return(infectivity_values)
  }

  infectious <- infectivity_values > 0
  if (!any(infectious)) {
    return(infectivity_values)
  }

  infectivity_values[infectious] <- infectivity_values[infectious] * scale
  infectivity_values
}

initial_immunity <- function(
  parameter,
  age,
  groups = NULL,
  eq = NULL,
  parameters = NULL,
  eq_name = NULL,
  hyp = NULL
  ) {
  if (!is.null(eq)) {
    age <- age / 365
    return(vnapply(
      seq_along(age),
      function(i) {
        g <- groups[[i]]
        a <- age[[i]]
        if(parameters$parasite == "falciparum"){
          eq[[g]][which.max(a < eq[[g]][, 'age']), eq_name]
        } else if (parameters$parasite == "vivax"){
          h <- hyp[[i]]
          eq[[eq_name]][which.max(a < eq$Age[-1]), g, h+1]
        }
      }
    ))
  }
  rep(parameter, length(age))
}

initial_state <- function(parameters, age, groups, eq, states) {
  ibm_states <- states
  if (!is.null(eq)) {
    eq_states <- c('S', 'D', 'A', 'U', 'T')
    age <- age / 365
    return(vcapply(
      seq_along(age),
      function(i) {
        g <- groups[[i]]
        a <- age[[i]]
        sample(
          ibm_states,
          size = 1,
          prob = eq[[g]][which.max(a < eq[[g]][, 'age']), eq_states]
        )
      }
    ))
  }
  rep(ibm_states, times = calculate_initial_counts(parameters))
}

initial_state_vivax <- function(parameters, age, groups, eq, states) {
  # vivax human states and hypnozoites must be calculated over a combined probability distribution
  ibm_states <- states
  if (!is.null(eq)) {
    eq_states <- c('S', 'D', 'A', 'U', 'T')
    age <- age / 365
    human_states_pop <- sapply(
      seq_along(age),
      function(i) {
        g <- groups[[i]]
        a <- age[[i]]
        human_states <- c(expand.grid("hyp" = 0:parameters$kmax, "human_state" = eq_states))
        probs <- c(sapply(eq_states, function(state){eq[[state]][which.max(a < eq$Age[-1]), g, ]}))
        smp <- sample(x = 1:length(probs), size = 1, replace = T, prob = probs)
        return(c(human_states$human_state[smp],human_states$hyp[smp]))
      }
    )
    return(list(human_states = states[human_states_pop[1,]],
                hypnozoites_v = human_states_pop[2,]))
  }
  return(list(human_states = rep(ibm_states, calculate_initial_counts(parameters)),
              hypnozoites_v = rep(parameters$init_hyp, parameters$human_population)))
}

calculate_initial_counts <- function(parameters) {
  pop <- get_human_population(parameters, 0)
  initial_counts <- round(
    c(
      parameters$s_proportion,
      parameters$d_proportion,
      parameters$a_proportion,
      parameters$u_proportion,
      parameters$t_proportion
    ) * pop
  )
  left_over <- pop - sum(initial_counts)
  initial_counts[[1]] <- initial_counts[[1]] + left_over
  initial_counts
}

calculate_eq <- function(het_nodes, parameters) {
  ft <- equilibrium_treatment_coverage(parameters)
  if(parameters$parasite == "falciparum"){
    lapply(
      het_nodes,
      function(n) {
        malariaEquilibrium::human_equilibrium_no_het(
          parameters$init_EIR * calculate_zeta(n, parameters),
          ft,
          parameters$eq_params,
          EQUILIBRIUM_AGES
        )
      }
    )
  } else if (parameters$parasite == "vivax"){
    eq <- malariaEquilibriumVivax::vivax_equilibrium(
      EIR = parameters$init_EIR,
      ft = equilibrium_treatment_coverage(parameters),
      age = EQUILIBRIUM_AGES,
      p = translate_vivax_parameters(parameters)
    )$states
  }
}

calculate_zeta <- function(zeta_norm, parameters) {
  exp(
    zeta_norm * sqrt(parameters$sigma_squared) - parameters$sigma_squared/2
  )
}

calculate_initial_ages <- function(parameters) {
  n_pop <- get_human_population(parameters, 0)
  # check if we've set up a custom demography
  if (!parameters$custom_demography) {
    return(round(rtexp(
      n_pop,
      1 / parameters$average_age,
      max(EQUILIBRIUM_AGES)*365
    )))
  }

  deathrates <- parameters$deathrates[1, , drop = FALSE]
  age_high <- parameters$deathrate_agegroups
  age_width <- diff(c(0, age_high))
  age_low <- age_high - age_width
  n_age <- length(age_high)
  birthrate <- find_birthrates(parameters$human_population, age_high, deathrates)
  deathrates <- parameters$deathrates[1,]

  eq_pop <- get_equilibrium_population(age_high, birthrate, deathrates)

  group <- sample.int(
    n_age,
    n_pop,
    replace = TRUE,
    prob = eq_pop
  )

  # sample truncated exponential for each age group
  ages <- rep(NA, n_pop)
  for (g in seq(n_age)) {
    in_group <- group == g
    group_ages <- rtexp(sum(in_group), deathrates[[g]], age_width[[g]])
    ages[in_group] <- age_low[[g]] + group_ages
  }

  ages
}
