#' @title Calculate equilibrium solution for vector counts
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param species the index of the species to find the equilibrium for
#' @param foim equilibrium foim
#' @param m the total number of female adult mosquitos
#' @noRd
initial_mosquito_counts <- function(parameters, species, foim, m) {
  omega <- calculate_omega(parameters, species)
  traits <- equilibrium_species_traits(parameters, species)
  mum <- traits$muF
  n_E <- 2 * omega * mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_L <- 2 * mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_P <- 2 * parameters$dpl * mum * m

  n_Sm <- m * mum / (foim + mum)

  incubation_survival <- exp(-mum * parameters$dem)

  n_Pm <- m * foim / (foim + mum) * (
    1. - incubation_survival
  )

  n_Im <- m * foim / (foim + mum) * incubation_survival

  c(n_E, n_L, n_P, n_Sm, n_Pm, n_Im)
}

#' @title Calculate omega value
#' @description useful value for calculating equilibrium solutions for vectors
#' taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param species the index of the species to calculate for
#' @noRd
calculate_omega <- function(parameters, species) {
  sub_omega <- parameters$gamma * parameters$ml / parameters$me - (
    parameters$del / parameters$dl
  ) + (
    (parameters$gamma - 1) * parameters$ml * parameters$del
  )

  traits <- equilibrium_species_traits(parameters, species)
  mum <- traits$muF

  beta <- eggs_laid(
    species_beta_value(parameters, species),
    mum,
    traits$f
  )

  -.5 * sub_omega + sqrt(
    .25 * sub_omega**2 +
      .5 * parameters$gamma * beta * parameters$ml * parameters$del /
      (parameters$me * mum * parameters$dl * (
        1. + parameters$dpl * parameters$mup
      ))
  )
}

#' @title Calculate the vector carrying capacity
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param m number of adult mosquitoes
#' @param species index of the species to calculate for
calculate_carrying_capacity <- function(parameters, m, species) {
  omega <- calculate_omega(parameters, species)
  traits <- equilibrium_species_traits(parameters, species)

  m * 2 * parameters$dl * traits$muF * (
    1. + parameters$dpl * parameters$mup
  ) * parameters$gamma * (omega + 1) / (
    omega / (parameters$ml * parameters$del) - (
      1. / (parameters$ml * parameters$dl)
    ) - 1.
  )
}

#' @title Calculate the mean rainfall throughout the year
#' @param parameters model parameters
#' @noRd
calculate_R_bar <- function(parameters) {
  mean(vnapply(1:365, function(t) rainfall(
		t,
    parameters$g0,
    parameters$g,
    parameters$h,
    parameters$rainfall_floor
	)))
}

#' @title Calculate equilibrium total_M from parameters
#'
#' @param parameters to work from
#' @param EIR equilibrium to use, bites per person per year
#' @importFrom stats weighted.mean
#' @noRd
equilibrium_total_M <- function(parameters, EIR) {
  if (EIR == 0) {
    return(0)
  }
  if (parameters$init_foim == 0) {
    stop('init_foim must be > 0 to calculate a non-zero equilibrium total_M')
  }
  total_daily_eir <- EIR * parameters$human_population / 365

  bednet_background <- equilibrium_bednet_background(parameters)
  if (!is.null(bednet_background)) {
    lifetime <- parameters$init_foim * exp(-bednet_background$muF * parameters$dem) / (
      parameters$init_foim + bednet_background$muF
    )
    return(total_daily_eir / sum(
      parameters$species_proportions * bednet_background$a * lifetime
    ))
  }

  mum <- weighted.mean(parameters$mum, parameters$species_proportions)
  lifetime <- parameters$init_foim * exp(-mum * parameters$dem) / (
    parameters$init_foim + mum
  )
  total_daily_eir / sum(
    parameters$species_proportions * parameters$blood_meal_rates * parameters$Q0 * lifetime
  )
}

equilibrium_bednet_background <- function(parameters) {
  if (!isTRUE(parameters$bednets) ||
      isTRUE(parameters$spraying) ||
      is.null(parameters$initial_bednet_coverage) ||
      is.null(parameters$initial_bednet_time)) {
    return(NULL)
  }

  coverage <- as.numeric(parameters$initial_bednet_coverage)
  if (!is.finite(coverage) || length(coverage) != 1L || coverage <= 0 || coverage > 1) {
    return(NULL)
  }

  bednet_time <- as.integer(round(parameters$initial_bednet_time))
  match_idx <- match(bednet_time, as.integer(parameters$bednet_timesteps))
  if (is.na(match_idx)) {
    return(NULL)
  }

  n_species <- length(parameters$species)
  muF <- a <- f <- numeric(n_species)
  for (species in seq_len(n_species)) {
    phi_bednets <- parameters$phi_bednets[[species]]
    rn <- prob_repelled_bednets(match_idx, 0, species, parameters)
    sn <- prob_survives_bednets(rn, match_idx, 0, species, parameters)

    mean_prob_bitten_survives <- (
      (1 - coverage) +
      coverage * (1 - phi_bednets + phi_bednets * sn)
    )
    mean_prob_repelled <- coverage * phi_bednets * rn

    Q0 <- parameters$Q0[[species]]
    W <- (1 - Q0) + Q0 * mean_prob_bitten_survives
    Z <- Q0 * mean_prob_repelled
    f[[species]] <- blood_meal_rate(species, Z, parameters)
    a[[species]] <- .human_blood_meal_rate(f[[species]], species, W, parameters)
    muF[[species]] <- death_rate(f[[species]], W, Z, species, parameters)
  }

  list(
    coverage = coverage,
    muF = muF,
    a = a,
    f = f
  )
}

equilibrium_species_traits <- function(parameters, species) {
  bg <- equilibrium_bednet_background(parameters)
  if (is.null(bg)) {
    return(list(
      muF = parameters$mum[[species]],
      a = parameters$blood_meal_rates[[species]] * parameters$Q0[[species]],
      f = parameters$blood_meal_rates[[species]]
    ))
  }

  list(
    muF = bg$muF[[species]],
    a = bg$a[[species]],
    f = bg$f[[species]]
  )
}

#' @title Calculate the yearly offset (in timesteps) for the peak mosquito
#' season
#'
#' @param parameters to work from
#' @export
peak_season_offset <- function(parameters) {
  if (!parameters$model_seasonality) {
    return(0)
  }
  which.max(vnapply(seq(365), function(t) {
    rainfall(
      t,
      parameters$g0,
      parameters$g,
      parameters$h,
      parameters$rainfall_floor
    )
  }))[[1]]
}

#' @title Calculate the death rate of mosquitoes given interventions
#'
#' @param f the feeding rate for this species of mosquito
#' @param W the mean probability that a mosquito feeds and survives
#' @param Z the mean probability that a mosquito is repelled
#' @param Z the mean probability that a mosquito is repelled
#' @noRd
death_rate <- function(f, W, Z, species, parameters) {
  mum <- parameters$mum[[species]]
  p1_0 <- exp(-mum * parameters$foraging_time[[species]])
  gonotrophic_cycle <- get_gonotrophic_cycle(species, parameters)
  p2 <- exp(-mum * gonotrophic_cycle)
  p1 <- p1_0 * W / (1 - Z * p1_0)
  -f * log(p1 * p2)
}

get_gonotrophic_cycle <- function(v, parameters) {
  f <- parameters$blood_meal_rates[[v]]
  gonotrophic_cycle <- 1 / f - parameters$foraging_time[[v]]
}

#' @title Update the individual mosquito model after biting
#'
#' @param variables a list of variables in this simulation
#' @param foim force of infection towards mosquitoes
#' @param events events in the simulation
#' @param species the index of the species to calculate for
#' @param susceptible_species the indices of susceptible mosquitos of the
#' species
#' @param adult_species the indices of adult mosquitos of the species
#' @param mu the death rate of the current species
#' @param parameters the model parameters
#' @param timestep the current timestep
#' @noRd
biting_effects_individual <- function(
    variables,
    foim,
    events,
    species,
    susceptible_species,
    adult_species,
    mu,
    parameters,
    timestep,
    omega_by_genotype = NULL
  ) {
  # deal with mosquito infections
  target <- sample_bitset(susceptible_species, foim)
  variables$mosquito_state$queue_update('Pm', target)
  events$mosquito_infection$schedule(
    target,
    log_uniform(target$size(), parameters$dem)
  )

  # deal with mosquito deaths
  if (is.null(omega_by_genotype) || length(omega_by_genotype) == 0L) {
    died <- sample_bitset(adult_species, mu)
  } else if (length(omega_by_genotype) == 1L || is.null(variables$geno_id)) {
    died <- sample_bitset(adult_species, mu * omega_by_genotype[[1]])
  } else {
    adult_genotypes <- variables$geno_id$get_values(adult_species)
    died_idx <- integer(0)
    for (g in seq_along(omega_by_genotype)) {
      g_pos <- which(adult_genotypes == g)
      if (length(g_pos) == 0L) {
        next
      }
      died_g <- sample_bitset(
        bitset_at(adult_species, g_pos),
        mu * omega_by_genotype[[g]]
      )$to_vector()
      if (length(died_g) > 0L) {
        died_idx <- c(died_idx, died_g)
      }
    }
    died <- individual::Bitset$new(parameters$mosquito_limit)
    if (length(died_idx) > 0L) {
      died$insert(died_idx)
    }
  }

  events$mosquito_death$schedule(died, 0)
}

cube_genotype_info <- function(cube) {
  if (is.null(cube)) {
    return(list(
      G = 1L,
      genotypesID = "WT",
      wild_type_index = 1L
    ))
  }

  if (is.null(cube$ih)) {
    stop("cube$ih must be provided")
  }
  ih_dim <- dim(cube$ih)
  if (length(ih_dim) != 3) {
    stop("cube$ih must be a 3D array")
  }
  G <- ih_dim[[3]]
  genotypes_id <- cube$genotypesID
  if (is.null(genotypes_id)) {
    genotypes_id <- as.character(seq_len(G))
  }
  if (length(genotypes_id) != G) {
    stop("length(cube$genotypesID) must match dim(cube$ih)[3]")
  }

  list(
    G = as.integer(G),
    genotypesID = genotypes_id,
    wild_type_index = cube_wild_type_index(cube)
  )
}

cube_wild_type_index <- function(cube) {
  if (is.null(cube) || is.null(cube$wildType)) {
    return(1L)
  }
  wt <- cube$wildType
  if (is.numeric(wt) && length(wt) >= 1) {
    return(as.integer(wt[[1]]))
  }
  if (is.character(wt) && length(wt) >= 1) {
    if (is.null(cube$genotypesID)) {
      stop("cube$genotypesID is required when cube$wildType is a character")
    }
    idx <- match(wt[[1]], cube$genotypesID)
    if (is.na(idx)) {
      stop("cube$wildType was not found in cube$genotypesID")
    }
    return(as.integer(idx))
  }
  1L
}

sample_genotype_counts <- function(n, p) {
  if (n <= 0) {
    return(integer(length(p)))
  }
  if (length(p) == 1L) {
    # exact and RNG-free in the trivial one-genotype case
    out <- rep.int(0L, length(p))
    out[[which.max(p)]] <- as.integer(n)
    return(out)
  }
  as.integer(stats::rmultinom(1, size = n, prob = p)[, 1])
}

sample_genotype_ids <- function(n, p) {
  counts <- sample_genotype_counts(n, p)
  rep.int(seq_along(counts), counts)
}

adult_female_genotype_counts_by_species <- function(variables, species_name, G) {
  adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
  species_index <- variables$species$get_index_of(species_name)$and(adult_index)
  if (species_index$size() == 0) {
    return(rep.int(0, G))
  }
  tabulate(variables$geno_id$get_values(species_index), nbins = G)
}

uses_genotype_resolved_aquatic_states <- function(solver_states, cube) {
  if (is.null(cube)) {
    return(FALSE)
  }
  cube_info <- cube_genotype_info(cube)
  length(solver_states) == 3L * cube_info$G && cube_info$G > 1L
}

aquatic_genotype_stage_indices <- function(G, stage) {
  stage_offsets <- c(E = 1L, L = 2L, P = 3L)
  if (!(stage %in% names(stage_offsets))) {
    stop("stage must be one of 'E', 'L', or 'P'")
  }
  seq.int(stage_offsets[[stage]], by = 3L, length.out = G)
}

aquatic_stage_values_by_genotype <- function(solver_states, cube, stage) {
  if (is.null(cube)) {
    stop("cube is required for genotype-resolved aquatic stage values")
  }
  cube_info <- cube_genotype_info(cube)
  if (!uses_genotype_resolved_aquatic_states(solver_states, cube)) {
    out <- rep.int(0, cube_info$G)
    out[[cube_info$wild_type_index]] <- solver_states[[ODE_INDICES[[stage]]]]
    names(out) <- cube_info$genotypesID
    return(out)
  }
  out <- solver_states[aquatic_genotype_stage_indices(cube_info$G, stage)]
  names(out) <- cube_info$genotypesID
  out
}

aquatic_stage_total <- function(solver_states, cube, stage) {
  if (is.null(cube) || !uses_genotype_resolved_aquatic_states(solver_states, cube)) {
    return(solver_states[[ODE_INDICES[[stage]]]])
  }
  sum(aquatic_stage_values_by_genotype(solver_states, cube, stage))
}

aquatic_stage_totals <- function(solver_states, cube) {
  c(
    E = aquatic_stage_total(solver_states, cube, "E"),
    L = aquatic_stage_total(solver_states, cube, "L"),
    P = aquatic_stage_total(solver_states, cube, "P")
  )
}

cube_phi_vector <- function(cube, G) {
  phi <- NULL
  if (!is.null(cube)) {
    phi <- cube$phi
  }
  if (is.null(phi)) {
    return(rep(0.5, G))
  }
  phi <- as.numeric(phi)
  if (length(phi) == 1L) {
    phi <- rep(phi, G)
  }
  if (length(phi) != G) {
    stop("cube$phi must have length 1 or length(cube$genotypesID)")
  }
  if (any(is.na(phi)) || any(phi < 0) || any(phi > 1)) {
    stop("cube$phi entries must be in [0, 1]")
  }
  phi
}

cube_omega_vector <- function(cube, G, genotypesID = NULL) {
  omega <- NULL
  if (!is.null(cube)) {
    omega <- cube$omega
  }
  if (is.null(omega)) {
    return(rep(1, G))
  }

  if (is.null(genotypesID) && !is.null(cube)) {
    genotypesID <- cube_genotype_info(cube)$genotypesID
  }

  if (length(omega) == 1L) {
    out <- rep(as.numeric(omega[[1]]), G)
    if (!is.null(genotypesID) && length(genotypesID) == G) {
      names(out) <- genotypesID
    }
  } else {
    nm <- names(omega)
    if (is.null(nm) || anyNA(nm) || any(nm == "")) {
      stop("cube$omega must be named with cube$genotypesID when length > 1")
    }
    if (anyDuplicated(nm)) {
      dup <- unique(nm[duplicated(nm)])
      stop(sprintf("cube$omega has duplicate names: %s", paste(dup, collapse = ", ")))
    }
    if (is.null(genotypesID)) {
      stop("cube$omega requires genotype names from cube$genotypesID")
    }
    missing_names <- setdiff(genotypesID, nm)
    if (length(missing_names) > 0) {
      stop(sprintf(
        "cube$omega is missing genotype names in cube$genotypesID: %s",
        paste(missing_names, collapse = ", ")
      ))
    }
    extra_names <- setdiff(nm, genotypesID)
    if (length(extra_names) > 0) {
      warning(sprintf(
        "cube$omega has extra genotype names that will be ignored: %s",
        paste(extra_names, collapse = ", ")
      ))
    }
    out <- as.numeric(omega[genotypesID])
    names(out) <- genotypesID
  }

  if (length(out) != G) {
    stop("cube$omega must have length 1 or length(cube$genotypesID)")
  }
  if (any(!is.finite(out)) || any(out < 0)) {
    stop("cube$omega entries must be finite and >= 0")
  }

  out
}

#' @title Calculate offspring genotype proportions and viability from cube
#' @param cube MGDrivE-style inheritance cube
#' @param female_counts adult female counts by genotype
#' @param male_counts adult male counts by genotype
#' @noRd
calc_pg_V_from_cube <- function(cube, female_counts, male_counts) {
  G <- length(female_counts)
  if (length(male_counts) != G) {
    stop("female_counts and male_counts must have the same length")
  }
  wt <- cube_wild_type_index(cube)
  p_fallback <- rep.int(0, G)
  p_fallback[[wt]] <- 1
  B_zero <- rep.int(0, G)

  if (is.null(cube) || is.null(cube$ih)) {
    return(list(p = p_fallback, V = 1, B = B_zero))
  }

  ih <- cube$ih
  ih_dim <- dim(ih)
  if (length(ih_dim) != 3 || any(ih_dim != c(G, G, G))) {
    stop("cube$ih dimensions must be G x G x G")
  }
  tau <- cube$tau
  if (is.null(tau)) {
    tau <- array(1, dim = ih_dim)
  } else if (!all(dim(tau) == ih_dim)) {
    stop("cube$tau dimensions must match cube$ih")
  }
  eta <- cube$eta
  if (is.null(eta)) {
    eta <- matrix(1, nrow = G, ncol = G)
  }
  if (!all(dim(eta) == c(G, G))) {
    stop("cube$eta dimensions must be G x G")
  }

  total_males <- sum(male_counts)
  if (total_males <= 0) {
    return(list(p = p_fallback, V = 1, B = B_zero))
  }

  Q <- outer(female_counts, male_counts / total_males)
  B <- vnapply(seq_len(G), function(g) {
    sum(Q * ih[, , g] * tau[, , g] * eta)
  })
  total_B <- sum(B)
  if (total_B > 0) {
    p <- B / total_B
  } else {
    p <- p_fallback
  }

  den <- sum(Q * apply(ih, c(1, 2), sum) * eta)
  V <- if (den > 0) total_B / den else 1

  list(p = p, V = V, B = B)
}

#' @title Mosquito emergence process
#' @description Move mosquitos from NonExistent to Sm in line with the number of
#' pupals in the ODE models
#'
#' @param solvers a list of solver objects for each species of mosquito
#' @param models mosquito model objects (used for genotype state when cube is set)
#' @param variables simulation variables (mosquito_state, species, geno_id)
#' @param parameters model parameters
#' @noRd
create_mosquito_emergence_process <- function(
  solvers,
  models,
  variables,
  parameters,
  dpl = NULL
  ) {
  if (!is.list(parameters) && !is.null(dpl)) {
    state <- models
    species <- variables
    species_names <- parameters
    rate <- .5 * 1 / dpl

    return(function(timestep) {
      p_counts <- vnapply(seq_along(solvers), function(i) solvers[[i]]$get_states()[[ODE_INDICES[["P"]]]])
      n <- sum(p_counts * rate)
      available <- state$get_size_of("NonExistent")
      if (n > available) {
        stop(paste0(
          "Not enough mosquitoes (short by ",
          n - available,
          "). Please raise parameters$mosquito_limit. ",
          "If you have used parameterise_mosquito_equilibrium,",
          "your seasonality parameters lead to more mosquitoes than expected."
        ))
      }
      non_existent <- state$get_index_of("NonExistent")
      latest <- 1
      for (i in seq_along(species_names)) {
        to_hatch <- p_counts[[i]] * rate
        hatched <- bitset_at(non_existent, seq(latest, latest + to_hatch))
        state$queue_update("Sm", hatched)
        species$queue_update(species_names[[i]], hatched)
        latest <- latest + to_hatch + 1
      }
    })
  }

  if (native_mosquito_backend_enabled(parameters)) {
    return(function(timestep) {})
  }
  state <- variables$mosquito_state
  species <- variables$species
  geno_id <- variables$geno_id
  species_names <- parameters$species
  rate <- .5 * 1 / parameters$dpl
  function(timestep) {
    solver_rows <- lapply(solvers, function(solver) solver$get_states())
    p_counts <- vnapply(
      seq_along(solvers),
      function(i) {
        aquatic_stage_total(solver_rows[[i]], models[[i]]$cube, "P")
      }
    )
    n <- 0
    for (i in seq_along(species_names)) {
      if (!is.null(models[[i]]$cube) &&
          uses_genotype_resolved_aquatic_states(solver_rows[[i]], models[[i]]$cube)) {
        cube_info <- cube_genotype_info(models[[i]]$cube)
        p_by_g <- aquatic_stage_values_by_genotype(solver_rows[[i]], models[[i]]$cube, "P")
        phi <- cube_phi_vector(models[[i]]$cube, cube_info$G)
        n_total_by_g <- pmax(0L, as.integer(round(as.numeric(p_by_g) / parameters$dpl)))
        n_female_by_g <- pmax(0L, pmin(n_total_by_g, as.integer(round(n_total_by_g * phi))))
        n <- n + sum(n_female_by_g)
      } else {
        n <- n + p_counts[[i]] * rate
      }
    }
    available <- state$get_size_of('NonExistent')
    if (n > available) {
      stop(paste0(
        'Not enough mosquitoes (short by ',
        n - available,
        '). Please raise parameters$mosquito_limit. ',
        'If you have used parameterise_mosquito_equilibrium,',
        'your seasonality parameters lead to more mosquitoes than expected.'
      ))
    }
    non_existent <- state$get_index_of('NonExistent')
    latest <- 1
    for (i in seq_along(species_names)) {
      solver_states <- solver_rows[[i]]
      if (!is.null(models[[i]]$cube) &&
          uses_genotype_resolved_aquatic_states(solver_states, models[[i]]$cube)) {
        cube_info <- cube_genotype_info(models[[i]]$cube)
        p_by_g <- aquatic_stage_values_by_genotype(solver_states, models[[i]]$cube, "P")
        phi <- cube_phi_vector(models[[i]]$cube, cube_info$G)
        n_total_by_g <- pmax(0L, as.integer(round(as.numeric(p_by_g) / parameters$dpl)))
        n_female_by_g <- pmax(0L, pmin(n_total_by_g, as.integer(round(n_total_by_g * phi))))
        n_male_by_g <- n_total_by_g - n_female_by_g

        for (g in seq_len(cube_info$G)) {
          n_female <- n_female_by_g[[g]]
          if (n_female > 0L) {
            hatched <- bitset_at(non_existent, seq.int(latest, latest + n_female - 1L))
            geno_id$queue_update(g, hatched)
            state$queue_update('Sm', hatched)
            species$queue_update(species_names[[i]], hatched)
            latest <- latest + n_female
          }
        }
        models[[i]]$genotype_state$male_counts <- models[[i]]$genotype_state$male_counts + n_male_by_g
      } else {
        to_hatch <- p_counts[[i]] * rate
        hatched <- bitset_at(non_existent, seq(latest, latest + to_hatch))
        n_hatched <- hatched$size()

        if (!is.null(models[[i]]$cube) && n_hatched > 0) {
        cube_info <- cube_genotype_info(models[[i]]$cube)
        female_counts <- adult_female_genotype_counts_by_species(
          variables,
          species_names[[i]],
          cube_info$G
        )
        names(female_counts) <- cube_info$genotypesID
        male_counts <- models[[i]]$genotype_state$male_counts
        if (is.null(names(male_counts))) {
          names(male_counts) <- cube_info$genotypesID
        }
        pgv <- calc_pg_V_from_cube(models[[i]]$cube, female_counts, male_counts)
        names(pgv$p) <- cube_info$genotypesID
        female_assigned <- sample_genotype_counts(n_hatched, pgv$p)
        names(female_assigned) <- cube_info$genotypesID
        male_assigned <- sample_genotype_counts(n_hatched, pgv$p)
        names(male_assigned) <- cube_info$genotypesID
        geno_id$queue_update(rep.int(seq_along(female_assigned), female_assigned), hatched)
        models[[i]]$genotype_state$male_counts <- male_counts + male_assigned
        if (genotype_debug_enabled(parameters, timestep)) {
          genotype_debug_log(
            parameters,
            timestep,
            "EMERGE_PG",
            species_names[[i]],
            sprintf(
              "p={%s} n_new(total=%d,female=%d,male=%d) female_assigned={%s} male_assigned={%s}",
              genotype_debug_fmt_counts(pgv$p),
              as.integer(2 * n_hatched),
              as.integer(n_hatched),
              as.integer(n_hatched),
              genotype_debug_fmt_counts(female_assigned),
              genotype_debug_fmt_counts(male_assigned)
            )
          )
          genotype_debug_log_counts(
            parameters,
            timestep,
            "AFTER_EMERGE_VISIBLE",
            species_names[[i]],
            genotype_debug_species_counts(variables, models, parameters, i),
            extra = "female_assigned queued (not yet visible this timestep)"
          )
        }
        } else if (n_hatched > 0) {
          geno_id$queue_update(1L, hatched)
        }

        state$queue_update('Sm', hatched)
        species$queue_update(species_names[[i]], hatched)
        latest <- latest + to_hatch + 1
      }
    }
  }
}
