mosquito_standalone_fields <- c(
  "species", "species_proportions", "total_M",
  "del", "dl", "dpl", "me", "ml", "mup", "mum",
  "beta", "blood_meal_rates", "Q0",
  "g0", "g", "h", "gamma", "model_seasonality", "rainfall_floor",
  "carrying_capacity", "carrying_capacity_timesteps",
  "carrying_capacity_values", "carrying_capacity_scalers",
  "native_mosquito_nE", "native_mosquito_nL", "native_mosquito_nP",
  "native_mosquito_nEIP", "native_mosquito_nu", "mosquito_tau_step",
  "cube", "move_probs", "move_rates", "mosquito_move_probs", "mosquito_move_rates",
  "r_tol", "a_tol", "ode_max_steps",
  "bednets", "spraying"
)

mosquito_default_extras <- function() {
  list(
    carrying_capacity_scalers = NULL,
    native_mosquito_nE = NULL,
    native_mosquito_nL = NULL,
    native_mosquito_nP = NULL,
    native_mosquito_nEIP = NULL,
    native_mosquito_nu = NULL,
    mosquito_tau_step = NULL
  )
}

mosquito_scalar <- function(x, name, positive = FALSE, nonnegative = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("`%s` must be a finite numeric scalar.", name), call. = FALSE)
  }
  if (positive && x <= 0) {
    stop(sprintf("`%s` must be > 0.", name), call. = FALSE)
  }
  if (nonnegative && x < 0) {
    stop(sprintf("`%s` must be >= 0.", name), call. = FALSE)
  }
  as.numeric(x)
}

mosquito_validate_parameters <- function(parameters) {
  for (name in c("del", "dl", "dpl", "total_M", "rainfall_floor", "r_tol", "a_tol")) {
    parameters[[name]] <- mosquito_scalar(parameters[[name]], name, positive = TRUE)
  }
  for (name in c("me", "ml", "mup", "gamma")) {
    parameters[[name]] <- mosquito_scalar(parameters[[name]], name, nonnegative = TRUE)
  }

  if (isTRUE(parameters$bednets) || isTRUE(parameters$spraying)) {
    stop(
      "Standalone mosquito simulations do not support human-derived vector-control schedules; ",
      "provide already-adjusted mosquito mortality/fecundity values instead.",
      call. = FALSE
    )
  }

  parameters$species <- as.character(parameters$species)
  n_species <- length(parameters$species)
  if (n_species < 1L || anyNA(parameters$species) || any(parameters$species == "")) {
    stop("`species` must contain at least one non-empty species name.", call. = FALSE)
  }
  if (anyDuplicated(parameters$species)) {
    stop("`species` names must be unique.", call. = FALSE)
  }

  parameters$species_proportions <- as.numeric(parameters$species_proportions)
  if (length(parameters$species_proportions) != n_species ||
      anyNA(parameters$species_proportions) ||
      any(!is.finite(parameters$species_proportions)) ||
      any(parameters$species_proportions < 0) ||
      sum(parameters$species_proportions) <= 0) {
    stop("`species_proportions` must be finite non-negative values matching `species`.", call. = FALSE)
  }
  parameters$species_proportions <- parameters$species_proportions / sum(parameters$species_proportions)

  for (name in c("mum", "blood_meal_rates", "Q0", "beta")) {
    value <- as.numeric(parameters[[name]])
    if (!(length(value) %in% c(1L, n_species)) ||
        anyNA(value) || any(!is.finite(value)) || any(value < 0)) {
      stop(sprintf("`%s` must be a non-negative scalar or one value per species.", name), call. = FALSE)
    }
    parameters[[name]] <- value
  }
  if (any(parameters$mum <= 0) || any(parameters$beta <= 0)) {
    stop("`mum` and `beta` must be > 0 for standalone mosquito dynamics.", call. = FALSE)
  }

  cfg <- native_mosquito_stage_config(parameters)
  for (name in c("nE", "nL", "nP")) {
    if (!is.finite(cfg[[name]]) || cfg[[name]] < 1L) {
      stop("Native mosquito aquatic stage counts must be positive integers.", call. = FALSE)
    }
  }
  if (!is.finite(cfg$nu) || cfg$nu < 0) {
    stop("`native_mosquito_nu` must be a finite non-negative scalar.", call. = FALSE)
  }
  if (!is.finite(cfg$dt_stoch) || cfg$dt_stoch <= 0 || cfg$dt_stoch > 1) {
    stop("`mosquito_tau_step` must be in (0, 1].", call. = FALSE)
  }

  cube_genotype_info(parameters$cube)
  class(parameters) <- unique(c("mosquito_parameters", setdiff(class(parameters), "mosquito_parameters")))
  parameters
}

#' Get standalone mosquito lifecycle parameters
#'
#' Creates a mosquito-only parameter object for `run_mosquito_simulation()`.
#' Defaults are copied from [get_parameters()] for mosquito biology, genetics,
#' movement, carrying capacity, and solver controls, but human and parasite
#' dynamics are not used by the standalone runner.
#'
#' @param overrides named list of parameter values overriding defaults.
#' @param cube optional MGDrivE-style inheritance cube.
#'
#' @return A validated object of class `mosquito_parameters`.
#' @export
get_mosquito_parameters <- function(overrides = list(), cube = NULL) {
  if (!is.list(overrides)) {
    stop("`overrides` must be a named list.", call. = FALSE)
  }
  if (length(overrides) > 0L && (is.null(names(overrides)) || any(names(overrides) == ""))) {
    stop("`overrides` must be a named list.", call. = FALSE)
  }
  if (!is.null(cube)) {
    overrides$cube <- cube
  }
  base_full <- get_parameters()
  accepted <- union(names(base_full), mosquito_standalone_fields)
  unknown <- setdiff(names(overrides), accepted)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown standalone mosquito parameter override(s): %s.", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  full_overrides <- overrides[intersect(names(overrides), names(base_full))]
  base <- get_parameters(full_overrides)
  keep <- intersect(mosquito_standalone_fields, names(base))
  parameters <- base[keep]
  parameters <- c(parameters, mosquito_default_extras()[setdiff(names(mosquito_default_extras()), names(parameters))])
  mosquito_overrides <- overrides[intersect(names(overrides), mosquito_standalone_fields)]
  for (name in names(mosquito_overrides)) {
    parameters[[name]] <- mosquito_overrides[[name]]
  }
  for (name in setdiff(mosquito_standalone_fields, names(parameters))) {
    parameters[[name]] <- NULL
  }
  mosquito_validate_parameters(parameters)
}

#' Extract standalone mosquito parameters
#'
#' Converts a full malaria simulation parameter list into the mosquito-only
#' parameter object consumed by `run_mosquito_simulation()`.
#'
#' @param parameters A `get_parameters()` result, a `mosquito_parameters`
#' object, or a named list of mosquito fields.
#'
#' @return A validated object of class `mosquito_parameters`.
#' @export
as_mosquito_parameters <- function(parameters) {
  if (inherits(parameters, "mosquito_parameters")) {
    return(mosquito_validate_parameters(parameters))
  }
  if (!is.list(parameters)) {
    stop("`parameters` must be a parameter list.", call. = FALSE)
  }
  keep <- intersect(mosquito_standalone_fields, names(parameters))
  out <- parameters[keep]
  defaults <- get_mosquito_parameters()
  for (name in setdiff(names(defaults), names(out))) {
    out[[name]] <- defaults[[name]]
  }
  mosquito_validate_parameters(out)
}

mosquito_normalize_network <- function(parameters) {
  if (inherits(parameters, "mosquito_parameters") || is.null(parameters[[1L]])) {
    parameters <- list(node1 = as_mosquito_parameters(parameters))
  } else {
    parameters <- lapply(parameters, as_mosquito_parameters)
    if (is.null(names(parameters)) || any(names(parameters) == "")) {
      names(parameters) <- paste0("node", seq_along(parameters))
    }
  }
  native_assert_identical_configs(parameters, native_mosquito_stage_config(parameters[[1L]]))
  ref_species <- parameters[[1L]]$species
  ref_cube <- parameters[[1L]]$cube
  for (i in seq_along(parameters)) {
    if (!identical(parameters[[i]]$species, ref_species)) {
      stop("All standalone mosquito nodes must use identical species ordering.", call. = FALSE)
    }
    if (!isTRUE(all.equal(parameters[[i]]$cube, ref_cube, check.attributes = FALSE))) {
      stop("All standalone mosquito nodes must use the same inheritance cube.", call. = FALSE)
    }
  }
  parameters
}

build_ento_indices <- function(G, nNodes, nE, nL, nP) {
  idx <- 1L
  take_block <- function(len, dims) {
    block <- array(seq.int(idx, length.out = len), dim = dims)
    idx <<- idx + len
    block
  }
  egg_ix <- take_block(nE * G * nNodes, c(nE, G, nNodes))
  larv_ix <- take_block(nL * G * nNodes, c(nL, G, nNodes))
  pup_ix <- take_block(nP * G * nNodes, c(nP, G, nNodes))
  male_ix <- take_block(G * nNodes, c(G, nNodes))
  unm_ix <- take_block(G * nNodes, c(G, nNodes))
  fem_ix <- take_block(G * G * nNodes, c(G * G, nNodes))
  list(
    egg_ix = egg_ix, larv_ix = larv_ix, pup_ix = pup_ix,
    male_ix = male_ix, unm_ix = unm_ix, fem_ix = fem_ix,
    total_state_len = idx - 1L
  )
}

mosquito_species_value <- function(parameters, name, species_i) {
  value <- parameters[[name]]
  if (length(value) == 1L) {
    return(as.numeric(value[[1L]]))
  }
  if (!is.null(names(value))) {
    return(as.numeric(value[[parameters$species[[species_i]]]]))
  }
  as.numeric(value[[species_i]])
}

mosquito_exact_equilibrium_node <- function(parameters, species_i, cube_info, timestep = 0L) {
  cfg <- native_mosquito_stage_config(parameters)
  G <- cube_info$G
  wt <- cube_info$wild_type_index
  gids <- cube_info$genotypesID
  wt_pair <- wt + (wt - 1L) * G
  m <- parameters$total_M * parameters$species_proportions[[species_i]]

  B_mat <- native_build_birth_matrix(parameters$cube, cube_info)
  birth_weights <- as.numeric(B_mat[wt_pair, ])
  expected_birth <- rep(0, G)
  expected_birth[[wt]] <- 1
  if (!isTRUE(all.equal(birth_weights, expected_birth, tolerance = 1e-10))) {
    stop(
      "Default standalone mosquito equilibrium requires wild-type self-cross to produce only wild-type offspring; ",
      "provide `initial_state` for this cube.",
      call. = FALSE
    )
  }

  omega <- cube_omega_vector(parameters$cube, G, gids)
  omega_inv <- ifelse(omega == 0, 1e3, 1 / omega)
  phi <- native_align_cube_vec(if (is.null(parameters$cube)) NULL else parameters$cube$phi, gids, 0.5)
  xiF <- native_align_cube_vec(if (is.null(parameters$cube)) NULL else parameters$cube$xiF, gids, 1)
  xiM <- native_align_cube_vec(if (is.null(parameters$cube)) NULL else parameters$cube$xiM, gids, 1)
  rE <- cfg$nE / parameters$del
  rL <- cfg$nL / parameters$dl
  rP <- cfg$nP / parameters$dpl
  muF <- mosquito_species_value(parameters, "mum", species_i) * omega_inv[[wt]]
  muM <- mosquito_species_value(parameters, "mum", species_i) * omega_inv[[wt]]

  if (m <= 0) {
    return(list(k0 = 0, egg = rep(0, cfg$nE), larv = rep(0, cfg$nL),
                pup = rep(0, cfg$nP), male = 0, female = 0))
  }
  if (any(c(muF, muM, phi[[wt]], xiF[[wt]], xiM[[wt]], rE, rL, rP) <= 0) ||
      any(!is.finite(c(muF, muM, phi[[wt]], xiF[[wt]], xiM[[wt]], rE, rL, rP)))) {
    stop("Default standalone mosquito equilibrium has invalid wild-type traits; provide `initial_state`.", call. = FALSE)
  }

  beta_eff <- mosquito_species_value(parameters, "beta", species_i)
  aE <- rE / (parameters$me + rE)
  egg_1 <- beta_eff * m / (parameters$me + rE)
  egg <- egg_1 * aE^(0:(cfg$nE - 1L))
  egg_last <- egg[[cfg$nE]]

  aP <- rP / (parameters$mup + rP)
  pup_last <- muF * m / (phi[[wt]] * xiF[[wt]] * rP)
  pup_1 <- pup_last / aP^(cfg$nP - 1L)
  pup <- pup_1 * aP^(0:(cfg$nP - 1L))

  larv_last <- (parameters$mup + rP) * pup_1 / rL
  dd_rate <- (rE * egg_last * rL^(cfg$nL - 1L) / larv_last)^(1 / cfg$nL) - rL
  if (!is.finite(dd_rate) || dd_rate <= parameters$ml) {
    stop("Default standalone mosquito equilibrium could not derive positive density dependence; provide `initial_state`.", call. = FALSE)
  }
  aL <- rL / (dd_rate + rL)
  larv_1 <- rE * egg_last / (dd_rate + rL)
  larv <- larv_1 * aL^(0:(cfg$nL - 1L))
  K_eff <- sum(larv) / (dd_rate / parameters$ml - 1)
  scale_t0 <- native_carrying_capacity_scale(parameters, species_i, timestep)
  if (!is.finite(K_eff) || K_eff <= 0 || !is.finite(scale_t0) || scale_t0 <= 0) {
    stop("Default standalone mosquito equilibrium could not derive carrying capacity; provide `initial_state`.", call. = FALSE)
  }

  list(
    k0 = K_eff / scale_t0,
    egg = egg,
    larv = larv,
    pup = pup,
    male = ((1 - phi[[wt]]) * xiM[[wt]] * rP * pup_last) / muM,
    female = m
  )
}

mosquito_build_carrying_lookup <- function(parameters, species_i, cube_info) {
  lapply(parameters, function(p) {
    eq <- mosquito_exact_equilibrium_node(p, species_i, cube_info, timestep = 0L)
    list(
      k0 = eq$k0,
      model_seasonality = p$model_seasonality,
      g0 = p$g0,
      g = p$g,
      h = p$h,
      R_bar = calculate_R_bar(p),
      rainfall_floor = p$rainfall_floor,
      carrying_capacity_timesteps = if (isTRUE(p$carrying_capacity) &&
        !is.null(p$carrying_capacity_timesteps)) as.integer(p$carrying_capacity_timesteps) else integer(0),
      carrying_capacity_scalers = if (isTRUE(p$carrying_capacity) &&
        !is.null(p$carrying_capacity_scalers)) as.numeric(p$carrying_capacity_scalers[, species_i]) else numeric(0)
    )
  })
}

mosquito_carrying_at <- function(lookup, timestep) {
  vapply(lookup, function(node) {
    k0 <- node$k0
    if (length(node$carrying_capacity_timesteps) > 0L) {
      changes <- which(node$carrying_capacity_timesteps <= timestep)
      if (length(changes) > 0L) {
        k0 <- k0 * node$carrying_capacity_scalers[[max(changes)]]
      }
    }
    carrying_capacity(timestep, node$model_seasonality, node$g0, node$g, node$h,
                      k0, node$R_bar, node$rainfall_floor)
  }, numeric(1))
}

mosquito_initial_state <- function(parameters, species_i, cube_info, index, sampler) {
  cfg <- native_mosquito_stage_config(parameters[[1L]])
  state <- numeric(index$total_state_len)
  G <- cube_info$G
  wt <- cube_info$wild_type_index
  for (node in seq_along(parameters)) {
    eq <- mosquito_exact_equilibrium_node(parameters[[node]], species_i, cube_info, 0L)
    state[index$egg_ix[, wt, node]] <- eq$egg
    state[index$larv_ix[, wt, node]] <- eq$larv
    state[index$pup_ix[, wt, node]] <- eq$pup
    state[index$male_ix[wt, node]] <- eq$male
    state[index$fem_ix[wt + (wt - 1L) * G, node]] <- eq$female
  }
  if (sampler == "tau") {
    state <- native_integerise_counts_preserve_total(state)
  }
  state
}

mosquito_config_signature <- function(parameters, sampler, cube_info, cfg, index) {
  list(
    sampler = sampler,
    nodes = names(parameters),
    species = parameters[[1L]]$species,
    genotypes = cube_info$genotypesID,
    stage_config = cfg[c("nE", "nL", "nP", "nu", "dt_stoch")],
    state_len = index$total_state_len
  )
}

mosquito_build_species_backend <- function(parameters, species_i, sampler, timesteps) {
  cfg <- native_mosquito_stage_config(parameters[[1L]])
  cube <- parameters[[1L]]$cube
  cube_info <- cube_genotype_info(cube)
  index <- build_ento_indices(cube_info$G, length(parameters), cfg$nE, cfg$nL, cfg$nP)
  move <- native_collect_mosquito_movement(parameters)
  gids <- cube_info$genotypesID
  omega <- cube_omega_vector(cube, cube_info$G, gids)
  phi <- native_align_cube_vec(if (is.null(cube)) NULL else cube$phi, gids, 0.5)
  xiF <- native_align_cube_vec(if (is.null(cube)) NULL else cube$xiF, gids, 1)
  xiM <- native_align_cube_vec(if (is.null(cube)) NULL else cube$xiM, gids, 1)
  eta <- native_align_cube_mat(if (is.null(cube)) NULL else cube$eta, gids, 1)
  B_mat <- native_build_birth_matrix(cube, cube_info)
  K_lookup <- mosquito_build_carrying_lookup(parameters, species_i, cube_info)
  init_K <- mosquito_carrying_at(K_lookup, 0L)
  beta_vec <- vapply(parameters, function(p) mosquito_species_value(p, "beta", species_i), numeric(1))
  mu_vec <- vapply(parameters, function(p) mosquito_species_value(p, "mum", species_i), numeric(1))
  muE <- vapply(parameters, function(p) p$me, numeric(1))
  muL <- vapply(parameters, function(p) p$ml, numeric(1))
  muP <- vapply(parameters, function(p) p$mup, numeric(1))
  gamma_dd <- rep(0, length(parameters))

  engine_ptr <- ento_engine_create(
    nNodes = length(parameters), nG = cube_info$G,
    nE = cfg$nE, nL = cfg$nL, nP = cfg$nP,
    egg_ix = index$egg_ix, larv_ix = index$larv_ix, pup_ix = index$pup_ix,
    male_ix = index$male_ix, unm_ix = index$unm_ix, fem_ix = index$fem_ix,
    rE = cfg$nE / parameters[[1L]]$del,
    rL = cfg$nL / parameters[[1L]]$dl,
    rP = cfg$nP / parameters[[1L]]$dpl,
    muE = muE, muL = muL, muP = muP, muM = mu_vec, muF = mu_vec,
    log_dd = TRUE, K = init_K, gamma_dd = gamma_dd,
    beta = beta_vec[[1L]], nu = cfg$nu,
    omega = omega, phi_cube = phi, xiF = xiF, xiM = xiM,
    eta = eta, B_mat = B_mat,
    has_move = move$has_move, move_probs_dense = move$move_probs,
    move_rates = move$move_rates,
    muM_node_base = rep(1, length(parameters)),
    muF_node_base = rep(1, length(parameters)),
    tol = 1e-8, nState = index$total_state_len
  )
  ento_engine_set_runtime(engine_ptr, beta_vec, mu_vec, mu_vec, init_K, gamma_dd)

  state <- mosquito_initial_state(parameters, species_i, cube_info, index, sampler)
  solver_ptr <- if (sampler == "ode") {
    create_ento_solver(engine_ptr, state,
                       parameters[[1L]]$r_tol, parameters[[1L]]$a_tol,
                       as.integer(parameters[[1L]]$ode_max_steps))
  } else {
    NULL
  }

  list(
    engine_ptr = engine_ptr, solver_ptr = solver_ptr, state = state, t = 0,
    species_i = species_i, species_name = parameters[[1L]]$species[[species_i]],
    index = index, cube_info = cube_info, cfg = cfg, move = move,
    beta_vec = beta_vec, mu_vec = mu_vec, muE = muE, muL = muL, muP = muP,
    gamma_dd = gamma_dd, omega = omega, phi = phi, xiF = xiF, xiM = xiM,
    eta = eta, B_mat = B_mat, K_lookup = K_lookup,
    signature = mosquito_config_signature(parameters, sampler, cube_info, cfg, index)
  )
}

mosquito_state_vector <- function(backend) {
  if (is.null(backend$solver_ptr)) {
    backend$state
  } else {
    as.numeric(solver_get_states(backend$solver_ptr))
  }
}

mosquito_set_state_vector <- function(backend, state, t) {
  backend$state <- as.numeric(state)
  backend$t <- as.numeric(t)
  if (!is.null(backend$solver_ptr)) {
    solver_set_states(backend$solver_ptr, as.numeric(t), as.numeric(state))
  }
  backend
}

mosquito_state_counts_one <- function(state, backend, node_names, timestep) {
  G <- backend$cube_info$G
  gids <- backend$cube_info$genotypesID
  rows <- list()
  idx <- 1L
  add <- function(node, species, stage, genotype, count, substage = NA_integer_,
                  mating_status = NA_character_, mate_genotype = NA_character_) {
    rows[[idx]] <<- data.frame(
      timestep = timestep, node = node, species = species, stage = stage,
      genotype = genotype, substage = substage, mating_status = mating_status,
      mate_genotype = mate_genotype, count = as.numeric(count),
      stringsAsFactors = FALSE
    )
    idx <<- idx + 1L
  }
  for (node_i in seq_along(node_names)) {
    for (g in seq_len(G)) {
      for (e in seq_len(backend$cfg$nE)) {
        add(node_names[[node_i]], backend$species_name, "egg", gids[[g]],
            state[backend$index$egg_ix[e, g, node_i]], e)
      }
      for (l in seq_len(backend$cfg$nL)) {
        add(node_names[[node_i]], backend$species_name, "larva", gids[[g]],
            state[backend$index$larv_ix[l, g, node_i]], l)
      }
      for (p in seq_len(backend$cfg$nP)) {
        add(node_names[[node_i]], backend$species_name, "pupa", gids[[g]],
            state[backend$index$pup_ix[p, g, node_i]], p)
      }
      add(node_names[[node_i]], backend$species_name, "adult_male", gids[[g]],
          state[backend$index$male_ix[g, node_i]])
      add(node_names[[node_i]], backend$species_name, "adult_female", gids[[g]],
          state[backend$index$unm_ix[g, node_i]], mating_status = "unmated")
      for (gm in seq_len(G)) {
        add(node_names[[node_i]], backend$species_name, "adult_female", gids[[g]],
            state[backend$index$fem_ix[g + (gm - 1L) * G, node_i]],
            mating_status = "mated", mate_genotype = gids[[gm]])
      }
    }
  }
  do.call(rbind, rows)
}

mosquito_aggregate_counts <- function(state_counts) {
  stats::aggregate(count ~ timestep + node + species + stage + genotype,
                   data = state_counts, FUN = sum)
}

mosquito_normalize_releases <- function(releases, parameters, cube_info) {
  if (is.null(releases) || NROW(releases) == 0L) {
    return(data.frame(
      timestep = integer(0), node = character(0), species = character(0),
      stage = character(0), genotype = character(0), count = numeric(0),
      substage = integer(0), mate_genotype = character(0),
      stringsAsFactors = FALSE
    ))
  }
  releases <- as.data.frame(releases, stringsAsFactors = FALSE)
  required <- c("timestep", "node", "species", "stage", "genotype", "count")
  missing <- setdiff(required, names(releases))
  if (length(missing) > 0L) {
    stop(sprintf("`releases` is missing required columns: %s.", paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!("substage" %in% names(releases))) releases$substage <- NA_integer_
  if (!("mate_genotype" %in% names(releases))) releases$mate_genotype <- NA_character_
  releases$timestep <- as.integer(releases$timestep)
  releases$count <- as.numeric(releases$count)
  if (anyNA(releases$timestep) || any(releases$timestep < 0L) ||
      anyNA(releases$count) || any(!is.finite(releases$count)) || any(releases$count < 0)) {
    stop("Release timesteps and counts must be finite non-negative values.", call. = FALSE)
  }
  releases$node <- as.character(releases$node)
  releases$species <- as.character(releases$species)
  releases$stage <- as.character(releases$stage)
  releases$genotype <- as.character(releases$genotype)
  releases$mate_genotype <- as.character(releases$mate_genotype)
  releases$mate_genotype[is.na(releases$mate_genotype) | releases$mate_genotype == "NA"] <- NA_character_
  valid_stages <- c("egg", "larva", "pupa", "adult_male", "adult_female")
  if (any(!(releases$stage %in% valid_stages))) {
    stop("Release `stage` must be one of egg, larva, pupa, adult_male, adult_female.", call. = FALSE)
  }
  if (any(!(releases$node %in% names(parameters)))) {
    stop("Release `node` values must match standalone mosquito node names.", call. = FALSE)
  }
  if (any(!(releases$species %in% parameters[[1L]]$species))) {
    stop("Release `species` values must match configured species names.", call. = FALSE)
  }
  if (any(!(releases$genotype %in% cube_info$genotypesID))) {
    stop("Release `genotype` values must match cube genotypes.", call. = FALSE)
  }
  bad_mate <- !is.na(releases$mate_genotype) & !(releases$mate_genotype %in% cube_info$genotypesID)
  if (any(bad_mate)) {
    stop("Release `mate_genotype` values must match cube genotypes when provided.", call. = FALSE)
  }
  releases
}

mosquito_apply_release <- function(state, backend, release, node_i, stochastic) {
  G <- backend$cube_info$G
  genotype_i <- match(release$genotype, backend$cube_info$genotypesID)
  count <- as.numeric(release$count)
  if (count <= 0) return(state)
  if (stochastic) count <- round(count)
  stage <- release$stage
  if (stage == "egg") {
    substage <- if (is.na(release$substage)) 1L else as.integer(release$substage)
    if (is.na(substage) || substage < 1L || substage > backend$cfg$nE) {
      stop("Egg releases require a valid `substage`.", call. = FALSE)
    }
    state[backend$index$egg_ix[substage, genotype_i, node_i]] <- state[backend$index$egg_ix[substage, genotype_i, node_i]] + count
  } else if (stage == "larva") {
    substage <- if (is.na(release$substage)) 1L else as.integer(release$substage)
    if (is.na(substage) || substage < 1L || substage > backend$cfg$nL) {
      stop("Larva releases require a valid `substage`.", call. = FALSE)
    }
    state[backend$index$larv_ix[substage, genotype_i, node_i]] <- state[backend$index$larv_ix[substage, genotype_i, node_i]] + count
  } else if (stage == "pupa") {
    substage <- if (is.na(release$substage)) 1L else as.integer(release$substage)
    if (is.na(substage) || substage < 1L || substage > backend$cfg$nP) {
      stop("Pupa releases require a valid `substage`.", call. = FALSE)
    }
    state[backend$index$pup_ix[substage, genotype_i, node_i]] <- state[backend$index$pup_ix[substage, genotype_i, node_i]] + count
  } else if (stage == "adult_male") {
    state[backend$index$male_ix[genotype_i, node_i]] <- state[backend$index$male_ix[genotype_i, node_i]] + count
  } else {
    mate <- release$mate_genotype
    if (!is.na(mate)) {
      mate_i <- match(mate, backend$cube_info$genotypesID)
      state[backend$index$fem_ix[genotype_i + (mate_i - 1L) * G, node_i]] <-
        state[backend$index$fem_ix[genotype_i + (mate_i - 1L) * G, node_i]] + count
    } else {
      male_counts <- state[backend$index$male_ix[, node_i]]
      weights <- as.numeric(backend$eta[genotype_i, ]) * male_counts
      if (sum(weights) > 0) {
        probs <- weights / sum(weights)
        alloc <- if (stochastic) sample_genotype_counts(count, probs) else count * probs
        for (gm in seq_len(G)) {
          state[backend$index$fem_ix[genotype_i + (gm - 1L) * G, node_i]] <-
            state[backend$index$fem_ix[genotype_i + (gm - 1L) * G, node_i]] + alloc[[gm]]
        }
      } else {
        state[backend$index$unm_ix[genotype_i, node_i]] <- state[backend$index$unm_ix[genotype_i, node_i]] + count
      }
    }
  }
  state
}

mosquito_apply_releases_at <- function(backends, releases, parameters, timestep, sampler) {
  if (NROW(releases) == 0L) return(backends)
  due <- releases[releases$timestep == timestep, , drop = FALSE]
  if (NROW(due) == 0L) return(backends)
  for (species_i in seq_along(backends)) {
    backend <- backends[[species_i]]
    state <- mosquito_state_vector(backend)
    species_due <- due[due$species == backend$species_name, , drop = FALSE]
    for (i in seq_len(NROW(species_due))) {
      node_i <- match(species_due$node[[i]], names(parameters))
      state <- mosquito_apply_release(state, backend, species_due[i, , drop = FALSE], node_i, sampler == "tau")
    }
    backend <- mosquito_set_state_vector(backend, state, timestep)
    backends[[species_i]] <- backend
  }
  backends
}

mosquito_step_backend <- function(backend, timestep, sampler, parameters) {
  K <- mosquito_carrying_at(backend$K_lookup, timestep - 1L)
  ento_engine_set_runtime(backend$engine_ptr, backend$beta_vec, backend$mu_vec, backend$mu_vec, K, backend$gamma_dd)
  if (sampler == "ode") {
    solver_step(backend$solver_ptr)
    backend$t <- timestep
    backend$state <- as.numeric(solver_get_states(backend$solver_ptr))
  } else {
    backend$state <- stoch_ento_step_native(
      state_r = backend$state, t = timestep - 1,
      dt_stoch = backend$cfg$dt_stoch, dt_out = 1,
      nNodes = length(parameters), nG = backend$cube_info$G,
      nE = backend$cfg$nE, nL = backend$cfg$nL, nP = backend$cfg$nP,
      egg_ix_r = as.integer(backend$index$egg_ix),
      larv_ix_r = as.integer(backend$index$larv_ix),
      pup_ix_r = as.integer(backend$index$pup_ix),
      male_ix_r = as.integer(backend$index$male_ix),
      unm_ix_r = as.integer(backend$index$unm_ix),
      fem_ix_r = as.integer(backend$index$fem_ix),
      rE = backend$cfg$nE / parameters[[1L]]$del,
      rL = backend$cfg$nL / parameters[[1L]]$dl,
      rP = backend$cfg$nP / parameters[[1L]]$dpl,
      muE_r = backend$muE, muL_r = backend$muL, muP_r = backend$muP,
      muM_r = backend$mu_vec, muF_r = backend$mu_vec,
      log_dd = TRUE, K_r = K, gamma_dd_r = backend$gamma_dd,
      beta_vec_r = backend$beta_vec, nu = backend$cfg$nu,
      omega_inv_r = ifelse(backend$omega == 0, 1e3, 1 / backend$omega),
      phi_r = backend$phi, xiF_r = backend$xiF, xiM_r = backend$xiM,
      eta_r = backend$eta, B_mat_r = backend$B_mat,
      has_move = backend$move$has_move,
      move_probs_r = backend$move$move_probs,
      move_rates_r = backend$move$move_rates,
      muM_node_base_r = rep(1, length(parameters)),
      muF_node_base_r = rep(1, length(parameters)),
      tol = 1e-8
    )
    backend$t <- timestep
  }
  backend
}

mosquito_checkpoint <- function(backends, sampler, timestep) {
  list(
    type = "standalone_mosquito_state",
    sampler = sampler,
    timestep = as.integer(timestep),
    species = lapply(backends, function(b) {
      list(
        species = b$species_name,
        state = mosquito_state_vector(b),
        time = b$t,
        signature = b$signature
      )
    }),
    rng_state = if (sampler == "tau" && exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
  )
}

mosquito_restore_checkpoint <- function(backends, checkpoint, sampler) {
  if (!is.list(checkpoint) || !identical(checkpoint$type, "standalone_mosquito_state")) {
    return(NULL)
  }
  if (!identical(checkpoint$sampler, sampler)) {
    stop("Standalone mosquito checkpoint sampler does not match requested sampler.", call. = FALSE)
  }
  if (length(checkpoint$species) != length(backends)) {
    stop("Standalone mosquito checkpoint species layout does not match parameters.", call. = FALSE)
  }
  for (i in seq_along(backends)) {
    if (!identical(checkpoint$species[[i]]$signature, backends[[i]]$signature)) {
      stop("Standalone mosquito checkpoint configuration does not match parameters.", call. = FALSE)
    }
    backends[[i]] <- mosquito_set_state_vector(
      backends[[i]],
      checkpoint$species[[i]]$state,
      checkpoint$species[[i]]$time
    )
  }
  if (sampler == "tau" && !is.null(checkpoint$rng_state)) {
    assign(".Random.seed", checkpoint$rng_state, envir = .GlobalEnv)
  }
  list(backends = backends, timestep = as.integer(checkpoint$timestep))
}

mosquito_normalize_initial_state_table <- function(initial_state, parameters, cube_info) {
  initial_state <- as.data.frame(initial_state, stringsAsFactors = FALSE)
  required <- c("node", "species", "stage", "genotype", "count")
  missing <- setdiff(required, names(initial_state))
  if (length(missing) > 0L) {
    stop(sprintf("`initial_state` is missing required columns: %s.", paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!("substage" %in% names(initial_state))) initial_state$substage <- rep(NA_integer_, NROW(initial_state))
  if (!("mating_status" %in% names(initial_state))) initial_state$mating_status <- rep(NA_character_, NROW(initial_state))
  if (!("mate_genotype" %in% names(initial_state))) initial_state$mate_genotype <- rep(NA_character_, NROW(initial_state))
  initial_state$node <- as.character(initial_state$node)
  initial_state$species <- as.character(initial_state$species)
  initial_state$stage <- as.character(initial_state$stage)
  initial_state$genotype <- as.character(initial_state$genotype)
  initial_state$mating_status <- as.character(initial_state$mating_status)
  initial_state$mate_genotype <- as.character(initial_state$mate_genotype)
  initial_state$mating_status[is.na(initial_state$mating_status) | initial_state$mating_status == "NA"] <- NA_character_
  initial_state$mate_genotype[is.na(initial_state$mate_genotype) | initial_state$mate_genotype == "NA"] <- NA_character_
  initial_state$count <- as.numeric(initial_state$count)
  if (anyNA(initial_state$count) || any(!is.finite(initial_state$count)) || any(initial_state$count < 0)) {
    stop("`initial_state$count` must contain finite non-negative values.", call. = FALSE)
  }
  valid_stages <- c("egg", "larva", "pupa", "adult_male", "adult_female")
  if (any(!(initial_state$stage %in% valid_stages))) {
    stop("`initial_state$stage` contains unsupported stages.", call. = FALSE)
  }
  if (any(!(initial_state$node %in% names(parameters)))) {
    stop("`initial_state$node` values must match standalone mosquito node names.", call. = FALSE)
  }
  if (any(!(initial_state$species %in% parameters[[1L]]$species))) {
    stop("`initial_state$species` values must match configured species names.", call. = FALSE)
  }
  if (any(!(initial_state$genotype %in% cube_info$genotypesID))) {
    stop("`initial_state$genotype` values must match cube genotypes.", call. = FALSE)
  }
  bad_mate <- !is.na(initial_state$mate_genotype) & !(initial_state$mate_genotype %in% cube_info$genotypesID)
  if (any(bad_mate)) {
    stop("`initial_state$mate_genotype` values must match cube genotypes.", call. = FALSE)
  }
  initial_state
}

mosquito_state_from_table_one <- function(table, backend, parameters, sampler) {
  state <- numeric(backend$index$total_state_len)
  G <- backend$cube_info$G
  for (i in seq_len(NROW(table))) {
    row <- table[i, , drop = FALSE]
    node_i <- match(row$node[[1L]], names(parameters))
    genotype_i <- match(row$genotype[[1L]], backend$cube_info$genotypesID)
    count <- row$count[[1L]]
    if (sampler == "tau") {
      if (abs(count - round(count)) > 1e-8) {
        stop("Tau initial state counts must be integer-valued.", call. = FALSE)
      }
      count <- round(count)
    }
    if (row$stage[[1L]] == "egg") {
      substage <- if (is.na(row$substage[[1L]])) 1L else as.integer(row$substage[[1L]])
      if (is.na(substage) || substage < 1L || substage > backend$cfg$nE) {
        stop("Egg initial state rows require a valid `substage`.", call. = FALSE)
      }
      state[backend$index$egg_ix[substage, genotype_i, node_i]] <- state[backend$index$egg_ix[substage, genotype_i, node_i]] + count
    } else if (row$stage[[1L]] == "larva") {
      substage <- if (is.na(row$substage[[1L]])) 1L else as.integer(row$substage[[1L]])
      if (is.na(substage) || substage < 1L || substage > backend$cfg$nL) {
        stop("Larva initial state rows require a valid `substage`.", call. = FALSE)
      }
      state[backend$index$larv_ix[substage, genotype_i, node_i]] <- state[backend$index$larv_ix[substage, genotype_i, node_i]] + count
    } else if (row$stage[[1L]] == "pupa") {
      substage <- if (is.na(row$substage[[1L]])) 1L else as.integer(row$substage[[1L]])
      if (is.na(substage) || substage < 1L || substage > backend$cfg$nP) {
        stop("Pupa initial state rows require a valid `substage`.", call. = FALSE)
      }
      state[backend$index$pup_ix[substage, genotype_i, node_i]] <- state[backend$index$pup_ix[substage, genotype_i, node_i]] + count
    } else if (row$stage[[1L]] == "adult_male") {
      state[backend$index$male_ix[genotype_i, node_i]] <- state[backend$index$male_ix[genotype_i, node_i]] + count
    } else {
      status <- row$mating_status[[1L]]
      if (is.na(status) || status == "unmated") {
        state[backend$index$unm_ix[genotype_i, node_i]] <- state[backend$index$unm_ix[genotype_i, node_i]] + count
      } else if (status == "mated") {
        mate_i <- match(row$mate_genotype[[1L]], backend$cube_info$genotypesID)
        if (is.na(mate_i)) {
          stop("Mated adult-female initial state rows require `mate_genotype`.", call. = FALSE)
        }
        state[backend$index$fem_ix[genotype_i + (mate_i - 1L) * G, node_i]] <-
          state[backend$index$fem_ix[genotype_i + (mate_i - 1L) * G, node_i]] + count
      } else {
        stop("`initial_state$mating_status` must be `unmated`, `mated`, or NA.", call. = FALSE)
      }
    }
  }
  state
}

mosquito_restore_initial_state_table <- function(backends, initial_state, parameters, cube_info, sampler) {
  table <- mosquito_normalize_initial_state_table(initial_state, parameters, cube_info)
  for (species_i in seq_along(backends)) {
    species_table <- table[table$species == backends[[species_i]]$species_name, , drop = FALSE]
    state <- mosquito_state_from_table_one(species_table, backends[[species_i]], parameters, sampler)
    backends[[species_i]] <- mosquito_set_state_vector(backends[[species_i]], state, 0)
  }
  backends
}

#' Run a standalone native mosquito lifecycle simulation
#'
#' Simulates mosquito eggs, larvae, pupae, adult males, unmated adult females,
#' and mated adult females by genotype without creating human or mosquito
#' infection compartments. `sampler = "ode"` uses the native deterministic
#' entomology engine; `sampler = "tau"` uses a native entomology-only tau step.
#'
#' @param timesteps final timestep to simulate.
#' @param parameters a `mosquito_parameters` object for one node, or a named
#' list of such objects for a network.
#' @param sampler `"ode"` or `"tau"`.
#' @param releases optional data frame with `timestep`, `node`, `species`,
#' `stage`, `genotype`, `count`, and optional `substage`, `mate_genotype`.
#' @param initial_state `NULL`, or a checkpoint returned in `state`.
#' @param return_state whether to include a resumable state checkpoint.
#' @param seed optional R RNG seed used for `sampler = "tau"`.
#'
#' @return A list with `counts`, `state_counts`, `releases`, `state`, and
#' `metadata`.
#' @export
run_mosquito_simulation <- function(
  timesteps,
  parameters = get_mosquito_parameters(),
  sampler = "ode",
  releases = NULL,
  initial_state = NULL,
  return_state = TRUE,
  seed = NULL
) {
  timesteps <- as.integer(timesteps)
  if (length(timesteps) != 1L || is.na(timesteps) || timesteps < 0L) {
    stop("`timesteps` must be a non-negative integer scalar.", call. = FALSE)
  }
  sampler <- match.arg(sampler, c("ode", "tau"))
  if (!is.null(seed)) {
    if (sampler != "tau") {
      warning("`seed` is ignored for deterministic standalone mosquito ODE runs.", call. = FALSE)
    } else {
      set.seed(seed)
    }
  }
  parameters <- mosquito_normalize_network(parameters)
  cube_info <- cube_genotype_info(parameters[[1L]]$cube)
  releases <- mosquito_normalize_releases(releases, parameters, cube_info)

  backends <- lapply(seq_along(parameters[[1L]]$species), function(species_i) {
    mosquito_build_species_backend(parameters, species_i, sampler, timesteps)
  })

  start_t <- 0L
  restored <- mosquito_restore_checkpoint(backends, initial_state, sampler)
  if (!is.null(restored)) {
    backends <- restored$backends
    start_t <- restored$timestep
    if (timesteps < start_t) {
      stop("`timesteps` must be at least the checkpoint timestep.", call. = FALSE)
    }
  } else if (!is.null(initial_state)) {
    backends <- mosquito_restore_initial_state_table(backends, initial_state, parameters, cube_info, sampler)
  }

  state_counts <- list()
  output_i <- 1L
  record <- function(timestep) {
    dfs <- lapply(backends, function(b) {
      mosquito_state_counts_one(mosquito_state_vector(b), b, names(parameters), timestep)
    })
    state_counts[[output_i]] <<- do.call(rbind, dfs)
    output_i <<- output_i + 1L
  }

  backends <- mosquito_apply_releases_at(backends, releases, parameters, start_t, sampler)
  record(start_t)
  if (timesteps > start_t) {
    for (t in seq.int(start_t + 1L, timesteps)) {
      for (species_i in seq_along(backends)) {
        backends[[species_i]] <- mosquito_step_backend(backends[[species_i]], t, sampler, parameters)
      }
      backends <- mosquito_apply_releases_at(backends, releases, parameters, t, sampler)
      record(t)
    }
  }

  state_counts <- do.call(rbind, state_counts)
  counts <- mosquito_aggregate_counts(state_counts)
  list(
    counts = counts,
    state_counts = state_counts,
    releases = releases,
    state = if (return_state) mosquito_checkpoint(backends, sampler, timesteps) else NULL,
    metadata = list(
      sampler = sampler,
      nodes = names(parameters),
      species = parameters[[1L]]$species,
      genotypes = cube_info$genotypesID,
      stage_units = "mosquito count",
      timing = "releases are applied at timestep boundaries before output"
    )
  )
}
