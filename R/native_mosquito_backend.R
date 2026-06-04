native_mosquito_backend_enabled <- function(parameters) {
  if (is.null(parameters$native_mosquito_backend)) {
    return(FALSE)
  }
  isTRUE(parameters$native_mosquito_backend)
}

native_backend_warning_state <- new.env(parent = emptyenv())

legacy_individual_mosquito_backend_enabled <- function(parameters) {
  isTRUE(parameters$individual_mosquitoes) && !native_mosquito_backend_enabled(parameters)
}

legacy_individual_mosquito_events_enabled <- function(parameters, events = NULL) {
  enabled <- legacy_individual_mosquito_backend_enabled(parameters)
  if (is.null(events)) {
    return(enabled)
  }

  has_infection <- !is.null(events$mosquito_infection)
  has_death <- !is.null(events$mosquito_death)
  if (xor(has_infection, has_death)) {
    stop(
      paste(
        "Mosquito event wiring is inconsistent:",
        "`mosquito_infection` and `mosquito_death` must be created together."
      )
    )
  }
  if ((has_infection && has_death) != enabled) {
    stop("Mosquito event wiring is inconsistent with the selected mosquito backend.")
  }
  enabled
}

species_beta_value <- function(parameters, species_i) {
  beta <- parameters$beta
  if (!is.atomic(beta) || !is.numeric(beta) || length(beta) == 0L) {
    stop("parameters$beta must be a numeric scalar or one value per mosquito species.")
  }

  if (length(beta) == 1L) {
    return(as.numeric(beta[[1L]]))
  }

  species_names <- as.character(parameters$species)
  if (!is.null(names(beta))) {
    missing_species <- setdiff(species_names, names(beta))
    if (length(missing_species) == 0L) {
      return(as.numeric(beta[[species_names[[species_i]]]]))
    }
  }

  if (length(beta) != length(species_names)) {
    stop("parameters$beta must be a numeric scalar or have one value per mosquito species.")
  }
  as.numeric(beta[[species_i]])
}

native_warn_tau_leap_backend <- function(parameters) {
  if (!native_mosquito_backend_enabled(parameters) || !isTRUE(parameters$individual_mosquitoes)) {
    return(invisible(NULL))
  }
  if (isTRUE(native_backend_warning_state$tau_leap_backend)) {
    return(invisible(NULL))
  }
  warning(
    paste(
      "`native_mosquito_backend = TRUE` with `individual_mosquitoes = TRUE`",
      "uses the native count-based tau-leap mosquito engine,",
      "not the legacy individual mosquito event backend."
    ),
    call. = FALSE
  )
  native_backend_warning_state$tau_leap_backend <- TRUE
  invisible(NULL)
}

native_metapop_backend_enabled <- function(parameters) {
  enabled <- vapply(parameters, native_mosquito_backend_enabled, logical(1))
  if (all(enabled)) {
    return(TRUE)
  }
  if (any(enabled)) {
    stop("All metapopulation parameter sets must agree on native_mosquito_backend.")
  }
  FALSE
}

native_mosquito_stage_config <- function(parameters) {
  list(
    nE = as.integer(if (is.null(parameters$native_mosquito_nE)) 1L else parameters$native_mosquito_nE),
    nL = as.integer(if (is.null(parameters$native_mosquito_nL)) 1L else parameters$native_mosquito_nL),
    nP = as.integer(if (is.null(parameters$native_mosquito_nP)) 1L else parameters$native_mosquito_nP),
    nEIP = as.integer(if (is.null(parameters$native_mosquito_nEIP)) 1L else parameters$native_mosquito_nEIP),
    nu = as.numeric(if (is.null(parameters$native_mosquito_nu)) 1 else parameters$native_mosquito_nu),
    dt_stoch = as.numeric(if (is.null(parameters$mosquito_tau_step)) 0.1 else parameters$mosquito_tau_step)
  )
}

native_assert_identical_scalar <- function(values, label) {
  base <- values[[1]]
  if (length(values) > 1L && any(!vapply(values[-1L], function(x) isTRUE(all.equal(x, base)), logical(1)))) {
    stop(sprintf("Native mosquito backend requires identical `%s` across shared mosquito nodes.", label))
  }
  base
}

native_assert_identical_configs <- function(parameters, cfg) {
  for (i in seq_along(parameters)) {
    cfg_i <- native_mosquito_stage_config(parameters[[i]])
    if (!isTRUE(all.equal(cfg_i, cfg))) {
      stop("All populations must use the same native mosquito stage configuration.")
    }
  }
}

native_extract_move_value <- function(parameters, primary, fallback) {
  value <- parameters[[primary]]
  if (is.null(value)) {
    value <- parameters[[fallback]]
  }
  value
}

native_collect_shared_value <- function(parameters, primary, fallback = NULL) {
  out <- NULL
  for (i in seq_along(parameters)) {
    value <- parameters[[i]][[primary]]
    if (is.null(value) && !is.null(fallback)) {
      value <- parameters[[i]][[fallback]]
    }
    if (is.null(value)) {
      next
    }
    if (is.null(out)) {
      out <- value
      next
    }
    if (!isTRUE(all.equal(value, out, check.attributes = FALSE))) {
      stop(sprintf("All populations must share the same `%s` when native mosquito movement is enabled.", primary))
    }
  }
  out
}

human_mobility_enabled_any <- function(parameters) {
  any(vapply(parameters, function(p) isTRUE(p$human_mobility_enabled), logical(1)))
}

human_mobility_value <- function(parameters, name, default) {
  value <- parameters[[name]]
  if (is.null(value)) {
    return(default)
  }
  value
}

human_mobility_validate_scalar_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name))
  }
}

human_mobility_validate_trip_duration <- function(parameters) {
  duration_type <- human_mobility_value(parameters, "human_trip_duration_type", "fixed")
  if (!is.character(duration_type) ||
      length(duration_type) != 1L ||
      !(duration_type %in% c("fixed", "geometric"))) {
    stop('`human_trip_duration_type` must be "fixed" or "geometric".')
  }

  duration_mean <- human_mobility_value(parameters, "human_trip_duration_mean", 1)
  if (!is.numeric(duration_mean) ||
      length(duration_mean) != 1L ||
      is.na(duration_mean) ||
      !is.finite(duration_mean)) {
    stop("`human_trip_duration_mean` must be a finite numeric scalar.")
  }

  if (identical(duration_type, "fixed") &&
      (duration_mean < 1 || !isTRUE(all.equal(duration_mean, as.integer(duration_mean))))) {
    stop('`human_trip_duration_mean` must be a positive integer >= 1 when `human_trip_duration_type = "fixed"`.')
  }

  if (identical(duration_type, "geometric") && duration_mean < 1) {
    stop('`human_trip_duration_mean` must be numeric and >= 1 when `human_trip_duration_type = "geometric"`.')
  }
}

human_mobility_validate_matrix <- function(move_probs, n_nodes) {
  if (!is.matrix(move_probs) || !is.numeric(move_probs)) {
    stop("`human_move_probs` must be a numeric matrix.")
  }
  if (nrow(move_probs) != ncol(move_probs)) {
    stop("`human_move_probs` must be square.")
  }
  if (!all(dim(move_probs) == c(n_nodes, n_nodes))) {
    stop("`human_move_probs` must have dimension n_nodes x n_nodes.")
  }
  if (anyNA(move_probs)) {
    stop("`human_move_probs` must not contain missing values.")
  }
  if (!all(is.finite(move_probs))) {
    stop("`human_move_probs` must contain only finite values.")
  }
  if (any(move_probs < 0)) {
    stop("`human_move_probs` must not contain negative entries.")
  }
  if (!all(vlapply(seq_len(n_nodes), function(i) approx_sum(move_probs[i, ], 1)))) {
    stop("Each row of `human_move_probs` must sum to 1.")
  }
}

human_mobility_collect_move_probs <- function(parameters, n_nodes) {
  move_probs <- Filter(
    Negate(is.null),
    lapply(parameters, function(p) p$human_move_probs)
  )
  if (length(move_probs) == 0L) {
    stop("`human_move_probs` must be provided when `human_mobility_enabled = TRUE`.")
  }

  for (matrix_i in seq_along(move_probs)) {
    human_mobility_validate_matrix(move_probs[[matrix_i]], n_nodes)
  }

  shared <- move_probs[[1L]]
  if (length(move_probs) > 1L) {
    for (matrix_i in seq_along(move_probs)[-1L]) {
      if (!isTRUE(all.equal(
        move_probs[[matrix_i]],
        shared,
        check.attributes = FALSE,
        tolerance = sqrt(.Machine$double.eps)
      ))) {
        stop("All non-NULL `human_move_probs` matrices must be identical.")
      }
    }
  }

  shared
}

human_mobility_resolve_move_probs <- function(parameters, n_nodes = length(parameters)) {
  if (!human_mobility_enabled_any(parameters)) {
    return(NULL)
  }
  human_mobility_collect_move_probs(parameters, n_nodes)
}

human_mobility_validate_parameters <- function(parameters, n_nodes = length(parameters)) {
  for (i in seq_along(parameters)) {
    human_mobility_validate_scalar_logical(
      human_mobility_value(parameters[[i]], "human_mobility_enabled", FALSE),
      "human_mobility_enabled"
    )
    if (!identical(human_mobility_value(parameters[[i]], "human_mobility_mode", "explicit"), "explicit")) {
      stop('`human_mobility_mode` must be "explicit".')
    }
    human_mobility_validate_trip_duration(parameters[[i]])
    if (isTRUE(human_mobility_value(parameters[[i]], "human_mobility_enabled", FALSE))) {
      de <- human_mobility_value(parameters[[i]], "de", NA_real_)
      if (!is.numeric(de) ||
          length(de) != 1L ||
          is.na(de) ||
          !is.finite(de) ||
          de <= 0) {
        stop("Explicit human mobility requires `de` to be a positive numeric value greater than 0.")
      }
    }
    if (!is.null(parameters[[i]]$human_move_rates)) {
      stop("`human_move_rates` is not supported for Stage 1 explicit human mobility.")
    }
  }

  if (human_mobility_enabled_any(parameters)) {
    human_mobility_collect_move_probs(parameters, n_nodes)
  }
}

human_mobility_identity_matrix <- function(matrix, n_nodes) {
  matrix <- as.matrix(matrix)
  all(dim(matrix) == c(n_nodes, n_nodes)) &&
    all(abs(matrix - diag(n_nodes)) < sqrt(.Machine$double.eps))
}

human_mobility_validate_identity_mixing <- function(mixing, n_nodes, label) {
  for (i in seq_along(mixing)) {
    if (!human_mobility_identity_matrix(mixing[[i]], n_nodes)) {
      stop(sprintf("Explicit human mobility requires identity `%s` matrices.", label))
    }
  }
}

human_mobility_validate_zero_capture <- function(p_captured, n_nodes, p_success) {
  for (i in seq_along(p_captured)) {
    p_captured_i <- as.matrix(p_captured[[i]])
    if (!all(dim(p_captured_i) == c(n_nodes, n_nodes))) {
      stop("Explicit human mobility requires `p_captured` matrices to be n_nodes x n_nodes.")
    }
    if (any(abs(p_captured_i) >= sqrt(.Machine$double.eps))) {
      if (identical(p_success, 0)) {
        stop("Explicit human mobility does not support nonzero `p_captured`.")
      }
      stop("Explicit human mobility does not support active border test-and-treat capture.")
    }
  }
}

validate_explicit_human_mobility_metapop <- function(
  parameters,
  export_mixing,
  import_mixing,
  p_captured,
  p_success
) {
  n_nodes <- length(parameters)
  human_mobility_validate_parameters(parameters, n_nodes)
  if (!human_mobility_enabled_any(parameters)) {
    return(invisible(NULL))
  }

  if (n_nodes < 2L || !all(vapply(parameters, native_mosquito_backend_enabled, logical(1)))) {
    stop("Explicit human mobility requires the native metapop mosquito backend.")
  }
  human_mobility_validate_identity_mixing(export_mixing, n_nodes, "export_mixing")
  human_mobility_validate_identity_mixing(import_mixing, n_nodes, "import_mixing")
  human_mobility_validate_zero_capture(p_captured, n_nodes, p_success)

  invisible(NULL)
}

native_validate_human_movement <- function(parameters) {
  human_mobility_validate_parameters(parameters)
  if (human_mobility_enabled_any(parameters) && length(parameters) < 2L) {
    stop("Explicit human mobility requires the native metapop mosquito backend.")
  }
}

native_collect_mosquito_movement <- function(parameters) {
  n_nodes <- length(parameters)
  move_probs <- native_collect_shared_value(parameters, "move_probs", "mosquito_move_probs")
  move_rates <- native_collect_shared_value(parameters, "move_rates", "mosquito_move_rates")

  if (n_nodes == 1L || is.null(move_probs) || is.null(move_rates)) {
    return(list(
      has_move = FALSE,
      move_probs = matrix(0, nrow = n_nodes, ncol = n_nodes),
      move_rates = rep(0, n_nodes)
    ))
  }

  move_probs <- as.matrix(move_probs)
  move_rates <- as.numeric(move_rates)
  if (length(move_rates) == 1L) {
    move_rates <- rep(move_rates, n_nodes)
  }
  if (!all(dim(move_probs) == c(n_nodes, n_nodes))) {
    stop("Mosquito movement matrix must be n_pop x n_pop for the shared native backend.")
  }
  if (length(move_rates) != n_nodes) {
    stop("Mosquito movement rates must have length n_pop for the shared native backend.")
  }

  list(
    has_move = any(move_rates > 0) && any(move_probs != 0),
    move_probs = move_probs,
    move_rates = move_rates
  )
}

build_native_epi_indices <- function(G, nNodes, nE, nL, nP, nEIP) {
  idx <- 1L
  take_block <- function(len, dims) {
    block <- array(seq.int(idx, length.out = len), dim = dims)
    idx <<- idx + len
    block
  }

  nStages <- nEIP + 2L
  egg_ix <- take_block(nE * G * nNodes, c(nE, G, nNodes))
  larv_ix <- take_block(nL * G * nNodes, c(nL, G, nNodes))
  pup_ix <- take_block(nP * G * nNodes, c(nP, G, nNodes))
  male_ix <- take_block(G * nNodes, c(G, nNodes))
  unm_ix <- take_block(G * nNodes, c(G, nNodes))
  fem_ix <- take_block(G * G * nStages * nNodes, c(G, G, nStages, nNodes))

  hS_ix <- integer(nNodes)
  hI_ix <- integer(nNodes)
  for (node in seq_len(nNodes)) {
    hS_ix[[node]] <- idx
    idx <- idx + 1L
    hI_ix[[node]] <- idx
    idx <- idx + 1L
  }

  list(
    egg_ix = egg_ix,
    larv_ix = larv_ix,
    pup_ix = pup_ix,
    male_ix = male_ix,
    unm_ix = unm_ix,
    fem_ix = fem_ix,
    hS_ix = hS_ix,
    hI_ix = hI_ix,
    nStages = nStages,
    total_state_len = idx - 1L
  )
}

native_align_cube_vec <- function(x, key, default) {
  if (is.null(x)) {
    return(rep(default, length(key)))
  }
  if (!is.null(names(x))) {
    return(as.numeric(x[key]))
  }
  x <- as.numeric(x)
  if (length(x) == 1L) {
    return(rep(x, length(key)))
  }
  x
}

native_align_cube_mat <- function(x, key, default = 1) {
  if (is.null(x)) {
    return(matrix(default, nrow = length(key), ncol = length(key)))
  }
  if (!is.null(rownames(x)) && !is.null(colnames(x))) {
    return(as.matrix(x[key, key, drop = FALSE]))
  }
  as.matrix(x)
}

native_align_cube_arr3 <- function(x, key) {
  if (is.null(x)) {
    return(NULL)
  }
  dn <- dimnames(x)
  if (!is.null(dn) && length(dn) >= 3) {
    return(x[key, key, key, drop = FALSE])
  }
  x
}

native_build_birth_matrix <- function(cube, cube_info) {
  G <- cube_info$G
  g <- cube_info$genotypesID
  if (is.null(cube)) {
    B_mat <- matrix(0, nrow = G * G, ncol = G)
    B_mat[seq_len(G), cube_info$wild_type_index] <- 1
    return(B_mat)
  }

  s_fec <- native_align_cube_vec(cube$s, g, 1)
  ih <- native_align_cube_arr3(cube$ih, g)
  tau <- native_align_cube_arr3(cube$tau, g)
  if (is.null(ih)) {
    B_mat <- matrix(1 / G, nrow = G * G, ncol = G)
    for (f in seq_len(G)) {
      B_mat[((f - 1L) * G + 1L):(f * G), f] <- s_fec[[f]] / G
    }
    return(B_mat)
  }
  if (is.null(tau)) {
    tau <- array(1, dim = dim(ih))
  }
  o_prob <- ih * tau
  B_base <- array(0, dim = c(G, G, G))
  for (f in seq_len(G)) {
    B_base[f, , ] <- s_fec[[f]] * o_prob[f, , ]
  }
  matrix(B_base, nrow = G * G, ncol = G)
}

native_warn_equilibrium_fallback <- function(message) {
  if (isTRUE(native_backend_warning_state$equilibrium_fallback)) {
    return(invisible(NULL))
  }
  warning(message, call. = FALSE)
  native_backend_warning_state$equilibrium_fallback <- TRUE
  invisible(NULL)
}

native_carrying_capacity_scale <- function(parameters, species_i, timestep = 0L) {
  scale <- 1
  if (isTRUE(parameters$carrying_capacity) &&
      !is.null(parameters$carrying_capacity_timesteps) &&
      length(parameters$carrying_capacity_timesteps) > 0L) {
    changes <- which(parameters$carrying_capacity_timesteps <= timestep)
    if (length(changes) > 0L) {
      scale <- scale * parameters$carrying_capacity_scalers[max(changes), species_i]
    }
  }

  scale * carrying_capacity(
    timestep,
    parameters$model_seasonality,
    parameters$g0,
    parameters$g,
    parameters$h,
    1,
    calculate_R_bar(parameters),
    parameters$rainfall_floor
  )
}

native_exact_equilibrium_node <- function(parameters, species_i, cube_info, timestep = 0L) {
  cfg <- native_mosquito_stage_config(parameters)
  G <- cube_info$G
  wt <- cube_info$wild_type_index
  gids <- cube_info$genotypesID
  wt_pair <- wt + (wt - 1L) * G

  m <- parameters$total_M * parameters$species_proportions[[species_i]]
  if (!is.finite(m) || m <= 0) {
    return(list(
      valid = TRUE,
      k0 = 0,
      egg = rep(0, cfg$nE),
      larv = rep(0, cfg$nL),
      pup = rep(0, cfg$nP),
      male = 0,
      female_S = 0,
      female_E = rep(0, cfg$nEIP),
      female_I = 0
    ))
  }

  B_mat <- native_build_birth_matrix(parameters$cube, cube_info)
  birth_weights <- as.numeric(B_mat[wt_pair, ])
  expected_birth <- rep(0, G)
  expected_birth[[wt]] <- 1
  if (!isTRUE(all.equal(birth_weights, expected_birth, tolerance = 1e-10))) {
    native_warn_equilibrium_fallback(
      paste(
        "Native mosquito equilibrium fallback:",
        "wild-type self-cross does not produce only wild-type offspring,",
        "so the backend is using the legacy approximate initialization."
      )
    )
    return(list(valid = FALSE))
  }

  omega <- cube_omega_vector(parameters$cube, G, gids)
  omega_inv <- ifelse(omega == 0, 1e3, 1 / omega)
  phi <- native_align_cube_vec(if (is.null(parameters$cube)) NULL else parameters$cube$phi, gids, 0.5)
  xiF <- native_align_cube_vec(if (is.null(parameters$cube)) NULL else parameters$cube$xiF, gids, 1)
  xiM <- native_align_cube_vec(if (is.null(parameters$cube)) NULL else parameters$cube$xiM, gids, 1)
  c_vec <- native_align_cube_vec(if (is.null(parameters$cube)) NULL else parameters$cube$c, gids, 1)

  eq_traits <- equilibrium_species_traits(parameters, species_i)
  muF_eff <- eq_traits$muF * omega_inv[[wt]]
  muM_eff <- parameters$mum[[species_i]] * omega_inv[[wt]]
  if (!is.finite(muF_eff) || muF_eff <= 0 ||
      !is.finite(muM_eff) || muM_eff <= 0 ||
      phi[[wt]] <= 0 || xiF[[wt]] <= 0 || xiM[[wt]] <= 0) {
    native_warn_equilibrium_fallback(
      "Native mosquito equilibrium fallback: invalid wild-type adult traits for exact initialization."
    )
    return(list(valid = FALSE))
  }

  rE <- cfg$nE / parameters$del
  rL <- cfg$nL / parameters$dl
  rP <- cfg$nP / parameters$dpl
  rEIP <- cfg$nEIP / parameters$dem
  if (!all(is.finite(c(rE, rL, rP))) || any(c(rE, rL, rP) <= 0)) {
    native_warn_equilibrium_fallback(
      "Native mosquito equilibrium fallback: invalid native stage progression rates."
    )
    return(list(valid = FALSE))
  }

  beta_eff <- eggs_laid(
    species_beta_value(parameters, species_i),
    eq_traits$muF,
    eq_traits$f
  )

  aE <- rE / (parameters$me + rE)
  egg_1 <- beta_eff * m / (parameters$me + rE)
  egg <- egg_1 * aE^(0:(cfg$nE - 1L))
  egg_last <- egg[[cfg$nE]]

  aP <- rP / (parameters$mup + rP)
  pup_last <- muF_eff * m / (phi[[wt]] * xiF[[wt]] * rP)
  pup_1 <- pup_last / aP^(cfg$nP - 1L)
  pup <- pup_1 * aP^(0:(cfg$nP - 1L))

  larv_last <- (parameters$mup + rP) * pup_1 / rL
  dd_rate <- (rE * egg_last * rL^(cfg$nL - 1L) / larv_last)^(1 / cfg$nL) - rL
  if (!is.finite(dd_rate) || dd_rate <= parameters$ml) {
    native_warn_equilibrium_fallback(
      "Native mosquito equilibrium fallback: failed to derive a positive density-dependent equilibrium."
    )
    return(list(valid = FALSE))
  }
  aL <- rL / (dd_rate + rL)
  larv_1 <- rE * egg_last / (dd_rate + rL)
  larv <- larv_1 * aL^(0:(cfg$nL - 1L))
  K_eff <- sum(larv) / (dd_rate / parameters$ml - 1)
  scale_t0 <- native_carrying_capacity_scale(parameters, species_i, timestep)
  if (!is.finite(K_eff) || K_eff <= 0 || !is.finite(scale_t0) || scale_t0 <= 0) {
    native_warn_equilibrium_fallback(
      "Native mosquito equilibrium fallback: failed to derive a valid carrying capacity."
    )
    return(list(valid = FALSE))
  }

  lambda <- c_vec[[wt]] * parameters$init_foim
  female_emerge <- muF_eff * m
  female_S <- female_emerge / (lambda + muF_eff)
  if (cfg$nEIP > 0L) {
    aI <- rEIP / (muF_eff + rEIP)
    female_E1 <- lambda * female_S / (muF_eff + rEIP)
    female_E <- female_E1 * aI^(0:(cfg$nEIP - 1L))
    female_I <- rEIP * female_E[[cfg$nEIP]] / muF_eff
  } else {
    female_E <- numeric(0)
    female_I <- lambda * female_S / muF_eff
  }

  list(
    valid = TRUE,
    k0 = K_eff / scale_t0,
    egg = egg,
    larv = larv,
    pup = pup,
    male = ((1 - phi[[wt]]) * xiM[[wt]] * rP * pup_last) / muM_eff,
    female_S = female_S,
    female_E = female_E,
    female_I = female_I
  )
}

native_effective_carrying_capacity <- function(parameters, species_i, timestep) {
  cube_info <- cube_genotype_info(parameters$cube)
  eq <- native_exact_equilibrium_node(parameters, species_i, cube_info, timestep = 0L)
  if (isTRUE(eq$valid)) {
    k0 <- eq$k0
  } else {
    p <- parameters$species_proportions[[species_i]]
    m <- p * parameters$total_M
    k0 <- calculate_carrying_capacity(parameters, m, species_i)
  }

  if (isTRUE(parameters$carrying_capacity) &&
      !is.null(parameters$carrying_capacity_timesteps) &&
      length(parameters$carrying_capacity_timesteps) > 0L) {
    scaler <- 1
    changes <- which(parameters$carrying_capacity_timesteps <= timestep)
    if (length(changes) > 0L) {
      scaler <- parameters$carrying_capacity_scalers[max(changes), species_i]
    }
    k0 <- k0 * scaler
  }

  carrying_capacity(
    timestep,
    parameters$model_seasonality,
    parameters$g0,
    parameters$g,
    parameters$h,
    k0,
    calculate_R_bar(parameters),
    parameters$rainfall_floor
  )
}

native_prepare_carrying_capacity_node <- function(parameters, species_i, cube_info) {
  eq <- native_exact_equilibrium_node(parameters, species_i, cube_info, timestep = 0L)
  if (isTRUE(eq$valid)) {
    k0 <- eq$k0
  } else {
    p <- parameters$species_proportions[[species_i]]
    m <- p * parameters$total_M
    k0 <- calculate_carrying_capacity(parameters, m, species_i)
  }

  list(
    k0 = k0,
    model_seasonality = parameters$model_seasonality,
    g0 = parameters$g0,
    g = parameters$g,
    h = parameters$h,
    R_bar = calculate_R_bar(parameters),
    rainfall_floor = parameters$rainfall_floor,
    carrying_capacity_timesteps = if (isTRUE(parameters$carrying_capacity) &&
      !is.null(parameters$carrying_capacity_timesteps) &&
      length(parameters$carrying_capacity_timesteps) > 0L) {
      as.integer(parameters$carrying_capacity_timesteps)
    } else {
      integer(0)
    },
    carrying_capacity_scalers = if (isTRUE(parameters$carrying_capacity) &&
      !is.null(parameters$carrying_capacity_scalers) &&
      NROW(parameters$carrying_capacity_scalers) > 0L) {
      as.numeric(parameters$carrying_capacity_scalers[, species_i])
    } else {
      numeric(0)
    }
  )
}

native_build_carrying_capacity_lookup <- function(parameters, species_i, timesteps) {
  cube_info <- cube_genotype_info(parameters[[1L]]$cube)
  list(
    timesteps = as.integer(timesteps),
    nodes = lapply(parameters, native_prepare_carrying_capacity_node, species_i = species_i, cube_info = cube_info)
  )
}

native_carrying_capacity_at <- function(lookup, timestep) {
  timestep <- max(0L, as.integer(timestep))
  vapply(
    lookup$nodes,
    function(node) {
      k0 <- node$k0
      if (length(node$carrying_capacity_timesteps) > 0L) {
        changes <- which(node$carrying_capacity_timesteps <= timestep)
        if (length(changes) > 0L) {
          k0 <- k0 * node$carrying_capacity_scalers[[max(changes)]]
        }
      }

      carrying_capacity(
        timestep,
        node$model_seasonality,
        node$g0,
        node$g,
        node$h,
        k0,
        node$R_bar,
        node$rainfall_floor
      )
    },
    numeric(1)
  )
}

native_initial_state <- function(parameters, species_i, cube_info, index) {
  state <- numeric(index$total_state_len)
  wt <- cube_info$wild_type_index
  s_stage <- 1L
  e_stage <- if (index$nStages > 2L) 2L else NA_integer_
  i_stage <- index$nStages

  for (node in seq_along(parameters)) {
    eq <- native_exact_equilibrium_node(parameters[[node]], species_i, cube_info, timestep = 0L)
    if (isTRUE(eq$valid)) {
      state[index$egg_ix[, wt, node]] <- eq$egg
      state[index$larv_ix[, wt, node]] <- eq$larv
      state[index$pup_ix[, wt, node]] <- eq$pup
      state[index$male_ix[wt, node]] <- eq$male
      state[index$fem_ix[wt, wt, s_stage, node]] <- eq$female_S
      if (!is.na(e_stage) && length(eq$female_E) > 0L) {
        for (j in seq_along(eq$female_E)) {
          state[index$fem_ix[wt, wt, s_stage + j, node]] <- eq$female_E[[j]]
        }
      }
      state[index$fem_ix[wt, wt, i_stage, node]] <- eq$female_I
    } else {
      init_counts <- initial_mosquito_counts(
        parameters[[node]],
        species_i,
        parameters[[node]]$init_foim,
        parameters[[node]]$total_M * parameters[[node]]$species_proportions[[species_i]]
      )
      adult_total <- sum(init_counts[ADULT_ODE_INDICES])

      state[index$egg_ix[1L, wt, node]] <- init_counts[[ODE_INDICES[["E"]]]]
      state[index$larv_ix[1L, wt, node]] <- init_counts[[ODE_INDICES[["L"]]]]
      state[index$pup_ix[1L, wt, node]] <- init_counts[[ODE_INDICES[["P"]]]]
      state[index$male_ix[wt, node]] <- adult_total
      state[index$fem_ix[wt, wt, s_stage, node]] <- init_counts[[ADULT_ODE_INDICES[["Sm"]]]]
      if (!is.na(e_stage)) {
        state[index$fem_ix[wt, wt, e_stage, node]] <- init_counts[[ADULT_ODE_INDICES[["Pm"]]]]
      }
      state[index$fem_ix[wt, wt, i_stage, node]] <- init_counts[[ADULT_ODE_INDICES[["Im"]]]]
    }
    state[index$hS_ix[[node]]] <- 0
    state[index$hI_ix[[node]]] <- 1
  }

  state
}

native_extract_summary <- function(state, model, node = model$node) {
  index <- model$shared$index
  G <- model$cube_info$G
  nStages <- index$nStages
  s_stage <- 1L
  i_stage <- nStages
  eip_stages <- if (nStages > 2L) 2L:(nStages - 1L) else integer(0)

  E_g <- vapply(seq_len(G), function(g) sum(state[index$egg_ix[, g, node]]), numeric(1))
  L_g <- vapply(seq_len(G), function(g) sum(state[index$larv_ix[, g, node]]), numeric(1))
  P_g <- vapply(seq_len(G), function(g) sum(state[index$pup_ix[, g, node]]), numeric(1))
  M_g <- state[index$male_ix[, node]]
  U_g <- state[index$unm_ix[, node]]

  female_g <- numeric(G)
  infectious_g <- numeric(G)
  Sm_total <- sum(U_g)
  Pm_total <- 0
  Im_total <- 0

  for (gf in seq_len(G)) {
    female_g[[gf]] <- female_g[[gf]] + U_g[[gf]]
    for (gm in seq_len(G)) {
      Sm_ij <- state[index$fem_ix[gf, gm, s_stage, node]]
      female_g[[gf]] <- female_g[[gf]] + Sm_ij
      Sm_total <- Sm_total + Sm_ij
      if (length(eip_stages) > 0L) {
        Pm_ij <- sum(state[index$fem_ix[gf, gm, eip_stages, node]])
        female_g[[gf]] <- female_g[[gf]] + Pm_ij
        Pm_total <- Pm_total + Pm_ij
      }
      Im_ij <- state[index$fem_ix[gf, gm, i_stage, node]]
      female_g[[gf]] <- female_g[[gf]] + Im_ij
      infectious_g[[gf]] <- infectious_g[[gf]] + Im_ij
      Im_total <- Im_total + Im_ij
    }
  }

  names(E_g) <- model$cube_info$genotypesID
  names(L_g) <- model$cube_info$genotypesID
  names(P_g) <- model$cube_info$genotypesID
  names(M_g) <- model$cube_info$genotypesID
  names(U_g) <- model$cube_info$genotypesID
  names(female_g) <- model$cube_info$genotypesID
  names(infectious_g) <- model$cube_info$genotypesID

  list(
    totals = c(E = sum(E_g), L = sum(L_g), P = sum(P_g), Sm = Sm_total, Pm = Pm_total, Im = Im_total),
    aquatic = list(E = E_g, L = L_g, P = P_g),
    male = M_g,
    unmated = U_g,
    female = female_g,
    infectious = infectious_g
  )
}

native_solver_compat_states <- function(state, model, stochastic) {
  summary <- native_extract_summary(state, model)
  if (!stochastic) {
    return(unname(summary$totals[c("E", "L", "P", "Sm", "Pm", "Im")]))
  }
  if (model$cube_info$G == 1L) {
    return(unname(summary$totals[c("E", "L", "P")]))
  }
  as.vector(rbind(summary$aquatic$E, summary$aquatic$L, summary$aquatic$P))
}

native_release_mate_allocation <- function(count, probs, stochastic) {
  if (stochastic) {
    return(sample_genotype_counts(count, probs))
  }
  as.numeric(count) * probs
}

native_apply_release_to_state <- function(state, model, sex, genotype_idx, count, stochastic, node = model$node) {
  if (count <= 0L) {
    return(state)
  }
  if (!(sex %in% c("M", "F"))) {
    stop("Mosquito releases must use sex='M' or sex='F'.")
  }

  if (sex == "M") {
    state[model$shared$index$male_ix[genotype_idx, node]] <-
      state[model$shared$index$male_ix[genotype_idx, node]] + count
    return(state)
  }

  male_counts <- native_extract_summary(state, model, node = node)$male
  mate_weights <- as.numeric(model$shared$eta[genotype_idx, ]) * as.numeric(male_counts)
  if (sum(mate_weights) > 0) {
    probs <- mate_weights / sum(mate_weights)
    alloc <- native_release_mate_allocation(count, probs, stochastic)
    for (gm in seq_len(model$cube_info$G)) {
      state[model$shared$index$fem_ix[genotype_idx, gm, 1L, node]] <-
        state[model$shared$index$fem_ix[genotype_idx, gm, 1L, node]] + alloc[[gm]]
    }
  } else {
    state[model$shared$index$unm_ix[genotype_idx, node]] <-
      state[model$shared$index$unm_ix[genotype_idx, node]] + count
  }

  state
}

native_empty_pending_inputs <- function(n_nodes) {
  list(
    timestep = NULL,
    beta = rep(NA_real_, n_nodes),
    mu = rep(NA_real_, n_nodes),
    foim = rep(NA_real_, n_nodes),
    ready = rep(FALSE, n_nodes),
    stepped = rep(FALSE, n_nodes)
  )
}

native_shared_state <- function(shared) {
  if (shared$stochastic) {
    return(shared$state)
  }
  solver_get_states(shared$solver_ptr)
}

native_validate_shared_state <- function(shared, state) {
  expected_len <- shared$index$total_state_len
  actual_len <- length(state)
  if (actual_len != expected_len) {
    stop(
      sprintf(
        paste(
          "Native mosquito state length mismatch during restore/set:",
          "expected %d entries but received %d.",
          "This usually means the saved state was created with an incompatible",
          "stage/genotype/node configuration."
        ),
        expected_len,
        actual_len
      )
    )
  }
  invisible(NULL)
}

native_set_shared_state <- function(shared, state, t = shared$t) {
  native_validate_shared_state(shared, state)
  if (shared$stochastic) {
    shared$state <- state
  } else {
    solver_set_states(shared$solver_ptr, t, state)
  }
  shared$t <- t
  invisible(NULL)
}

native_execute_shared_step <- function(shared) {
  pending <- shared$pending_inputs
  if (is.null(pending$timestep) || !all(pending$ready)) {
    stop("Native mosquito inputs were not provided for every node before stepping.")
  }

  K_t <- native_carrying_capacity_at(shared$carrying_capacity_lookup, pending$timestep)
  gamma_t <- as.numeric(shared$gamma_dd)

  if (shared$stochastic) {
    stoch_epi_engine_set_runtime(
      shared$engine_ptr,
      pending$beta,
      pending$mu,
      pending$mu,
      K_t,
      gamma_t,
      pending$foim
    )
    shared$state <- stoch_epi_step_native(
      shared$engine_ptr,
      shared$state,
      shared$t,
      shared$dt_stoch,
      1
    )
  } else {
    epi_engine_set_runtime(
      shared$engine_ptr,
      pending$beta,
      pending$mu,
      pending$mu,
      K_t,
      gamma_t,
      pending$foim
    )
    solver_step(shared$solver_ptr)
  }

  shared$t <- shared$t + 1
  shared$last_completed_timestep <- pending$timestep
  shared$pending_inputs <- native_empty_pending_inputs(shared$n_nodes)
  invisible(NULL)
}

NativeMosquitoModel <- R6::R6Class(
  "NativeMosquitoModel",
  public = list(
    shared = NULL,
    node = NULL,
    species_i = NULL,
    species_name = NULL,
    species_beta = NULL,
    cube = NULL,
    cube_info = NULL,
    initialize = function(shared, node, species_i, species_name, species_beta, cube, cube_info) {
      self$shared <- shared
      self$node <- node
      self$species_i <- species_i
      self$species_name <- species_name
      self$species_beta <- species_beta
      self$cube <- cube
      self$cube_info <- cube_info
    },
    save_state = function() {
      NULL
    },
    restore_state = function(t, state) {
      invisible(NULL)
    }
  )
)

NativeMosquitoSolver <- R6::R6Class(
  "NativeMosquitoSolver",
  public = list(
    model = NULL,
    initialize = function(model) {
      self$model <- model
    },
    get_native_state = function() {
      native_shared_state(self$model$shared)
    },
    set_native_state = function(state, t = self$model$shared$t) {
      native_set_shared_state(self$model$shared, state, t = t)
      invisible(NULL)
    },
    step = function() {
      pending <- self$model$shared$pending_inputs
      if (is.null(pending$timestep) || !isTRUE(pending$ready[[self$model$node]])) {
        stop("Native mosquito inputs were not set before stepping.")
      }
      if (isTRUE(pending$stepped[[self$model$node]])) {
        return(invisible(NULL))
      }
      pending$stepped[[self$model$node]] <- TRUE
      self$model$shared$pending_inputs <- pending
      if (all(pending$stepped)) {
        native_execute_shared_step(self$model$shared)
      }
      invisible(NULL)
    },
    get_states = function() {
      native_solver_compat_states(
        self$get_native_state(),
        self$model,
        self$model$shared$stochastic
      )
    },
    get_summary = function(node = self$model$node) {
      native_extract_summary(self$get_native_state(), self$model, node = node)
    },
    apply_release = function(sex, genotype_idx, count, timestep, node = self$model$node) {
      state <- native_apply_release_to_state(
        self$get_native_state(),
        self$model,
        sex,
        genotype_idx,
        count,
        self$model$shared$stochastic,
        node = node
      )
      self$set_native_state(state, t = self$model$shared$t)
      invisible(NULL)
    },
    save_state = function() {
      if (self$model$node != 1L) {
        return(NULL)
      }
      list(
        t = self$model$shared$t,
        state = self$get_native_state(),
        pending_inputs = self$model$shared$pending_inputs,
        last_completed_timestep = self$model$shared$last_completed_timestep
      )
    },
    restore_state = function(t, state) {
      if (is.null(state)) {
        return(invisible(NULL))
      }
      if (is.list(state) && !is.null(state$state)) {
        native_set_shared_state(self$model$shared, state$state, t = state$t)
        self$model$shared$pending_inputs <- state$pending_inputs
        self$model$shared$last_completed_timestep <- state$last_completed_timestep
      } else {
        native_set_shared_state(self$model$shared, state, t = t)
      }
      invisible(NULL)
    }
  )
)

native_build_shared_backend <- function(parameters, species_i, timesteps) {
  stopifnot(length(parameters) >= 1L)
  native_validate_human_movement(parameters)

  cfg <- native_mosquito_stage_config(parameters[[1]])
  native_assert_identical_configs(parameters, cfg)

  n_nodes <- length(parameters)
  stochastic <- isTRUE(parameters[[1]]$individual_mosquitoes)
  if (any(vapply(parameters, function(p) isTRUE(p$individual_mosquitoes) != stochastic, logical(1)))) {
    stop("All populations must agree on individual_mosquitoes when sharing the native mosquito backend.")
  }

  cube <- parameters[[1]]$cube
  if (length(parameters) > 1L &&
      any(!vapply(parameters[-1L], function(p) identical(p$cube, cube), logical(1)))) {
    stop("All populations must share the same cube for the shared native mosquito backend.")
  }
  cube_info <- cube_genotype_info(cube)

  del <- native_assert_identical_scalar(lapply(parameters, function(p) p$del), "del")
  dl <- native_assert_identical_scalar(lapply(parameters, function(p) p$dl), "dl")
  dpl <- native_assert_identical_scalar(lapply(parameters, function(p) p$dpl), "dpl")
  dem <- native_assert_identical_scalar(lapply(parameters, function(p) p$dem), "dem")
  nu <- native_assert_identical_scalar(lapply(parameters, function(p) cfg$nu), "native_mosquito_nu")
  r_tol <- native_assert_identical_scalar(lapply(parameters, function(p) p$r_tol), "r_tol")
  a_tol <- native_assert_identical_scalar(lapply(parameters, function(p) p$a_tol), "a_tol")
  ode_max_steps <- native_assert_identical_scalar(lapply(parameters, function(p) p$ode_max_steps), "ode_max_steps")

  move <- native_collect_mosquito_movement(parameters)
  index <- build_native_epi_indices(
    G = cube_info$G,
    nNodes = n_nodes,
    nE = cfg$nE,
    nL = cfg$nL,
    nP = cfg$nP,
    nEIP = cfg$nEIP
  )

  omega <- cube_omega_vector(cube, cube_info$G, cube_info$genotypesID)
  phi <- native_align_cube_vec(if (is.null(cube)) NULL else cube$phi, cube_info$genotypesID, 0.5)
  xiF <- native_align_cube_vec(if (is.null(cube)) NULL else cube$xiF, cube_info$genotypesID, 1)
  xiM <- native_align_cube_vec(if (is.null(cube)) NULL else cube$xiM, cube_info$genotypesID, 1)
  eta <- native_align_cube_mat(if (is.null(cube)) NULL else cube$eta, cube_info$genotypesID, default = 1)
  B_mat <- native_build_birth_matrix(cube, cube_info)

  muE <- vapply(parameters, function(p) p$me, numeric(1))
  muL <- vapply(parameters, function(p) p$ml, numeric(1))
  muP <- vapply(parameters, function(p) p$mup, numeric(1))
  mu0 <- vapply(parameters, function(p) p$mum[[species_i]], numeric(1))
  init_beta <- vapply(
    seq_along(parameters),
    function(node) eggs_laid(
      species_beta_value(parameters[[node]], species_i),
      mu0[[node]],
      parameters[[node]]$blood_meal_rates[[species_i]]
    ),
    numeric(1)
  )
  init_foim <- vapply(parameters, function(p) p$init_foim, numeric(1))
  time_offset <- native_assert_identical_scalar(
    lapply(
      parameters,
      function(p) {
        offset <- p$vector_control_time_offset
        if (is.null(offset)) {
          return(0L)
        }

        offset <- as.integer(offset)
        if (length(offset) != 1L || is.na(offset) || offset < 0L) {
          stop("vector_control_time_offset must be NULL or a single integer >= 0.")
        }
        offset
      }
    ),
    "vector_control_time_offset"
  )
  carrying_capacity_lookup <- native_build_carrying_capacity_lookup(parameters, species_i, timesteps)
  init_K <- native_carrying_capacity_at(carrying_capacity_lookup, 0L)

  model_type <- 1L
  # The native backend keeps human infection in the outer msimGD shell, so the
  # native human block is only used as a carrier for mosquito-side FOIM.
  b_vec <- rep(0, cube_info$G)
  c_vec <- native_align_cube_vec(if (is.null(cube)) NULL else cube$c, cube_info$genotypesID, 1)

  # --- Validate cube$b / cube$c / vector_infectivity_g consistency ---
  has_b <- !is.null(cube) && !is.null(cube$b)
  has_c <- !is.null(cube) && !is.null(cube$c)
  has_vi <- !is.null(parameters[[1]]$vector_infectivity_g)
  if (has_vi && has_b) {
    message("Note: vector_infectivity_g is set and will override cube$b for human-side transmission weighting.")
  }
  if (has_b && !has_c) {
    warning("cube$b is set but cube$c is NULL; c_vec defaults to 1 (all genotypes equally ",
            "susceptible). Set cube$c explicitly if you want genotype-specific mosquito infection.",
            call. = FALSE)
  }
  if (has_c && !has_b && !has_vi) {
    warning("cube$c is set but cube$b is NULL and vector_infectivity_g is not set; ",
            "all genotypes will transmit equally to humans. Set cube$b or vector_infectivity_g ",
            "if you want genotype-specific blocking.",
            call. = FALSE)
  }

  engine_ptr <- if (stochastic) {
    stoch_epi_engine_create(
      nNodes = n_nodes,
      nM = n_nodes,
      nG = cube_info$G,
      nE = cfg$nE,
      nL = cfg$nL,
      nP = cfg$nP,
      nEIP = cfg$nEIP,
      model_type = model_type,
      egg_ix_r = as.integer(index$egg_ix),
      larv_ix_r = as.integer(index$larv_ix),
      pup_ix_r = as.integer(index$pup_ix),
      male_ix_r = as.integer(index$male_ix),
      unm_ix_r = as.integer(index$unm_ix),
      fem_ix_r = as.integer(index$fem_ix),
      rE = cfg$nE / del,
      rL = cfg$nL / dl,
      rP = cfg$nP / dpl,
      rEIP = cfg$nEIP / dem,
      muE_r = muE,
      muL_r = muL,
      muP_r = muP,
      muM_r = mu0,
      muF_r = mu0,
      log_dd = TRUE,
      K_r = init_K,
      gamma_dd_r = rep(0, n_nodes),
      beta = init_beta[[1L]],
      nu = nu,
      omega_inv_r = ifelse(omega == 0, 1e3, 1 / omega),
      phi_r = phi,
      xiF_r = xiF,
      xiM_r = xiM,
      eta_r = eta,
      B_mat_r = B_mat,
      has_move = move$has_move,
      move_probs_r = move$move_probs,
      move_rates_r = move$move_rates,
      muM_node_base_r = rep(1, n_nodes),
      muF_node_base_r = rep(1, n_nodes),
      tol = 1e-8,
      c_vec_r = c_vec,
      b_vec_r = b_vec,
      a_vec_r = init_foim,
      muH_param = 0,
      r_param = 0,
      delta_param = 0,
      cT_vec_r = NULL,
      cD_vec_r = NULL,
      cU_vec_r = NULL,
      W_age_r = NULL,
      d1 = 0,
      fd_r = numeric(0),
      ID0 = 1,
      kd = 1,
      gamma1_imp = 0,
      mosy_nodes_r = seq_len(n_nodes),
      human_nodes_r = seq_len(n_nodes),
      hS_ix_r = as.integer(index$hS_ix),
      hI_ix_r = as.integer(index$hI_ix),
      hE_ix_r = as.integer(NA_integer_),
      hR_ix_r = as.integer(NA_integer_),
      has_hmove = FALSE,
      h_move_probs_r = matrix(0, nrow = n_nodes, ncol = n_nodes),
      h_move_rates_r = rep(0, n_nodes),
      nState = index$total_state_len
    )
  } else {
    epi_engine_create(
      nNodes = n_nodes,
      nM = n_nodes,
      nG = cube_info$G,
      nE = cfg$nE,
      nL = cfg$nL,
      nP = cfg$nP,
      nEIP = cfg$nEIP,
      model_type = model_type,
      egg_ix = index$egg_ix,
      larv_ix = index$larv_ix,
      pup_ix = index$pup_ix,
      male_ix = index$male_ix,
      unm_ix = index$unm_ix,
      fem_ix_flat = as.integer(index$fem_ix),
      rE = cfg$nE / del,
      rL = cfg$nL / dl,
      rP = cfg$nP / dpl,
      rEIP = cfg$nEIP / dem,
      muE = muE,
      muL = muL,
      muP = muP,
      muM = mu0,
      muF = mu0,
      log_dd = TRUE,
      K = init_K,
      gamma_dd = rep(0, n_nodes),
      beta = init_beta[[1L]],
      nu = nu,
      omega = omega,
      phi_cube = phi,
      xiF = xiF,
      xiM = xiM,
      eta = eta,
      B_mat = B_mat,
      has_move = move$has_move,
      move_probs_dense = move$move_probs,
      move_rates = move$move_rates,
      muM_node_base = rep(1, n_nodes),
      muF_node_base = rep(1, n_nodes),
      tol = 1e-8,
      c_vec_r = c_vec,
      b_vec_r = b_vec,
      a_vec = init_foim,
      muH_param = 0,
      r_param = 0,
      delta_param = 0,
      cT_vec_r = NULL,
      cD_vec_r = NULL,
      cU_vec_r = NULL,
      W_age_r = NULL,
      d1 = 0,
      fd_r = numeric(0),
      ID0 = 1,
      kd = 1,
      gamma1_imp = 0,
      mosy_nodes_r = seq_len(n_nodes),
      human_nodes_r = seq_len(n_nodes),
      hS_ix_r = as.integer(index$hS_ix),
      hI_ix_r = as.integer(index$hI_ix),
      hE_ix_r = as.integer(NA_integer_),
      hR_ix_r = as.integer(NA_integer_),
      has_hmove = FALSE,
      h_move_probs_dense = matrix(0, nrow = n_nodes, ncol = n_nodes),
      h_move_rates = rep(0, n_nodes),
      nState = index$total_state_len
    )
  }

  init_state <- native_initial_state(parameters, species_i, cube_info, index)
  if (stochastic) {
    stoch_epi_engine_set_runtime(
      engine_ptr,
      init_beta,
      mu0,
      mu0,
      init_K,
      rep(0, n_nodes),
      init_foim
    )
    solver_ptr <- NULL
  } else {
    epi_engine_set_runtime(
      engine_ptr,
      init_beta,
      mu0,
      mu0,
      init_K,
      rep(0, n_nodes),
      init_foim
    )
    solver_ptr <- create_epi_solver(
      engine_ptr,
      init_state,
      r_tol,
      a_tol,
      as.integer(ode_max_steps)
    )
  }

  shared <- new.env(parent = emptyenv())
  shared$engine_ptr <- engine_ptr
  shared$solver_ptr <- solver_ptr
  shared$stochastic <- stochastic
  shared$dt_stoch <- cfg$dt_stoch
  shared$state <- if (stochastic) init_state else NULL
  shared$index <- index
  shared$eta <- eta
  shared$carrying_capacity_lookup <- carrying_capacity_lookup
  shared$gamma_dd <- rep(0, n_nodes)
  shared$n_nodes <- n_nodes
  shared$time_offset <- as.integer(time_offset)
  shared$t <- 0
  shared$last_completed_timestep <- NULL
  shared$pending_inputs <- native_empty_pending_inputs(n_nodes)

  models <- lapply(
    seq_len(n_nodes),
    function(node) {
      NativeMosquitoModel$new(
        shared = shared,
        node = node,
        species_i = species_i,
        species_name = parameters[[node]]$species[[species_i]],
        species_beta = species_beta_value(parameters[[node]], species_i),
        cube = cube,
        cube_info = cube_info
      )
    }
  )
  solvers <- lapply(models, NativeMosquitoSolver$new)

  list(models = models, solvers = solvers, shared = shared)
}

create_native_mosquito_model <- function(parameters, species_i, timesteps) {
  native_build_shared_backend(list(parameters), species_i, timesteps)$models[[1L]]
}

create_native_mosquito_solver <- function(model, parameters) {
  NativeMosquitoSolver$new(model)
}

parameterise_native_metapop_backends <- function(parameters, timesteps) {
  native_warn_tau_leap_backend(parameters[[1]])
  n_pops <- length(parameters)
  models <- vector("list", n_pops)
  solvers <- vector("list", n_pops)
  for (i in seq_len(n_pops)) {
    models[[i]] <- vector("list", length(parameters[[i]]$species))
    solvers[[i]] <- vector("list", length(parameters[[i]]$species))
  }

  n_species <- length(parameters[[1]]$species)
  for (species_i in seq_len(n_species)) {
    backend <- native_build_shared_backend(parameters, species_i, timesteps)
    for (pop_i in seq_len(n_pops)) {
      models[[pop_i]][[species_i]] <- backend$models[[pop_i]]
      solvers[[pop_i]][[species_i]] <- backend$solvers[[pop_i]]
    }
  }

  list(models = models, solvers = solvers)
}

native_interleave_metapop_processes <- function(processes) {
  if (length(processes) <= 1L) {
    return(unlist(processes, recursive = FALSE))
  }

  reference_names <- names(processes[[1]])
  reference_len <- length(processes[[1]])
  for (i in seq_along(processes)) {
    if (length(processes[[i]]) != reference_len || !identical(names(processes[[i]]), reference_names)) {
      stop(
        paste(
          "Native metapop mosquito stepping requires identical process layouts across populations.",
          "The process lists differ, so synchronized shared-node stepping would be ambiguous."
        )
      )
    }
  }

  out <- vector("list", reference_len * length(processes))
  out_names <- character(length(out))
  idx <- 1L
  for (step_i in seq_len(reference_len)) {
    for (pop_i in seq_along(processes)) {
      out[[idx]] <- processes[[pop_i]][[step_i]]
      out_names[[idx]] <- sprintf("%s_pop%d", reference_names[[step_i]], pop_i)
      idx <- idx + 1L
    }
  }
  names(out) <- out_names
  out
}

native_mosquito_model_update <- function(model, timestep, mu, foim, f) {
  shared <- model$shared
  pending <- shared$pending_inputs
  time_offset <- if (!is.null(shared$time_offset)) as.integer(shared$time_offset) else 0L
  absolute_timestep <- as.integer(timestep) + time_offset
  if (!is.null(pending$timestep) && pending$timestep != absolute_timestep) {
    stop("Native mosquito backend received mixed timesteps before stepping.")
  }
  foim <- as.numeric(foim)
  if (length(foim) != 1L || anyNA(foim) || !is.finite(foim) || foim < 0) {
    stop("Native mosquito backend FOIM must be a single nonnegative finite value.", call. = FALSE)
  }

  pending$timestep <- absolute_timestep
  pending$beta[[model$node]] <- as.numeric(eggs_laid(model$species_beta, mu, f))
  pending$mu[[model$node]] <- as.numeric(mu)
  pending$foim[[model$node]] <- foim[[1L]]
  pending$ready[[model$node]] <- TRUE
  shared$pending_inputs <- pending
  invisible(NULL)
}
