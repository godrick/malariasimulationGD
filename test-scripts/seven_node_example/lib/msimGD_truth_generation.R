# ------------------------------------------------------------------------------
# msimGD_truth_generation.R
# ------------------------------------------------------------------------------
# Generate stochastic truth from msimGD metapop (individual-based epi) and
# extract infection/clinical incidence plus carrier frequency for downstream
# calibration and inference.
#
# Main functions:
#   1. msimGD_run_truth()            — run one stochastic metapop simulation
#   2. msimGD_extract_epi_obs()      — aggregate infection or clinical incidence
#                                      into windows
#   3. msimGD_extract_prevalence_obs() — extract age-specific prevalence
#   4. msimGD_extract_carrier_freq() — carrier freq from genotype counts
#
# Dependencies (must be loaded before sourcing):
#   - msimGD package    via devtools::load_all("scripts/msimGD")
#   - movement_mu.R     for mosquito_movement_from_mu()
#   - customMGDrive2    for calc_move_rate()
# ------------------------------------------------------------------------------


#' Build a base msimGD parameter list with mosquito biology overridden to match
#' the customMGDrive2 lifecycle parameters used by the deterministic forward model.
#'
#' @param NH human population per node
#' @param theta lifecycle parameter list from customMGDrive2
#' @param individual_mosquitoes whether to use the native stochastic mosquito backend
#' @return msimGD parameter list before equilibrium calibration / movement
msimGD_make_base_parameters <- function(NH, theta, individual_mosquitoes = TRUE) {
  bp <- get_parameters(list(
    human_population = NH,
    individual_mosquitoes = individual_mosquitoes,
    native_mosquito_backend = TRUE,
    progress_bar = FALSE
  ))

  # Match the mosquito life-history parameters used by the inference forward model.
  bp$del <- theta$nE / theta$qE
  bp$dl  <- theta$nL / theta$qL
  bp$dpl <- theta$nP / theta$qP
  bp$me  <- theta$muE
  bp$ml  <- theta$muL
  bp$mup <- theta$muP
  bp$mum <- theta$muF
  bp$beta <- theta$beta
  bp$native_mosquito_nE <- theta$nE
  bp$native_mosquito_nL <- theta$nL
  bp$native_mosquito_nP <- theta$nP
  if (!is.null(theta$nu)) {
    bp$native_mosquito_nu <- theta$nu
  }

  bp
}


#' Calibrate msimGD's equilibrium mosquito abundance from the requested EIR.
#'
#' @param NH human population per node
#' @param theta lifecycle parameter list from customMGDrive2
#' @param init_EIR endemic EIR passed into set_equilibrium()
#' @param individual_mosquitoes whether to build stochastic or deterministic-style parameters
#' @param parameter_modifier optional function `(parameters, node_index,
#'   warmup_days)` applied before equilibrium is solved
#' @param node_index node index passed into `parameter_modifier`
#' @param warmup_days warmup days passed into `parameter_modifier`
#' @return list(total_M, init_foim, parameters) after equilibrium calibration
msimGD_equilibrium_baseline <- function(NH, theta, init_EIR = 5,
                                        individual_mosquitoes = TRUE,
                                        parameter_modifier = NULL,
                                        node_index = 1L,
                                        warmup_days = 0L) {
  bp <- msimGD_make_base_parameters(
    NH = NH,
    theta = theta,
    individual_mosquitoes = individual_mosquitoes
  )
  if (!is.null(parameter_modifier)) {
    bp <- parameter_modifier(
      bp,
      node_index = as.integer(node_index),
      warmup_days = as.integer(warmup_days)
    )
    if (is.null(bp) || !is.list(bp)) {
      stop("parameter_modifier must return a parameter list.", call. = FALSE)
    }
  }
  bp <- set_equilibrium(bp, init_EIR = init_EIR)

  list(
    total_M = unname(as.numeric(bp$total_M)),
    init_foim = unname(as.numeric(bp$init_foim)),
    parameters = bp
  )
}


msimGD_stationary_human_library_snapshot_timesteps <- function(burnin_timesteps,
                                                               n_snapshots = 1L,
                                                               snapshot_spacing = 0L) {
  malariasimulationGD:::stationary_initializer_snapshot_timesteps(
    burnin_timesteps = burnin_timesteps,
    n_snapshots = n_snapshots,
    snapshot_spacing = snapshot_spacing,
    label = "Requested stationary human library snapshots"
  )
}


msimGD_snapshot_total_M_by_node <- function(result) {
  if (!is.list(result) || length(result) < 1L) {
    stop("result must be a non-empty list of node-level outputs.", call. = FALSE)
  }

  vapply(
    result,
    function(df) {
      total_cols <- grep("^total_M_", names(df), value = TRUE)
      if (length(total_cols) < 1L) {
        stop("result does not contain any exported `total_M_*` columns.", call. = FALSE)
      }
      sum(as.numeric(df[nrow(df), total_cols, drop = TRUE]))
    },
    numeric(1)
  )
}

msimGD_simulation_total_M_by_node <- function(sim) {
  if (is.list(sim) &&
      !is.null(sim$summary) &&
      is.list(sim$summary) &&
      !is.null(sim$summary$total_M_by_node)) {
    return(as.numeric(sim$summary$total_M_by_node))
  }

  msimGD_snapshot_total_M_by_node(sim$data)
}

msimGD_snapshot_init_foim_by_node <- function(snapshot, parameters) {
  if (!is.list(snapshot) || is.null(snapshot$node_libraries)) {
    stop("snapshot must be a metapop stationary snapshot.", call. = FALSE)
  }
  if (!is.list(parameters) || length(parameters) != length(snapshot$node_libraries)) {
    stop("parameters must match the snapshot node libraries.", call. = FALSE)
  }

  vapply(
    seq_along(snapshot$node_libraries),
    function(nd) {
      foim <- malariasimulationGD:::stationary_human_initializer_library_foim(
        snapshot$node_libraries[[nd]],
        parameters[[nd]]
      )
      if (!is.numeric(foim) || length(foim) != 1L || !is.finite(foim) || foim < 0) {
        foim <- as.numeric(parameters[[nd]]$init_foim)
      }
      if (!is.numeric(foim) || length(foim) != 1L || !is.finite(foim) || foim < 0) {
        foim <- NA_real_
      }
      foim
    },
    numeric(1)
  )
}


msimGD_build_stationary_checkpoint_library <- function(
    setup,
    cube,
    NF,
    NH,
    mu,
    p_move,
    theta,
    init_EIR = 5,
    parameter_modifier = NULL,
    burnin_timesteps = 5840L,
    n_snapshots = 1L,
    snapshot_spacing = 0L,
    seed = NULL,
    library_path = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  snapshot_timesteps <- malariasimulationGD:::stationary_initializer_snapshot_timesteps(
    burnin_timesteps = burnin_timesteps,
    n_snapshots = n_snapshots,
    snapshot_spacing = snapshot_spacing,
    label = "Requested stationary checkpoint snapshots"
  )
  prep <- .msimGD_prepare_truth_run(
    setup = setup,
    cube = cube,
    NF = NF,
    NH = NH,
    tmax = 1L,
    mu = mu,
    p_move = p_move,
    release = NULL,
    theta = theta,
    init_EIR = init_EIR,
    prevalence_band = list(min = integer(0), max = integer(0)),
    infection_band = list(min = integer(0), max = integer(0)),
    clinical_band = list(min = integer(0), max = integer(0)),
    warmup_days = 0L,
    parameter_modifier = parameter_modifier,
    baseline_checkpoint_state = NULL,
    baseline_checkpoint_metadata = NULL,
    human_initialization_library = NULL
  )

  current_state <- NULL
  checkpoints <- vector("list", length(snapshot_timesteps))
  snapshot_total_M_by_node <- vector("list", length(snapshot_timesteps))
  for (i in seq_along(snapshot_timesteps)) {
    message(sprintf(
      "[checkpoint library] snapshot %d/%d | day=%d",
      i,
      length(snapshot_timesteps),
      snapshot_timesteps[[i]]
    ))
    sim <- run_metapop_simulation(
      timesteps = snapshot_timesteps[[i]],
      parameters = prep$parameters,
      mixing_tt = 1,
      export_mixing = list(diag(prep$n_nodes)),
      import_mixing = list(diag(prep$n_nodes)),
      p_captured_tt = 1,
      p_captured = list(matrix(0, prep$n_nodes, prep$n_nodes)),
      p_success = 0,
      initial_state = current_state,
      restore_random_state = !is.null(current_state),
      return_state = TRUE,
      render_output = FALSE,
      return_summary = TRUE
    )
    current_state <- sim$state
    checkpoints[[i]] <- list(
      state = current_state,
      metadata = list(
        snapshot_index = as.integer(i),
        timesteps = as.integer(snapshot_timesteps[[i]]),
        seasonal_cycle_days = 365L,
        seasonal_phase_day = .msimGD_seasonal_phase_day(snapshot_timesteps[[i]]),
        baseline_time_dependent_signature = prep$baseline_time_dependent_signature,
        baseline_contact_signature = prep$baseline_contact_signature,
        contact_multiplier_by_node = prep$contact_multiplier_by_node,
        contact_covariates_by_node = prep$contact_covariates_by_node
      )
    )
    snapshot_total_M_by_node[[i]] <- msimGD_simulation_total_M_by_node(sim)
  }

  snapshot_total_M_by_node <- do.call(rbind, snapshot_total_M_by_node)
  colnames(snapshot_total_M_by_node) <- seq_len(ncol(snapshot_total_M_by_node))
  library <- malariasimulationGD:::new_stationary_checkpoint_library(
    checkpoints = checkpoints,
    metadata = list(
      checkpoint_type = "metapopulation_stationary_library",
      build_mode = "metapopulation",
      burnin_timesteps = as.integer(burnin_timesteps),
      n_snapshots = as.integer(n_snapshots),
      snapshot_spacing = as.integer(snapshot_spacing),
      snapshot_timesteps = as.integer(snapshot_timesteps),
      snapshot_phase_days = vapply(
        snapshot_timesteps,
        .msimGD_seasonal_phase_day,
        integer(1)
      ),
      seed = seed,
      init_EIR = as.numeric(init_EIR),
      seasonal_cycle_days = 365L,
      baseline_time_dependent_signature = prep$baseline_time_dependent_signature,
      baseline_contact_signature = prep$baseline_contact_signature,
      contact_multiplier_by_node = prep$contact_multiplier_by_node,
      contact_covariates_by_node = prep$contact_covariates_by_node,
      equilibrium_total_M = prep$equilibrium_total_M,
      effective_total_M = prep$effective_total_M,
      init_foim = prep$init_foim,
      parameter_total_M = prep$effective_total_M,
      snapshot_total_M_by_node = snapshot_total_M_by_node
    )
  )

  if (!is.null(library_path)) {
    saveRDS(library, library_path)
  }

  library
}


msimGD_build_node_conditioned_human_library <- function(
    setup,
    cube,
    NH,
    mu,
    p_move,
    theta,
    init_EIR = 5,
    parameter_modifier = NULL,
    burnin_timesteps = 5840L,
    n_snapshots = 3L,
    snapshot_spacing = 730L,
    seed = NULL,
    library_path = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  snapshot_timesteps <- msimGD_stationary_human_library_snapshot_timesteps(
    burnin_timesteps = burnin_timesteps,
    n_snapshots = n_snapshots,
    snapshot_spacing = snapshot_spacing
  )
  prep <- .msimGD_prepare_truth_run(
    setup = setup,
    cube = cube,
    NF = NULL,
    NH = NH,
    tmax = 1L,
    mu = mu,
    p_move = p_move,
    release = NULL,
    theta = theta,
    init_EIR = init_EIR,
    prevalence_band = list(min = integer(0), max = integer(0)),
    infection_band = list(min = integer(0), max = integer(0)),
    clinical_band = list(min = integer(0), max = integer(0)),
    warmup_days = 0L,
    parameter_modifier = parameter_modifier,
    baseline_checkpoint_state = NULL,
    human_initialization_library = NULL
  )

  current_state <- NULL
  snapshots <- vector("list", length(snapshot_timesteps))
  snapshot_total_M_by_node <- vector("list", length(snapshot_timesteps))
  for (i in seq_along(snapshot_timesteps)) {
    message(sprintf(
      "[human library] snapshot %d/%d | day=%d",
      i,
      length(snapshot_timesteps),
      snapshot_timesteps[[i]]
    ))
    sim <- run_metapop_simulation(
      timesteps = snapshot_timesteps[[i]],
      parameters = prep$parameters,
      mixing_tt = 1,
      export_mixing = list(diag(prep$n_nodes)),
      import_mixing = list(diag(prep$n_nodes)),
      p_captured_tt = 1,
      p_captured = list(matrix(0, prep$n_nodes, prep$n_nodes)),
      p_success = 0,
      initial_state = current_state,
      restore_random_state = !is.null(current_state),
      return_state = TRUE,
      render_output = FALSE,
      return_summary = TRUE
    )
    current_state <- sim$state
    snapshots[[i]] <- malariasimulationGD:::stationary_human_initializer_extract_metapop_library(
      current_state,
      prep$parameters
    )
    snapshots[[i]]$metadata$snapshot_index <- as.integer(i)
    snapshots[[i]]$metadata$snapshot_timestep <- as.integer(snapshot_timesteps[[i]])
    snapshot_total_M_by_node[[i]] <- msimGD_simulation_total_M_by_node(sim)
  }

  library <- malariasimulationGD:::stationary_human_initializer_bind_node_conditioned_libraries(snapshots)
  snapshot_total_M_by_node <- do.call(rbind, snapshot_total_M_by_node)
  colnames(snapshot_total_M_by_node) <- seq_len(ncol(snapshot_total_M_by_node))
  baseline_total_M_by_node <- colMeans(snapshot_total_M_by_node)
  baseline_init_foim_by_node <- vapply(
    seq_along(library$node_libraries),
    function(nd) {
      malariasimulationGD:::stationary_human_initializer_library_foim(
        library$node_libraries[[nd]],
        prep$parameters[[nd]]
      )
    },
    numeric(1)
  )
  library$metadata <- c(
    list(
      build_mode = "metapopulation",
      burnin_timesteps = as.integer(burnin_timesteps),
      n_snapshots = as.integer(n_snapshots),
      snapshot_spacing = as.integer(snapshot_spacing),
      snapshot_timesteps = as.integer(snapshot_timesteps),
      seed = seed,
      stationary_total_M_by_node = baseline_total_M_by_node,
      init_EIR = as.numeric(init_EIR),
      baseline_time_dependent_signature = prep$baseline_time_dependent_signature,
      baseline_contact_signature = prep$baseline_contact_signature,
      contact_multiplier_by_node = prep$contact_multiplier_by_node,
      contact_covariates_by_node = prep$contact_covariates_by_node,
      stationary_init_foim_by_node = baseline_init_foim_by_node,
      baseline_total_M_by_node = baseline_total_M_by_node,
      baseline_init_foim_by_node = baseline_init_foim_by_node,
      snapshot_total_M_by_node = snapshot_total_M_by_node,
      equilibrium_total_M = prep$equilibrium_total_M,
      effective_total_M = prep$effective_total_M,
      init_foim = prep$init_foim
    ),
    library$metadata
  )

  if (!is.null(library_path)) {
    saveRDS(library, library_path)
  }

  library
}

msimGD_build_metapop_stationary_human_library <- function(
    setup,
    cube,
    NH,
    mu,
    p_move,
    theta,
    init_EIR = 5,
    parameter_modifier = NULL,
    burnin_timesteps = 5840L,
    n_snapshots = 3L,
    snapshot_spacing = 730L,
    seed = NULL,
    library_path = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  snapshot_timesteps <- msimGD_stationary_human_library_snapshot_timesteps(
    burnin_timesteps = burnin_timesteps,
    n_snapshots = n_snapshots,
    snapshot_spacing = snapshot_spacing
  )
  prep <- .msimGD_prepare_truth_run(
    setup = setup,
    cube = cube,
    NF = NULL,
    NH = NH,
    tmax = 1L,
    mu = mu,
    p_move = p_move,
    release = NULL,
    theta = theta,
    init_EIR = init_EIR,
    prevalence_band = list(min = integer(0), max = integer(0)),
    infection_band = list(min = integer(0), max = integer(0)),
    clinical_band = list(min = integer(0), max = integer(0)),
    warmup_days = 0L,
    parameter_modifier = parameter_modifier,
    baseline_checkpoint_state = NULL,
    human_initialization_library = NULL
  )

  current_state <- NULL
  snapshots <- vector("list", length(snapshot_timesteps))
  snapshot_total_M_by_node <- vector("list", length(snapshot_timesteps))
  snapshot_init_foim_by_node <- vector("list", length(snapshot_timesteps))
  for (i in seq_along(snapshot_timesteps)) {
    message(sprintf(
      "[metapop stationary library] snapshot %d/%d | day=%d",
      i,
      length(snapshot_timesteps),
      snapshot_timesteps[[i]]
    ))
    sim <- run_metapop_simulation(
      timesteps = snapshot_timesteps[[i]],
      parameters = prep$parameters,
      mixing_tt = 1,
      export_mixing = list(diag(prep$n_nodes)),
      import_mixing = list(diag(prep$n_nodes)),
      p_captured_tt = 1,
      p_captured = list(matrix(0, prep$n_nodes, prep$n_nodes)),
      p_success = 0,
      initial_state = current_state,
      restore_random_state = !is.null(current_state),
      return_state = TRUE,
      render_output = FALSE,
      return_summary = TRUE
    )
    current_state <- sim$state
    snapshots[[i]] <- malariasimulationGD:::stationary_human_initializer_extract_metapop_snapshot(
      current_state,
      prep$parameters
    )
    snapshots[[i]]$metadata$snapshot_index <- as.integer(i)
    snapshots[[i]]$metadata$snapshot_timestep <- as.integer(snapshot_timesteps[[i]])
    snapshot_total_M_by_node[[i]] <- msimGD_simulation_total_M_by_node(sim)
    snapshot_init_foim_by_node[[i]] <- msimGD_snapshot_init_foim_by_node(
      snapshots[[i]],
      prep$parameters
    )
    snapshots[[i]]$metadata$stationary_total_M_by_node <- snapshot_total_M_by_node[[i]]
    snapshots[[i]]$metadata$stationary_init_foim_by_node <- snapshot_init_foim_by_node[[i]]
  }

  library <- malariasimulationGD:::stationary_human_initializer_bind_metapop_snapshots(snapshots)
  snapshot_total_M_by_node <- do.call(rbind, snapshot_total_M_by_node)
  snapshot_init_foim_by_node <- do.call(rbind, snapshot_init_foim_by_node)
  colnames(snapshot_total_M_by_node) <- seq_len(ncol(snapshot_total_M_by_node))
  colnames(snapshot_init_foim_by_node) <- seq_len(ncol(snapshot_init_foim_by_node))
  baseline_total_M_by_node <- colMeans(snapshot_total_M_by_node)
  baseline_init_foim_by_node <- colMeans(snapshot_init_foim_by_node)
  library$metadata <- c(
    list(
      build_mode = "metapopulation",
      library_kind = "full_metapop_snapshots",
      burnin_timesteps = as.integer(burnin_timesteps),
      n_snapshots = as.integer(n_snapshots),
      snapshot_spacing = as.integer(snapshot_spacing),
      snapshot_timesteps = as.integer(snapshot_timesteps),
      seed = seed,
      stationary_total_M_by_node = baseline_total_M_by_node,
      init_EIR = as.numeric(init_EIR),
      baseline_time_dependent_signature = prep$baseline_time_dependent_signature,
      baseline_contact_signature = prep$baseline_contact_signature,
      contact_multiplier_by_node = prep$contact_multiplier_by_node,
      contact_covariates_by_node = prep$contact_covariates_by_node,
      stationary_init_foim_by_node = baseline_init_foim_by_node,
      baseline_total_M_by_node = baseline_total_M_by_node,
      baseline_init_foim_by_node = baseline_init_foim_by_node,
      snapshot_total_M_by_node = snapshot_total_M_by_node,
      snapshot_init_foim_by_node = snapshot_init_foim_by_node,
      equilibrium_total_M = prep$equilibrium_total_M,
      effective_total_M = prep$effective_total_M,
      init_foim = prep$init_foim
    ),
    library$metadata
  )

  if (!is.null(library_path)) {
    saveRDS(library, library_path)
  }

  library
}


.msimGD_normalise_optional_age_band <- function(min_age, max_age, name) {
  if (is.null(min_age) && is.null(max_age)) {
    return(list(min = integer(0), max = integer(0)))
  }
  if (is.null(min_age) || is.null(max_age)) {
    stop(sprintf(
      "%s age bounds must either both be NULL or both be provided.",
      name
    ), call. = FALSE)
  }
  if (length(min_age) != length(max_age) || length(min_age) < 1L) {
    stop(sprintf(
      "%s age bounds must be equal-length day-value vectors.",
      name
    ), call. = FALSE)
  }

  min_age <- as.integer(min_age)
  max_age <- as.integer(max_age)
  if (anyNA(min_age) || anyNA(max_age) || any(min_age < 0L) || any(max_age < min_age)) {
    stop(sprintf(
      "%s age bounds must satisfy 0 <= min_age <= max_age.",
      name
    ), call. = FALSE)
  }

  list(min = min_age, max = max_age)
}


.msimGD_release_is_active <- function(release) {
  if (is.null(release)) {
    return(FALSE)
  }
  nodes <- if (!is.null(release$nodes)) release$nodes else integer(0)
  size <- if (!is.null(release$size)) release$size else 0L
  length(nodes) > 0L && as.integer(size) > 0L
}


.msimGD_find_epi_column <- function(first_df, endpoint, age_min = NULL, age_max = NULL) {
  endpoint <- match.arg(endpoint, c("infection", "clinical"))
  col_pattern <- switch(
    endpoint,
    infection = "^n_inc_[0-9]+_[0-9]+$",
    clinical = "^n_inc_clinical_[0-9]+_[0-9]+$"
  )
  col_prefix <- switch(
    endpoint,
    infection = "n_inc_",
    clinical = "n_inc_clinical_"
  )

  epi_cols <- grep(col_pattern, names(first_df), value = TRUE, perl = TRUE)
  if (length(epi_cols) == 0L) {
    stop(sprintf(
      "No %s incidence column found. Set the matching rendering ages in parameters.",
      endpoint
    ), call. = FALSE)
  }

  if (!is.null(age_min) || !is.null(age_max)) {
    if (is.null(age_min) || is.null(age_max)) {
      stop("age_min and age_max must both be supplied when selecting a specific endpoint column.",
           call. = FALSE)
    }
    target_suffix <- sprintf("%d_%d", as.integer(age_min), as.integer(age_max))
    target_col <- paste0(col_prefix, target_suffix)
    if (!(target_col %in% epi_cols)) {
      stop(sprintf(
        "Requested %s incidence column `%s` was not found.",
        endpoint, target_col
      ), call. = FALSE)
    }
    epi_col <- target_col
  } else {
    epi_col <- epi_cols[1]
  }

  age_suffix <- sub(paste0("^", col_prefix), "", epi_col)
  age_col <- paste0("n_age_", age_suffix)

  list(
    endpoint = endpoint,
    epi_col = epi_col,
    age_col = age_col,
    age_suffix = age_suffix
  )
}


msimGD_epi_rate_ppy <- function(epi_obs) {
  if (!all(c("Y", "offset_person_days") %in% names(epi_obs))) {
    stop("epi_obs must contain `Y` and `offset_person_days` columns.", call. = FALSE)
  }

  total_person_days <- sum(epi_obs$offset_person_days)
  if (!is.finite(total_person_days) || total_person_days <= 0) {
    stop("Total person-days must be positive to compute a rate.", call. = FALSE)
  }

  365 * sum(epi_obs$Y) / total_person_days
}


.msimGD_find_prevalence_columns <- function(first_df, method = c("lm", "pcr"),
                                            age_min, age_max) {
  method <- match.arg(method)
  age_suffix <- sprintf("%d_%d", as.integer(age_min), as.integer(age_max))
  detect_col <- sprintf("n_detect_%s_%s", method, age_suffix)
  age_col <- sprintf("n_age_%s", age_suffix)

  if (!(detect_col %in% names(first_df))) {
    stop(sprintf(
      "Requested prevalence column `%s` was not found.",
      detect_col
    ), call. = FALSE)
  }
  if (!(age_col %in% names(first_df))) {
    stop(sprintf(
      "Requested age-band population column `%s` was not found.",
      age_col
    ), call. = FALSE)
  }

  list(
    method = method,
    detect_col = detect_col,
    age_col = age_col,
    age_suffix = age_suffix
  )
}


.msimGD_crop_result <- function(result, warmup_days) {
  if (is.null(warmup_days) || warmup_days <= 0L) {
    return(result)
  }

  warmup_days <- as.integer(warmup_days)
  lapply(result, function(df) {
    keep <- seq.int(warmup_days + 1L, nrow(df))
    cropped <- df[keep, , drop = FALSE]
    if ("timestep" %in% names(cropped)) {
      cropped$timestep <- seq_len(nrow(cropped))
    }

    attr_names <- names(attributes(df))
    preserve_attrs <- setdiff(attr_names, c("names", "row.names", "class"))
    for (attr_name in preserve_attrs) {
      attr_val <- attr(df, attr_name)
      if (is.matrix(attr_val) && nrow(attr_val) == nrow(df)) {
        attr(cropped, attr_name) <- attr_val[keep, , drop = FALSE]
      } else if (is.numeric(attr_val) && length(attr_val) == nrow(df)) {
        attr(cropped, attr_name) <- attr_val[keep]
      } else {
        attr(cropped, attr_name) <- attr_val
      }
    }

    release_schedule <- attr(cropped, "mosquito_release_schedule")
    if (!is.null(release_schedule) && "timestep" %in% names(release_schedule)) {
      keep_release <- release_schedule$timestep > warmup_days
      release_schedule <- release_schedule[keep_release, , drop = FALSE]
      release_schedule$timestep <- release_schedule$timestep - warmup_days
      attr(cropped, "mosquito_release_schedule") <- release_schedule
    }

    cropped
  })
}


.msimGD_normalise_checkpoint <- function(baseline_checkpoint) {
  if (is.null(baseline_checkpoint)) {
    return(NULL)
  }

  checkpoint <- baseline_checkpoint
  if (is.character(checkpoint)) {
    if (length(checkpoint) != 1L) {
      stop("baseline_checkpoint path must be a single file path.", call. = FALSE)
    }
    if (!file.exists(checkpoint)) {
      stop(sprintf("Baseline checkpoint file not found: %s", checkpoint), call. = FALSE)
    }
    checkpoint <- readRDS(checkpoint)
  }

  state <- checkpoint
  if (is.list(checkpoint) && !is.null(checkpoint$state)) {
    state <- checkpoint$state
  }

  if (!is.list(state) ||
      is.null(state$timesteps) ||
      is.null(state$individual) ||
      is.null(state$malariasimulationGD)) {
    stop(
      paste(
        "baseline_checkpoint must be either a resumable state object or a list",
        "with a `$state` entry containing `timesteps`, `individual`, and",
        "`malariasimulationGD`."
      ),
      call. = FALSE
    )
  }

  list(
    state = state,
    metadata = if (is.list(checkpoint) && !is.null(checkpoint$metadata)) checkpoint$metadata else list()
  )
}


.msimGD_seasonal_phase_day <- function(timesteps, cycle_days = 365L) {
  timesteps <- as.integer(timesteps)
  cycle_days <- as.integer(cycle_days)
  if (length(timesteps) != 1L || is.na(timesteps) || timesteps < 0L) {
    stop("timesteps must be a single integer >= 0.", call. = FALSE)
  }
  if (length(cycle_days) != 1L || is.na(cycle_days) || cycle_days <= 0L) {
    stop("cycle_days must be a single integer > 0.", call. = FALSE)
  }

  as.integer(timesteps %% cycle_days)
}


.msimGD_first_nonnull <- function(...) {
  vals <- list(...)
  for (val in vals) {
    if (!is.null(val)) {
      return(val)
    }
  }
  NULL
}


.msimGD_capture_time_dependent_signature <- function(parameters) {
  if (!is.list(parameters) || length(parameters) < 1L) {
    stop("parameters must be a non-empty list of node parameter lists.", call. = FALSE)
  }

  node_signatures <- lapply(seq_along(parameters), function(node_index) {
    node_parameters <- parameters[[node_index]]
    cube_info <- malariasimulationGD:::cube_genotype_info(node_parameters$cube)
    species_names <- as.character(node_parameters$species)
    species_signatures <- lapply(
      seq_along(species_names),
      function(species_i) {
        node <- malariasimulationGD:::native_prepare_carrying_capacity_node(
          node_parameters,
          species_i = species_i,
          cube_info = cube_info
        )
        list(
          k0 = as.numeric(node$k0),
          model_seasonality = isTRUE(node$model_seasonality),
          g0 = as.numeric(node$g0),
          g = as.numeric(node$g),
          h = as.numeric(node$h),
          R_bar = as.numeric(node$R_bar),
          rainfall_floor = as.numeric(node$rainfall_floor),
          carrying_capacity_timesteps = as.integer(node$carrying_capacity_timesteps),
          carrying_capacity_scalers = as.numeric(node$carrying_capacity_scalers)
        )
      }
    )
    names(species_signatures) <- species_names

    vector_control_time_offset <- node_parameters$vector_control_time_offset
    if (is.null(vector_control_time_offset)) {
      vector_control_time_offset <- 0L
    }

    list(
      node_index = as.integer(node_index),
      vector_control_time_offset = as.integer(vector_control_time_offset),
      species = species_signatures
    )
  })
  names(node_signatures) <- seq_along(node_signatures)

  list(
    signature_type = "baseline_time_dependent_signature_v1",
    seasonal_cycle_days = 365L,
    nodes = node_signatures
  )
}


.msimGD_contact_multiplier_values <- function(node_parameters) {
  species_names <- as.character(node_parameters$species)
  if (length(species_names) < 1L) {
    stop("node parameter list must contain at least one species.", call. = FALSE)
  }

  contact_multiplier <- node_parameters$contact_multiplier
  if (is.null(contact_multiplier)) {
    values <- rep(1, length(species_names))
    names(values) <- species_names
    return(values)
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
      "`contact_multiplier` must be NULL or a single positive finite number.",
      call. = FALSE
    )
  }

  values <- rep(value, length(species_names))
  names(values) <- species_names
  values
}


.msimGD_contact_covariates <- function(node_parameters) {
  covariates <- node_parameters$contact_multiplier_covariates
  if (is.null(covariates)) {
    return(NULL)
  }

  if (is.atomic(covariates) && is.numeric(covariates)) {
    out <- as.numeric(covariates)
    names(out) <- names(covariates)
    return(out)
  }

  if (is.list(covariates) && !is.data.frame(covariates) &&
      all(vapply(covariates, function(x) {
        is.numeric(x) && length(x) == 1L && is.finite(x)
      }, logical(1)))) {
    out <- unlist(covariates, use.names = TRUE)
    out <- as.numeric(out)
    names(out) <- names(covariates)
    return(out)
  }

  covariates
}


.msimGD_contact_multiplier_by_node <- function(parameters) {
  out <- lapply(parameters, .msimGD_contact_multiplier_values)
  names(out) <- seq_along(out)
  if (all(lengths(out) == 1L)) {
    flattened <- vapply(out, function(x) as.numeric(x[[1L]]), numeric(1))
    names(flattened) <- names(out)
    return(flattened)
  }
  out
}


.msimGD_contact_covariates_by_node <- function(parameters) {
  out <- lapply(parameters, .msimGD_contact_covariates)
  names(out) <- seq_along(out)
  if (all(vapply(out, is.null, logical(1)))) {
    return(NULL)
  }
  out
}


.msimGD_capture_contact_signature <- function(parameters) {
  if (!is.list(parameters) || length(parameters) < 1L) {
    stop("parameters must be a non-empty list of node parameter lists.", call. = FALSE)
  }

  node_signatures <- lapply(seq_along(parameters), function(node_index) {
      node_parameters <- parameters[[node_index]]
      node_signature <- list(
        node_index = as.integer(node_index),
      contact_multiplier = .msimGD_contact_multiplier_values(node_parameters)
      )

    covariates <- .msimGD_contact_covariates(node_parameters)
    if (!is.null(covariates)) {
      node_signature$contact_multiplier_covariates <- covariates
    }

    label <- node_parameters$contact_multiplier_label
    if (!is.null(label)) {
      node_signature$contact_multiplier_label <- as.character(label)[[1L]]
    }

    source <- node_parameters$contact_multiplier_source
    if (!is.null(source)) {
      node_signature$contact_multiplier_source <- as.character(source)[[1L]]
    }

    node_signature
  })
  names(node_signatures) <- seq_along(node_signatures)

  hook <- .msimGD_first_nonnull(
    parameters[[1L]]$contact_multiplier_hook,
    "human_blood_meal_rate"
  )
  out <- list(
    signature_type = "baseline_contact_signature_v1",
    contact_hook = as.character(hook)[[1L]],
    nodes = node_signatures
  )

  label <- parameters[[1L]]$contact_multiplier_label
  if (!is.null(label)) {
    out$label <- as.character(label)[[1L]]
  }

  source <- parameters[[1L]]$contact_multiplier_source
  if (!is.null(source)) {
    out$source <- as.character(source)[[1L]]
  }

  effect_spec <- parameters[[1L]]$contact_multiplier_effect_spec
  if (!is.null(effect_spec)) {
    out$effect_spec <- effect_spec
  }

  normalize_flag <- parameters[[1L]]$contact_multiplier_normalize
  if (!is.null(normalize_flag)) {
    out$contact_multiplier_normalize <- isTRUE(normalize_flag)
  }

  out
}


.msimGD_signature_contact_values <- function(node_signature) {
  node_signature$contact_multiplier
}


.msimGD_contact_signature_is_identity <- function(signature, tolerance = 1e-12) {
  if (is.null(signature) || !is.list(signature) || is.null(signature$nodes)) {
    return(FALSE)
  }

  isTRUE(all(vapply(
    signature$nodes,
    function(node) {
      values <- as.numeric(.msimGD_signature_contact_values(node))
      length(values) > 0L &&
        all(is.finite(values)) &&
        all(abs(values - 1) <= tolerance)
    },
    logical(1)
  )))
}


.msimGD_canonicalize_contact_signature <- function(signature) {
  if (is.null(signature) || !is.list(signature) || is.null(signature$nodes)) {
    return(signature)
  }

  signature
}


.msimGD_time_dependent_signature_is_static <- function(signature) {
  if (is.null(signature) || !is.list(signature) || is.null(signature$nodes)) {
    return(FALSE)
  }

  isTRUE(all(vapply(
    signature$nodes,
    function(node) {
      offset <- as.integer(node$vector_control_time_offset)
      all(
        !vapply(node$species, function(species_signature) {
          isTRUE(species_signature$model_seasonality) ||
            length(species_signature$carrying_capacity_timesteps) > 0L
        }, logical(1))
      ) && identical(offset, 0L)
    },
    logical(1)
  )))
}


.msimGD_validate_checkpoint_time_dependent_signature <- function(checkpoint_metadata,
                                                                 current_signature,
                                                                 artifact_label = "baseline_checkpoint") {
  if (is.null(checkpoint_metadata)) {
    return(invisible(NULL))
  }

  checkpoint_signature <- checkpoint_metadata$baseline_time_dependent_signature
  if (is.null(checkpoint_signature)) {
    if (.msimGD_time_dependent_signature_is_static(current_signature)) {
      return(invisible(NULL))
    }

    stop(
      paste(
        artifact_label,
        "does not contain `baseline_time_dependent_signature`",
        "metadata required for seasonal or otherwise time-varying carrying",
        "capacity baselines. Rebuild the checkpoint or checkpoint library",
        "with the current seasonality-aware code."
      ),
      call. = FALSE
    )
  }

  signature_check <- all.equal(
    checkpoint_signature,
    current_signature,
    tolerance = 1e-8,
    check.attributes = FALSE
  )
  if (!isTRUE(signature_check)) {
    stop(
      paste(
        artifact_label,
        "baseline_time_dependent_signature does not match",
        "the current model setup:",
        paste(signature_check, collapse = " ")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}


.msimGD_validate_baseline_contact_signature <- function(metadata,
                                                        current_signature,
                                                        artifact_label = "baseline_checkpoint") {
  if (is.null(metadata)) {
    return(invisible(NULL))
  }

  stored_signature <- metadata$baseline_contact_signature
  if (is.null(stored_signature)) {
    if (.msimGD_contact_signature_is_identity(current_signature)) {
      return(invisible(NULL))
    }

    stop(
      paste(
        artifact_label,
        "does not contain `baseline_contact_signature` metadata required",
        "for non-identity contact-multiplier baselines. Rebuild the baseline",
        "artifact with the current contact-aware code."
      ),
      call. = FALSE
    )
  }

  signature_check <- all.equal(
    .msimGD_canonicalize_contact_signature(stored_signature),
    .msimGD_canonicalize_contact_signature(current_signature),
    tolerance = 1e-8,
    check.attributes = FALSE
  )
  if (!isTRUE(signature_check)) {
    stop(
      paste(
        artifact_label,
        "baseline_contact_signature does not match the current model setup:",
        paste(signature_check, collapse = " ")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}


.msimGD_validate_stationary_human_library_signatures <- function(
    human_initialization_library,
    current_time_signature,
    current_contact_signature
) {
  if (is.null(human_initialization_library)) {
    return(invisible(NULL))
  }

  library_object <- malariasimulationGD:::stationary_human_initializer_validate_library(
    malariasimulationGD:::stationary_human_initializer_load_library(human_initialization_library)
  )
  metadata <- library_object$metadata
  if (is.null(metadata)) {
    metadata <- list()
  }

  .msimGD_validate_checkpoint_time_dependent_signature(
    checkpoint_metadata = metadata,
    current_signature = current_time_signature,
    artifact_label = "human_initialization_library"
  )
  .msimGD_validate_baseline_contact_signature(
    metadata = metadata,
    current_signature = current_contact_signature,
    artifact_label = "human_initialization_library"
  )

  invisible(NULL)
}


.msimGD_prepare_truth_run <- function(setup, cube, NF, NH, tmax, mu, p_move, release,
                                      theta, init_EIR,
                                      prevalence_band,
                                      infection_band,
                                      clinical_band,
                                      warmup_days,
                                      parameter_modifier,
                                      baseline_checkpoint_state = NULL,
                                      baseline_checkpoint_metadata = NULL,
                                      human_initialization_library = NULL) {
  D <- as.matrix(setup$D)
  n_nodes <- nrow(D)
  warmup_days <- as.integer(warmup_days)
  if (is.na(warmup_days) || warmup_days < 0L) {
    stop("warmup_days must be a single integer >= 0.", call. = FALSE)
  }

  checkpoint_timesteps <- NULL
  if (!is.null(baseline_checkpoint_state)) {
    checkpoint_timesteps <- as.integer(baseline_checkpoint_state$timesteps)
    if (is.na(checkpoint_timesteps) || checkpoint_timesteps < 0L) {
      stop("baseline_checkpoint$timesteps must be a single integer >= 0.", call. = FALSE)
    }
    if (warmup_days > 0L && checkpoint_timesteps != warmup_days) {
      stop(
        sprintf(
          paste(
            "baseline_checkpoint timesteps (%d) do not match requested warmup_days (%d).",
            "Use the matching checkpoint or set warmup_days accordingly."
          ),
          checkpoint_timesteps,
          warmup_days
        ),
        call. = FALSE
      )
    }
  }

  effective_warmup_days <- if (!is.null(checkpoint_timesteps)) checkpoint_timesteps else warmup_days
  total_timesteps <- as.integer(tmax) + effective_warmup_days
  if (total_timesteps <= 0L) {
    stop("tmax + effective warmup must be positive.", call. = FALSE)
  }

  if (length(NH) == 1L) NH <- rep(NH, n_nodes)
  stopifnot(length(NH) == n_nodes)
  if (!is.null(NF) && length(NF) == 1L) NF <- rep(NF, n_nodes)
  if (!is.null(NF)) stopifnot(length(NF) == n_nodes)

  muF <- theta$muF
  # Inlined from customMGDrive2::calc_move_rate to avoid the customMGDrive2
  # dependency (formula: calc_move_rate(mu, P) = (P * mu) / (1 - P)).
  move_rate <- (p_move * muF) / (1 - p_move)
  mov <- mosquito_movement_from_mu(
    D = D, mu = mu,
    move_rates = rep(move_rate, n_nodes),
    attractiveness = rep(1, n_nodes),
    verbose = FALSE
  )

  parameters <- vector("list", n_nodes)
  equilibrium_total_M <- numeric(n_nodes)
  effective_total_M <- numeric(n_nodes)
  init_foim <- numeric(n_nodes)

  for (nd in seq_len(n_nodes)) {
    bp <- msimGD_make_base_parameters(
      NH = NH[nd],
      theta = theta,
      individual_mosquitoes = TRUE
    )
    bp$move_probs <- mov$mosquito_move_probs
    bp$move_rates <- mov$mosquito_move_rates
    if (length(prevalence_band$min) > 0L) {
      bp$prevalence_rendering_min_ages <- prevalence_band$min
      bp$prevalence_rendering_max_ages <- prevalence_band$max
    }
    if (length(infection_band$min) > 0L) {
      bp$incidence_rendering_min_ages <- infection_band$min
      bp$incidence_rendering_max_ages <- infection_band$max
    }
    if (length(clinical_band$min) > 0L) {
      bp$clinical_incidence_rendering_min_ages <- clinical_band$min
      bp$clinical_incidence_rendering_max_ages <- clinical_band$max
    }
    if (!is.null(parameter_modifier)) {
      bp <- parameter_modifier(bp, node_index = nd, warmup_days = effective_warmup_days)
      if (is.null(bp) || !is.list(bp)) {
        stop("parameter_modifier must return a parameter list.", call. = FALSE)
      }
    }
    bp <- set_equilibrium(bp, init_EIR = init_EIR)
    equilibrium_total_M[nd] <- unname(as.numeric(bp$total_M))
    init_foim[nd] <- unname(as.numeric(bp$init_foim))

    if (!is.null(NF)) {
      bp <- parameterise_total_M(bp, NF[nd])
    }
    if (!is.null(human_initialization_library)) {
      bp <- malariasimulationGD:::set_stationary_human_initialization_library(
        bp,
        library = human_initialization_library,
        node_index = nd
      )
    }
    bp$cube <- cube
    effective_total_M[nd] <- unname(as.numeric(bp$total_M))
    parameters[[nd]] <- bp
  }

  if (.msimGD_release_is_active(release)) {
    release_sex <- if (!is.null(release$sex)) release$sex else release$stage
    for (i in release$nodes) {
      parameters[[i]] <- set_releases(parameters[[i]], list(
        releasesStart = effective_warmup_days + as.integer(release$time),
        releasesNumber = 1L,
        releaseCount = as.integer(release$size),
        releaseSex = release_sex,
        releaseGenotype = release$genotype,
        releasesInterval = 0L
      ))
    }
  }

  names(equilibrium_total_M) <- seq_len(n_nodes)
  names(effective_total_M) <- seq_len(n_nodes)
  names(init_foim) <- seq_len(n_nodes)
  baseline_time_dependent_signature <- .msimGD_capture_time_dependent_signature(parameters)
  baseline_contact_signature <- .msimGD_capture_contact_signature(parameters)
  contact_multiplier_by_node <- .msimGD_contact_multiplier_by_node(parameters)
  contact_covariates_by_node <- .msimGD_contact_covariates_by_node(parameters)
  .msimGD_validate_checkpoint_time_dependent_signature(
    checkpoint_metadata = baseline_checkpoint_metadata,
    current_signature = baseline_time_dependent_signature
  )
  .msimGD_validate_baseline_contact_signature(
    metadata = baseline_checkpoint_metadata,
    current_signature = baseline_contact_signature
  )
  .msimGD_validate_stationary_human_library_signatures(
    human_initialization_library = human_initialization_library,
    current_time_signature = baseline_time_dependent_signature,
    current_contact_signature = baseline_contact_signature
  )

  list(
    n_nodes = n_nodes,
    parameters = parameters,
    total_timesteps = total_timesteps,
    effective_warmup_days = effective_warmup_days,
    equilibrium_total_M = equilibrium_total_M,
    effective_total_M = effective_total_M,
    init_foim = init_foim,
    baseline_time_dependent_signature = baseline_time_dependent_signature,
    baseline_contact_signature = baseline_contact_signature,
    contact_multiplier_by_node = contact_multiplier_by_node,
    contact_covariates_by_node = contact_covariates_by_node
  )
}


#' Build and save a resumable stationary metapop baseline state for later truth runs.
#'
#' @param setup list with distance matrix `D`
#' @param cube inheritance cube
#' @param NF adult females per node
#' @param NH humans per node
#' @param mu mean jump distance
#' @param p_move P(move before death)
#' @param theta lifecycle parameters
#' @param init_EIR endemic EIR passed into `set_equilibrium()`
#' @param warmup_days number of baseline days to bake into the checkpoint
#' @param parameter_modifier optional baseline modifier
#' @param seed RNG seed
#' @param checkpoint_path optional path to write the checkpoint RDS
#' @return checkpoint object with `$state` and `$metadata`
msimGD_build_baseline_checkpoint <- function(setup, cube, NF, NH, mu, p_move,
                                             theta, init_EIR = 5,
                                             warmup_days = 0L,
                                             parameter_modifier = NULL,
                                             seed = NULL,
                                             checkpoint_path = NULL) {
  if (!is.null(seed)) set.seed(seed)

  warmup_days <- as.integer(warmup_days)
  if (is.na(warmup_days) || warmup_days <= 0L) {
    stop("warmup_days must be a single integer > 0 to build a checkpoint.", call. = FALSE)
  }

  prep <- .msimGD_prepare_truth_run(
    setup = setup,
    cube = cube,
    NF = NF,
    NH = NH,
    tmax = 0L,
    mu = mu,
    p_move = p_move,
    release = NULL,
    theta = theta,
    init_EIR = init_EIR,
    prevalence_band = list(min = integer(0), max = integer(0)),
    infection_band = list(min = integer(0), max = integer(0)),
    clinical_band = list(min = integer(0), max = integer(0)),
    warmup_days = warmup_days,
    parameter_modifier = parameter_modifier,
    baseline_checkpoint_state = NULL,
    baseline_checkpoint_metadata = NULL
  )

  sim <- run_metapop_simulation(
    timesteps = prep$effective_warmup_days,
    parameters = prep$parameters,
    mixing_tt = 1,
    export_mixing = list(diag(prep$n_nodes)),
    import_mixing = list(diag(prep$n_nodes)),
    p_captured_tt = 1,
    p_captured = list(matrix(0, prep$n_nodes, prep$n_nodes)),
    p_success = 0,
    return_state = TRUE
  )

  checkpoint <- list(
    state = sim$state,
    metadata = list(
      seed = seed,
      init_EIR = init_EIR,
      timesteps = prep$effective_warmup_days,
      warmup_days = prep$effective_warmup_days,
      seasonal_cycle_days = 365L,
      seasonal_phase_day = .msimGD_seasonal_phase_day(prep$effective_warmup_days),
      baseline_time_dependent_signature = prep$baseline_time_dependent_signature,
      baseline_contact_signature = prep$baseline_contact_signature,
      contact_multiplier_by_node = prep$contact_multiplier_by_node,
      contact_covariates_by_node = prep$contact_covariates_by_node,
      equilibrium_total_M = prep$equilibrium_total_M,
      effective_total_M = prep$effective_total_M,
      init_foim = prep$init_foim
    )
  )

  if (!is.null(checkpoint_path)) {
    saveRDS(checkpoint, checkpoint_path)
  }

  checkpoint
}


#' Run one stochastic msimGD metapop simulation
#'
#' @param setup  list with D (distance matrix) from load_busia_landscape()
#' @param cube   inheritance cube (TP13)
#' @param NF     adult females per node. If `NULL`, use the abundance implied by
#'               `init_EIR` separately for each node. If supplied, may be either
#'               a scalar or a length-`n_nodes` vector and overrides the local
#'               equilibrium abundance after human/FOIM equilibrium is computed.
#' @param NH     human population per node (scalar or length-`n_nodes` vector)
#' @param tmax   simulation duration (days)
#' @param mu     mean jump distance (km)
#' @param p_move P(move before death)
#' @param release list(nodes, time, size, sex/stage, genotype). If `NULL` or
#'   empty, no release is applied.
#' @param theta  lifecycle parameter list from customMGDrive2 — MUST be provided
#'               so that the msimGD stochastic ento dynamics match the
#'               customMGDrive2 deterministic forward model (mean = stochastic mean)
#' @param init_EIR  initial EIR for human equilibrium (default 5)
#' @param prevalence_rendering_min_age minimum age in days for prevalence
#'   outputs. Use `NULL` to disable prevalence rendering.
#' @param prevalence_rendering_max_age maximum age in days for prevalence outputs.
#' @param infection_incidence_min_age minimum age in days for infection
#'   incidence outputs. Use `NULL` to disable infection incidence rendering.
#' @param infection_incidence_max_age maximum age in days for infection
#'   incidence outputs.
#' @param clinical_incidence_min_age minimum age in days for clinical incidence
#'   outputs. Use `NULL` to disable clinical incidence rendering.
#' @param clinical_incidence_max_age maximum age in days for clinical incidence outputs
#' @param warmup_days number of pre-observation simulation days to run and then
#'   discard before returning outputs. Release times are interpreted relative to
#'   the post-warmup timeline.
#' @param parameter_modifier optional function `(parameters, node_index,
#'   warmup_days)` applied to each node's parameter list before equilibrium
#'   initialisation. Use this to add routine treatment, bednets, or custom
#'   output renderers in one place.
#' @param human_initialization_library optional prebuilt stationary human library
#'   object (or `.rds` path) to use for true day-0 initialization. May be
#'   node-conditioned.
#' @param baseline_checkpoint optional resumable state object or RDS path created
#'   by `msimGD_build_baseline_checkpoint()`. When supplied, the truth run
#'   starts from that saved stationary state instead of replaying warmup.
#' @param seed   RNG seed
#' @return list of per-node data.frames from run_metapop_simulation()
msimGD_run_truth <- function(setup, cube, NF, NH, tmax, mu, p_move, release,
                             theta, init_EIR = 5,
                             prevalence_rendering_min_age = NULL,
                             prevalence_rendering_max_age = NULL,
                             infection_incidence_min_age = NULL,
                             infection_incidence_max_age = NULL,
                             clinical_incidence_min_age = 182L,
                             clinical_incidence_max_age = 5475L,
                             warmup_days = 0L,
                             parameter_modifier = NULL,
                             human_initialization_library = NULL,
                             baseline_checkpoint = NULL,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  prevalence_band <- .msimGD_normalise_optional_age_band(
    prevalence_rendering_min_age,
    prevalence_rendering_max_age,
    "Prevalence rendering"
  )
  infection_band <- .msimGD_normalise_optional_age_band(
    infection_incidence_min_age,
    infection_incidence_max_age,
    "Infection incidence"
  )
  clinical_band <- .msimGD_normalise_optional_age_band(
    clinical_incidence_min_age,
    clinical_incidence_max_age,
    "Clinical incidence"
  )

  baseline_checkpoint_info <- .msimGD_normalise_checkpoint(baseline_checkpoint)
  baseline_checkpoint_state <- if (!is.null(baseline_checkpoint_info)) baseline_checkpoint_info$state else NULL
  restore_random_state <- FALSE
  if (!is.null(baseline_checkpoint_info)) {
    checkpoint_seed <- baseline_checkpoint_info$metadata$seed
    if (!is.null(checkpoint_seed) && !is.null(seed)) {
      restore_random_state <- identical(as.integer(checkpoint_seed), as.integer(seed))
      if (!restore_random_state) {
        warning(
          sprintf(
            paste(
              "Using saved checkpoint state built with checkpoint seed (%d)",
              "and a different truth seed (%d);",
              "post-checkpoint stochastic continuation will use the new truth seed",
              "rather than replaying the original checkpoint RNG stream."
            ),
            as.integer(checkpoint_seed),
            as.integer(seed)
          ),
          call. = FALSE
        )
      }
    }
  }
  prep <- .msimGD_prepare_truth_run(
    setup = setup,
    cube = cube,
    NF = NF,
    NH = NH,
    tmax = tmax,
    mu = mu,
    p_move = p_move,
    release = release,
    theta = theta,
    init_EIR = init_EIR,
    prevalence_band = prevalence_band,
    infection_band = infection_band,
    clinical_band = clinical_band,
    warmup_days = warmup_days,
    parameter_modifier = parameter_modifier,
    baseline_checkpoint_state = baseline_checkpoint_state,
    baseline_checkpoint_metadata = if (!is.null(baseline_checkpoint_info)) baseline_checkpoint_info$metadata else NULL,
    human_initialization_library = human_initialization_library
  )

  # --- Run ---
  result <- run_metapop_simulation(
    timesteps = prep$total_timesteps,
    parameters = prep$parameters,
    mixing_tt = 1,
    export_mixing = list(diag(prep$n_nodes)),
    import_mixing = list(diag(prep$n_nodes)),
    p_captured_tt = 1,
    p_captured = list(matrix(0, prep$n_nodes, prep$n_nodes)),
    p_success = 0,
    initial_state = baseline_checkpoint_state,
    restore_random_state = restore_random_state
  )

  attr(result, "equilibrium_total_M") <- prep$equilibrium_total_M
  attr(result, "effective_total_M") <- prep$effective_total_M
  attr(result, "init_foim") <- prep$init_foim
  attr(result, "baseline_time_dependent_signature") <- prep$baseline_time_dependent_signature
  attr(result, "baseline_contact_signature") <- prep$baseline_contact_signature
  attr(result, "contact_multiplier_by_node") <- prep$contact_multiplier_by_node
  attr(result, "contact_covariates_by_node") <- prep$contact_covariates_by_node
  attr(result, "prevalence_rendering_min_age") <- prevalence_band$min
  attr(result, "prevalence_rendering_max_age") <- prevalence_band$max
  attr(result, "infection_incidence_min_age") <- infection_band$min
  attr(result, "infection_incidence_max_age") <- infection_band$max
  attr(result, "clinical_incidence_min_age") <- clinical_band$min
  attr(result, "clinical_incidence_max_age") <- clinical_band$max
  attr(result, "warmup_days") <- prep$effective_warmup_days

  if (is.null(baseline_checkpoint_state)) {
    return(.msimGD_crop_result(result, warmup_days = prep$effective_warmup_days))
  }

  result
}


#' Extract infection or clinical incidence from msimGD metapop output and
#' aggregate into observation windows, along with a matching age-band population
#' offset.
#'
#' @param result     output of msimGD_run_truth() (list of per-node data.frames)
#' @param n_nodes    number of nodes (58)
#' @param obs_every  window width in days (default 14)
#' @param NH         fallback human population per node if age-band counts are absent
#' @param obs_start  first window ends at this day (default obs_every)
#' @param endpoint which epi endpoint to extract: `"infection"` or `"clinical"`
#' @param age_min optional minimum age in days used to choose a specific
#'   endpoint rendering column.
#' @param age_max optional maximum age in days used to choose a specific
#'   endpoint rendering column.
#' @return data.frame(node, time, w, Y, offset_pop, offset_person_days)
#'   where `offset_pop` is the mean daily population in the rendered endpoint
#'   age band over the observation window, so that `w * offset_pop`
#'   equals the age-band person-days used by the Poisson offset.
msimGD_extract_epi_obs <- function(result, n_nodes, obs_every = 14L,
                                   NH = 770L, obs_start = NULL,
                                   endpoint = c("clinical", "infection"),
                                   age_min = NULL, age_max = NULL) {
  if (is.null(obs_start)) obs_start <- obs_every
  endpoint <- match.arg(endpoint)

  first_df <- result[[1]]
  epi_info <- .msimGD_find_epi_column(
    first_df = first_df,
    endpoint = endpoint,
    age_min = age_min,
    age_max = age_max
  )
  has_age_col <- epi_info$age_col %in% names(first_df)

  tmax <- nrow(first_df)

  # Window end times
  window_ends <- seq(from = obs_start, to = tmax, by = obs_every)

  # Pre-allocate
  n_windows <- length(window_ends)
  out <- data.frame(
    node = rep(seq_len(n_nodes), each = n_windows),
    time = rep(window_ends, times = n_nodes),
    w = NA_integer_,
    Y = NA_integer_,
    offset_pop = NA_real_,
    offset_person_days = NA_real_
  )

  for (nd in seq_len(n_nodes)) {
    daily_inc <- result[[nd]][[epi_info$epi_col]]
    # NA timesteps = renderer didn't fire -> no endpoint events; treat as 0
    daily_inc[is.na(daily_inc)] <- 0L
    if (has_age_col) {
      daily_pop <- result[[nd]][[epi_info$age_col]]
      if (anyNA(daily_pop)) {
        stop("Age-band population column ", epi_info$age_col,
             " contains NA values for node ", nd, ".", call. = FALSE)
      }
    } else {
      node_NH <- if (length(NH) > 1L) NH[nd] else NH
      daily_pop <- rep(as.numeric(node_NH), tmax)
    }
    for (wi in seq_along(window_ends)) {
      t_end <- window_ends[wi]
      t_start <- t_end - obs_every + 1L
      t_start <- max(t_start, 1L)
      t_end <- min(t_end, length(daily_inc))
      window_len <- t_end - t_start + 1L
      pop_person_days <- sum(daily_pop[t_start:t_end])
      row_idx <- (nd - 1L) * n_windows + wi
      out$w[row_idx] <- window_len
      out$Y[row_idx] <- sum(daily_inc[t_start:t_end])
      out$offset_person_days[row_idx] <- pop_person_days
      out$offset_pop[row_idx] <- pop_person_days / window_len
    }
  }

  attr(out, "epi_endpoint") <- endpoint
  attr(out, "epi_column") <- epi_info$epi_col
  attr(out, "age_population_col") <- if (has_age_col) epi_info$age_col else NA_character_
  attr(out, "age_suffix") <- epi_info$age_suffix

  out
}


#' Extract daily prevalence observations from msimGD metapop output.
#'
#' @param result output of msimGD_run_truth() (list of per-node data.frames)
#' @param n_nodes number of nodes
#' @param age_min minimum age in days
#' @param age_max maximum age in days
#' @param method detection method: `"lm"` or `"pcr"`
#' @param times optional integer vector of monitoring times to retain. By
#'   default, all daily rendered times are returned.
#' @return data.frame(node, time, n_detect, n_age, prevalence)
msimGD_extract_prevalence_obs <- function(result, n_nodes, age_min, age_max,
                                          method = c("lm", "pcr"),
                                          times = NULL) {
  method <- match.arg(method)

  first_df <- result[[1]]
  prev_info <- .msimGD_find_prevalence_columns(
    first_df = first_df,
    method = method,
    age_min = age_min,
    age_max = age_max
  )

  tmax <- nrow(first_df)
  out <- data.frame(
    node = rep(seq_len(n_nodes), each = tmax),
    time = rep(seq_len(tmax), times = n_nodes),
    n_detect = NA_real_,
    n_age = NA_real_,
    prevalence = NA_real_
  )

  for (nd in seq_len(n_nodes)) {
    df <- result[[nd]]
    n_detect <- df[[prev_info$detect_col]]
    n_age <- df[[prev_info$age_col]]
    n_detect[is.na(n_detect)] <- 0
    if (anyNA(n_age)) {
      stop(sprintf(
        "Age-band population column `%s` contains NA values for node %d.",
        prev_info$age_col, nd
      ), call. = FALSE)
    }

    idx <- ((nd - 1L) * tmax + 1L):(nd * tmax)
    out$n_detect[idx] <- n_detect
    out$n_age[idx] <- n_age
    out$prevalence[idx] <- ifelse(n_age > 0, n_detect / n_age, NA_real_)
  }

  attr(out, "prevalence_method") <- method
  attr(out, "prevalence_detect_col") <- prev_info$detect_col
  attr(out, "age_population_col") <- prev_info$age_col
  attr(out, "age_suffix") <- prev_info$age_suffix

  if (!is.null(times)) {
    times <- as.integer(times)
    if (length(times) < 1L || anyNA(times) || any(times < 1L)) {
      stop("times must be a non-empty integer vector with values >= 1.", call. = FALSE)
    }
    missing_times <- setdiff(times, seq_len(tmax))
    if (length(missing_times) > 0L) {
      stop(
        sprintf(
          "Requested prevalence times are outside the rendered range 1..%d. Examples: %s",
          tmax,
          paste(utils::head(missing_times, 10L), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    out <- out[out$time %in% unique(times), , drop = FALSE]
    rownames(out) <- NULL
  }

  out
}


msimGD_mean_prevalence <- function(prevalence_obs) {
  if (!all(c("n_detect", "n_age") %in% names(prevalence_obs))) {
    stop("prevalence_obs must contain `n_detect` and `n_age` columns.",
         call. = FALSE)
  }
  total_n_age <- sum(prevalence_obs$n_age)
  if (!is.finite(total_n_age) || total_n_age <= 0) {
    stop("Total age-band population must be positive to compute prevalence.",
         call. = FALSE)
  }
  sum(prevalence_obs$n_detect) / total_n_age
}


msimGD_mean_eir_ppy <- function(result, NH) {
  if (!is.list(result) || length(result) < 1L) {
    stop("result must be a non-empty list of node-level outputs.", call. = FALSE)
  }

  n_nodes <- length(result)
  NH <- as.numeric(NH)
  if (length(NH) == 1L) {
    NH <- rep.int(NH, n_nodes)
  }
  if (length(NH) != n_nodes || any(!is.finite(NH)) || any(NH <= 0)) {
    stop("NH must be a positive scalar or a positive vector matching result.",
         call. = FALSE)
  }

  tmax <- nrow(result[[1]])
  if (!is.numeric(tmax) || length(tmax) != 1L || tmax < 1L) {
    stop("Each node-level output must have at least one time step.", call. = FALSE)
  }

  total_eir_per_day <- numeric(tmax)
  total_nh <- sum(NH)

  for (nd in seq_len(n_nodes)) {
    df <- result[[nd]]
    if (!is.data.frame(df) || nrow(df) != tmax) {
      stop("All node-level outputs must be data frames with the same number of rows.",
           call. = FALSE)
    }
    eir_cols <- grep("^EIR_", names(df), value = TRUE)
    if (!length(eir_cols)) {
      next
    }
    total_eir_per_day <- total_eir_per_day + rowSums(df[, eir_cols, drop = FALSE])
  }

  mean(total_eir_per_day / total_nh * 365)
}


#' Extract carrier frequency time series from msimGD metapop output.
#'
#' @param result         output of msimGD_run_truth()
#' @param n_nodes        number of nodes
#' @param carrier_allele allele character to grep for (default "H")
#' @return data.frame(time, node, rep, carrier_freq) matching mg_extract_carrier_frequency() format
msimGD_extract_carrier_freq <- function(result, n_nodes, carrier_allele = "H") {
  dfs <- vector("list", n_nodes)

  for (nd in seq_len(n_nodes)) {
    fc <- attr(result[[nd]], "mosquito_genotype_counts_female")
    if (is.null(fc)) {
      stop(sprintf("Node %d has no mosquito_genotype_counts_female attribute. ", nd),
           "Ensure the simulation records genotype history.", call. = FALSE)
    }

    gnames <- colnames(fc)
    H_cols <- grep(carrier_allele, gnames)
    n_times <- nrow(fc)

    row_totals <- rowSums(fc)
    carrier_totals <- if (length(H_cols) > 0) rowSums(fc[, H_cols, drop = FALSE]) else rep(0, n_times)
    cf <- ifelse(row_totals > 0, carrier_totals / row_totals, 0)

    dfs[[nd]] <- data.frame(
      time = seq_len(n_times),
      node = nd,
      rep = 1L,
      carrier_freq = cf
    )
  }

  do.call(rbind, dfs)
}
