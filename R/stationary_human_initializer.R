.stationary_human_initializer_cache <- new.env(parent = emptyenv())

human_initialization_mode <- function(parameters) {
  mode <- parameters$human_initialization
  if (is.null(mode)) {
    return("equilibrium")
  }

  mode <- as.character(mode)
  if (length(mode) != 1L || is.na(mode)) {
    stop("human_initialization must be a single non-missing string.", call. = FALSE)
  }

  mode <- match.arg(mode, c("equilibrium", "stochastic_resample"))
  mode
}

set_stochastic_human_initialization <- function(parameters,
                                                burnin_timesteps,
                                                n_snapshots = 1L,
                                                snapshot_spacing = 0L,
                                                library_population = NULL,
                                                seed = 1L) {
  burnin_timesteps <- as.integer(burnin_timesteps)
  if (length(burnin_timesteps) != 1L || is.na(burnin_timesteps) || burnin_timesteps <= 0L) {
    stop("burnin_timesteps must be a single integer > 0.", call. = FALSE)
  }

  n_snapshots <- as.integer(n_snapshots)
  if (length(n_snapshots) != 1L || is.na(n_snapshots) || n_snapshots <= 0L) {
    stop("n_snapshots must be a single integer > 0.", call. = FALSE)
  }

  snapshot_spacing <- as.integer(snapshot_spacing)
  if (length(snapshot_spacing) != 1L || is.na(snapshot_spacing) || snapshot_spacing < 0L) {
    stop("snapshot_spacing must be a single integer >= 0.", call. = FALSE)
  }
  if (n_snapshots > 1L && snapshot_spacing <= 0L) {
    stop("snapshot_spacing must be > 0 when n_snapshots > 1.", call. = FALSE)
  }

  if (!is.null(library_population)) {
    library_population <- as.integer(library_population)
    if (length(library_population) != 1L || is.na(library_population) || library_population <= 0L) {
      stop("library_population must be NULL or a single integer > 0.", call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (length(seed) != 1L || is.na(seed)) {
      stop("seed must be NULL or a single integer.", call. = FALSE)
    }
  }

  parameters$human_initialization <- "stochastic_resample"
  parameters$human_initialization_burnin_timesteps <- burnin_timesteps
  parameters$human_initialization_n_snapshots <- n_snapshots
  parameters$human_initialization_snapshot_spacing <- snapshot_spacing
  parameters$human_initialization_library_population <- library_population
  parameters$human_initialization_seed <- seed
  parameters
}

clear_stationary_human_initializer_cache <- function() {
  rm(list = ls(envir = .stationary_human_initializer_cache), envir = .stationary_human_initializer_cache)
  invisible(NULL)
}

set_stationary_human_initialization_library <- function(parameters,
                                                        library,
                                                        node_index = NULL) {
  if (is.null(library)) {
    stop("library must be a stationary human library object or .rds path.", call. = FALSE)
  }

  if (!is.null(node_index)) {
    node_index <- as.integer(node_index)
    if (length(node_index) != 1L || is.na(node_index) || node_index <= 0L) {
      stop("node_index must be NULL or a single integer > 0.", call. = FALSE)
    }
  }

  parameters$human_initialization <- "stochastic_resample"
  parameters$human_initialization_library <- library
  parameters$human_initialization_node_index <- node_index
  parameters
}

stationary_human_initializer_load_library <- function(library) {
  if (is.null(library)) {
    return(NULL)
  }

  if (is.character(library)) {
    if (length(library) != 1L || is.na(library) || library == "") {
      stop("human_initialization_library path must be a single non-empty string.", call. = FALSE)
    }
    if (!file.exists(library)) {
      stop(sprintf("human_initialization_library file not found: %s", library), call. = FALSE)
    }
    library <- readRDS(library)
  }

  library
}

stationary_human_initializer_is_node_conditioned_library <- function(library) {
  is.list(library) && !is.null(library$node_libraries)
}

stationary_human_initializer_is_metapop_library <- function(library) {
  is.list(library) && !is.null(library$snapshots)
}

stationary_human_initializer_validate_library <- function(library) {
  if (!is.list(library)) {
    stop("human_initialization_library must be a stationary human library object or .rds path.",
         call. = FALSE)
  }

  if (!is.null(library$variables)) {
    return(library)
  }

  if (stationary_human_initializer_is_node_conditioned_library(library)) {
    return(library)
  }

  if (stationary_human_initializer_is_metapop_library(library)) {
    return(library)
  }

  stop(
    paste(
      "human_initialization_library must either contain `variables`",
      "for a shared donor pool, `node_libraries` for a node-conditioned pool,",
      "or `snapshots` for a full-metapop stationary library."
    ),
    call. = FALSE
  )
}

stationary_human_initializer_metadata_value <- function(metadata, field_names, index = NULL) {
  if (!is.list(metadata) || length(field_names) < 1L) {
    return(NULL)
  }

  for (field_name in field_names) {
    value <- metadata[[field_name]]
    if (is.null(value)) {
      next
    }

    if (!is.null(index)) {
      index <- as.integer(index)
      if (length(index) == 1L && !is.na(index) && index > 0L) {
        index_name <- as.character(index)
        if (!is.null(names(value)) && index_name %in% names(value)) {
          value <- value[[index_name]]
        } else if (length(value) >= index) {
          value <- value[[index]]
        }
      }
    }

    if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
      return(as.numeric(value))
    }
  }

  NULL
}

stationary_human_initializer_enrich_library_metadata <- function(library, node_index = NULL) {
  library <- stationary_human_initializer_validate_library(library)
  metadata <- library$metadata
  if (is.null(metadata)) {
    metadata <- list()
  }

  stationary_total_M <- stationary_human_initializer_metadata_value(
    metadata,
    field_names = c(
      "stationary_total_M",
      "stationary_total_M_by_node",
      "baseline_total_M",
      "baseline_total_M_by_node",
      "parameter_total_M",
      "effective_total_M"
    ),
    index = node_index
  )
  stationary_init_foim <- stationary_human_initializer_metadata_value(
    metadata,
    field_names = c(
      "stationary_init_foim",
      "stationary_init_foim_by_node",
      "baseline_init_foim",
      "baseline_init_foim_by_node",
      "init_foim"
    ),
    index = node_index
  )

  library$metadata <- c(
    metadata,
    list(
      stationary_total_M = stationary_total_M,
      stationary_init_foim = stationary_init_foim
    )
  )
  library
}

stationary_human_initializer_select_library <- function(library, parameters) {
  library <- stationary_human_initializer_validate_library(
    stationary_human_initializer_load_library(library)
  )

  if (stationary_human_initializer_is_metapop_library(library)) {
    stop(
      paste(
        "A full-metapop stationary human library must be resolved before",
        "node-level initialization."
      ),
      call. = FALSE
    )
  }

  if (!stationary_human_initializer_is_node_conditioned_library(library)) {
    return(stationary_human_initializer_enrich_library_metadata(library))
  }

  node_index <- parameters$human_initialization_node_index
  if (is.null(node_index)) {
    stop(
      paste(
        "A node-conditioned human_initialization_library was provided, but",
        "human_initialization_node_index is NULL."
      ),
      call. = FALSE
    )
  }

  node_index <- as.integer(node_index)
  if (length(node_index) != 1L || is.na(node_index) || node_index <= 0L) {
    stop("human_initialization_node_index must be a single integer > 0.", call. = FALSE)
  }

  node_libraries <- library$node_libraries
  node_name <- as.character(node_index)
  selected <- if (!is.null(names(node_libraries)) && node_name %in% names(node_libraries)) {
    node_libraries[[node_name]]
  } else if (length(node_libraries) >= node_index) {
    node_libraries[[node_index]]
  } else {
    NULL
  }

  if (is.null(selected)) {
    stop(
      sprintf(
        "Node-conditioned human_initialization_library has no entry for node %d.",
        node_index
      ),
      call. = FALSE
    )
  }

  selected$metadata <- c(selected$metadata, library$metadata)
  stationary_human_initializer_enrich_library_metadata(selected, node_index = node_index)
}

new_stationary_human_initializer_node_conditioned_library <- function(node_libraries,
                                                                      metadata = list()) {
  if (!is.list(node_libraries) || length(node_libraries) < 1L) {
    stop("node_libraries must be a non-empty list.", call. = FALSE)
  }

  if (is.null(names(node_libraries))) {
    names(node_libraries) <- to_char_vector(seq_along(node_libraries))
  }

  structure(
    list(
      library_type = "node_conditioned",
      node_libraries = node_libraries,
      metadata = metadata
    ),
    class = c("msimGD_stationary_human_library", "list")
  )
}

stationary_human_initializer_validate_metapop_snapshot <- function(snapshot) {
  if (!is.list(snapshot) || is.null(snapshot$node_libraries)) {
    stop("metapop snapshot must contain `node_libraries`.", call. = FALSE)
  }

  node_libraries <- snapshot$node_libraries
  if (!is.list(node_libraries) || length(node_libraries) < 1L) {
    stop("metapop snapshot `node_libraries` must be a non-empty list.", call. = FALSE)
  }

  invisible(lapply(node_libraries, stationary_human_initializer_validate_library))
  snapshot
}

new_stationary_human_initializer_metapop_library <- function(snapshots, metadata = list()) {
  if (!is.list(snapshots) || length(snapshots) < 1L) {
    stop("snapshots must be a non-empty list.", call. = FALSE)
  }

  invisible(lapply(snapshots, stationary_human_initializer_validate_metapop_snapshot))

  structure(
    list(
      library_type = "metapop_snapshot",
      snapshots = snapshots,
      metadata = metadata
    ),
    class = c("msimGD_stationary_human_library", "list")
  )
}

stationary_human_initializer_set_accepted_metapop_snapshots <- function(
    library,
    accepted_snapshot_indices
) {
  library <- stationary_human_initializer_validate_library(
    stationary_human_initializer_load_library(library)
  )
  if (!stationary_human_initializer_is_metapop_library(library)) {
    stop(
      "library must be a full-metapop stationary human library.",
      call. = FALSE
    )
  }

  accepted_snapshot_indices <- unique(as.integer(accepted_snapshot_indices))
  accepted_snapshot_indices <- accepted_snapshot_indices[
    !is.na(accepted_snapshot_indices) &
      accepted_snapshot_indices > 0L &
      accepted_snapshot_indices <= length(library$snapshots)
  ]

  library$metadata <- c(
    library$metadata,
    list(accepted_snapshot_indices = accepted_snapshot_indices)
  )
  library
}

new_stationary_checkpoint_library <- function(checkpoints, metadata = list()) {
  if (!is.list(checkpoints) || length(checkpoints) < 1L) {
    stop("checkpoints must be a non-empty list.", call. = FALSE)
  }

  invisible(lapply(checkpoints, stationary_checkpoint_library_validate))

  structure(
    list(
      library_type = "checkpoint",
      checkpoints = checkpoints,
      metadata = metadata
    ),
    class = c("msimGD_checkpoint_library", "list")
  )
}

stationary_checkpoint_library_load <- function(library) {
  if (is.null(library)) {
    return(NULL)
  }

  if (is.character(library)) {
    if (length(library) != 1L || is.na(library) || library == "") {
      stop("checkpoint library path must be a single non-empty string.", call. = FALSE)
    }
    if (!file.exists(library)) {
      stop(sprintf("checkpoint library file not found: %s", library), call. = FALSE)
    }
    library <- readRDS(library)
  }

  library
}

stationary_checkpoint_library_validate <- function(checkpoint) {
  if (is.null(checkpoint)) {
    stop("checkpoint must be a resumable state object or wrapper list.", call. = FALSE)
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
        "checkpoint must contain a resumable state with `timesteps`,",
        "`individual`, and `malariasimulationGD`."
      ),
      call. = FALSE
    )
  }

  checkpoint
}

stationary_checkpoint_library_select <- function(library,
                                                 snapshot_index = NULL,
                                                 timesteps = NULL) {
  library <- stationary_checkpoint_library_load(library)
  if (!is.list(library) || is.null(library$checkpoints)) {
    stop(
      "checkpoint library must be an msimGD checkpoint library object or .rds path.",
      call. = FALSE
    )
  }

  checkpoints <- library$checkpoints
  if (!is.list(checkpoints) || length(checkpoints) < 1L) {
    stop("checkpoint library has no checkpoints.", call. = FALSE)
  }

  if (!is.null(snapshot_index) && !is.null(timesteps)) {
    stop("Specify at most one of snapshot_index or timesteps.", call. = FALSE)
  }

  idx <- NULL
  if (!is.null(snapshot_index)) {
    snapshot_index <- as.integer(snapshot_index)
    if (length(snapshot_index) != 1L || is.na(snapshot_index) || snapshot_index <= 0L) {
      stop("snapshot_index must be a single integer > 0.", call. = FALSE)
    }
    idx <- snapshot_index
  } else if (!is.null(timesteps)) {
    timesteps <- as.integer(timesteps)
    if (length(timesteps) != 1L || is.na(timesteps) || timesteps <= 0L) {
      stop("timesteps must be a single integer > 0.", call. = FALSE)
    }
    idx <- which(vapply(
      checkpoints,
      function(x) {
        state <- if (!is.null(x$state)) x$state else x
        identical(as.integer(state$timesteps), timesteps)
      },
      logical(1)
    ))
    if (length(idx) != 1L) {
      stop(
        sprintf("checkpoint library has %d matches for timesteps %d.", length(idx), timesteps),
        call. = FALSE
      )
    }
  } else {
    idx <- length(checkpoints)
  }

  if (length(checkpoints) < idx) {
    stop(sprintf("checkpoint library has no snapshot_index %d.", idx), call. = FALSE)
  }

  stationary_checkpoint_library_validate(checkpoints[[idx]])
}

stationary_initializer_snapshot_timesteps <- function(burnin_timesteps,
                                                      n_snapshots = 1L,
                                                      snapshot_spacing = 0L,
                                                      label = "Requested stationary initialization snapshots") {
  burnin_timesteps <- as.integer(burnin_timesteps)
  if (length(burnin_timesteps) != 1L || is.na(burnin_timesteps) || burnin_timesteps <= 0L) {
    stop("burnin_timesteps must be a single integer > 0.", call. = FALSE)
  }

  n_snapshots <- as.integer(n_snapshots)
  if (length(n_snapshots) != 1L || is.na(n_snapshots) || n_snapshots <= 0L) {
    stop("n_snapshots must be a single integer > 0.", call. = FALSE)
  }

  snapshot_spacing <- as.integer(snapshot_spacing)
  if (length(snapshot_spacing) != 1L || is.na(snapshot_spacing) || snapshot_spacing < 0L) {
    stop("snapshot_spacing must be a single integer >= 0.", call. = FALSE)
  }

  if (n_snapshots == 1L) {
    return(burnin_timesteps)
  }

  if (snapshot_spacing <= 0L) {
    stop("snapshot_spacing must be > 0 when n_snapshots > 1.", call. = FALSE)
  }

  first <- burnin_timesteps - snapshot_spacing * (n_snapshots - 1L)
  if (first <= 0L) {
    stop(
      paste(
        label,
        "start before day 1. Increase burnin_timesteps or reduce",
        "n_snapshots / snapshot_spacing."
      ),
      call. = FALSE
    )
  }

  seq.int(from = first, by = snapshot_spacing, length.out = n_snapshots)
}

stationary_human_initializer_human_variable_names <- function(
    parameters,
    include_human_slot_contact_multiplier = TRUE
) {
  common <- c(
    "state",
    "birth",
    "last_boosted_ica",
    "icm",
    "ica",
    "zeta",
    "zeta_group",
    "infectivity",
    "progression_rates",
    "drug",
    "drug_time",
    "last_pev_timestep",
    "last_eff_pev_timestep",
    "pev_profile",
    "tbv_vaccinated",
    "net_time",
    "spray_time"
  )

  if (parameters$parasite == "falciparum") {
    out <- c(
      common,
      "last_boosted_ib",
      "last_boosted_iva",
      "last_boosted_id",
      "ivm",
      "ib",
      "iva",
      "id"
    )
  } else if (parameters$parasite == "vivax") {
    out <- c(
      common,
      "last_boosted_iaa",
      "iaa",
      "iam",
      "hypnozoites"
    )
    if (any(parameters$drug_hypnozoite_efficacy > 0)) {
      out <- c(out, "ls_drug", "ls_drug_time")
    }
  } else {
    stop(sprintf("Unsupported parasite for stochastic human initialization: %s", parameters$parasite), call. = FALSE)
  }

  if (isTRUE(include_human_slot_contact_multiplier)) {
    out <- c(out, "human_slot_contact_multiplier")
  }
  out
}

stationary_human_initializer_state_variable_names <- function(parameters, saved_variables) {
  current <- stationary_human_initializer_human_variable_names(parameters)
  legacy <- stationary_human_initializer_human_variable_names(
    parameters,
    include_human_slot_contact_multiplier = FALSE
  )
  saved_count <- length(saved_variables)
  if (saved_count < length(current) || identical(as.integer(saved_count), as.integer(length(legacy)))) {
    return(legacy)
  }
  birth_idx <- match("birth", legacy)
  size <- if (!is.na(birth_idx) && length(saved_variables) >= birth_idx) {
    length(saved_variables[[birth_idx]])
  } else {
    NA_integer_
  }
  candidate <- saved_variables[[length(current)]]
  candidate_is_current <- is.numeric(candidate) &&
    length(candidate) %in% c(1L, size) &&
    all(is.finite(candidate)) &&
    all(candidate > 0)
  if (!candidate_is_current) {
    return(legacy)
  }
  current
}

stationary_human_initializer_categorical_variables <- function(parameters) {
  out <- list(
    state = c("S", "D", "A", "U", "Tr"),
    zeta_group = to_char_vector(seq_len(parameters$n_heterogeneity_groups))
  )
  out
}

stationary_human_initializer_time_variables <- function(parameters) {
  out <- c(
    "birth",
    "drug_time",
    "last_pev_timestep",
    "last_eff_pev_timestep",
    "tbv_vaccinated"
  )

  if (parameters$parasite == "falciparum") {
    out <- c(out, "last_boosted_ica", "last_boosted_ib", "last_boosted_iva", "last_boosted_id")
  } else if (parameters$parasite == "vivax") {
    out <- c(out, "last_boosted_ica", "last_boosted_iaa")
    if (any(parameters$drug_hypnozoite_efficacy > 0)) {
      out <- c(out, "ls_drug_time")
    }
  }

  out
}

stationary_human_initializer_effective_library_population <- function(parameters) {
  configured <- parameters$human_initialization_library_population
  if (is.null(configured)) {
    return(as.integer(max(1000L, get_human_population(parameters, 0))))
  }

  configured <- as.integer(configured)
  if (length(configured) != 1L || is.na(configured) || configured <= 0L) {
    stop("human_initialization_library_population must be NULL or a single integer > 0.", call. = FALSE)
  }

  configured
}

stationary_human_initializer_parameters <- function(parameters) {
  params <- parameters
  params$human_population <- stationary_human_initializer_effective_library_population(parameters)
  params$human_initialization <- "equilibrium"
  params$human_initialization_library <- NULL
  params$human_initialization_node_index <- NULL
  params$progress_bar <- FALSE
  params$releases <- NULL
  params$releases_schedule <- NULL
  params$mosquito_genotype_history <- NULL
  params$mosquito_aquatic_genotype_history <- NULL
  params$mosquito_infectious_genotype_history <- NULL
  params$cube <- NULL
  params$vector_infectivity_g <- NULL
  params$vector_infectivity_g_by_species <- NULL
  params$debug_genotypes <- FALSE
  params$debug_genotype_timesteps <- NULL

  # Future campaign-style interventions should not influence the baseline
  # library used to construct a stationary day-0 human population.
  params$mda <- FALSE
  params$mda_drug <- 0
  params$mda_timesteps <- NULL
  params$mda_coverages <- NULL
  params$mda_min_ages <- -1
  params$mda_max_ages <- -1
  params$smc <- FALSE
  params$smc_drug <- 0
  params$smc_timesteps <- NULL
  params$smc_coverages <- NULL
  params$smc_min_ages <- -1
  params$smc_max_ages <- -1
  params$pmc <- FALSE
  params$pmc_drug <- 0
  params$pmc_timesteps <- NULL
  params$pmc_coverages <- NULL
  params$pcs_ages <- -1
  params$tbv_timesteps <- NULL
  params$tbv_coverages <- NULL
  params$tbv_ages <- NULL
  params$mass_pev_timesteps <- NULL
  params$mass_pev_coverages <- NULL
  params$mass_pev_min_wait <- NULL
  params$mass_pev_min_ages <- NULL
  params$mass_pev_max_ages <- NULL
  params$mass_pev_booster_spacing <- NULL
  params$mass_pev_booster_coverage <- NULL
  params$mass_pev_booster_profile <- NULL
  params$pev_epi_timesteps <- NULL
  params$pev_epi_coverage <- NULL
  params$pev_epi_waits <- NULL
  params$pev_epi_first_dose_min_ages <- NULL
  params$pev_epi_profile <- NULL
  params$pev_epi_booster_spacing <- NULL
  params$pev_epi_booster_coverage <- NULL
  params$pev_epi_booster_profile <- NULL

  if (!is.null(params$init_EIR)) {
    params <- set_equilibrium(
      params,
      init_EIR = params$init_EIR,
      eq_params = params$eq_params
    )
  }

  params
}

stationary_human_initializer_parameter_context <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.environment(x) || inherits(x, "externalptr") || is.function(x)) {
    return(NULL)
  }

  if (is.list(x)) {
    out <- x
    for (nm in names(out)) {
      out[[nm]] <- stationary_human_initializer_parameter_context(out[[nm]])
    }
    return(out)
  }

  x
}

stationary_human_initializer_reference_timestep <- function(library) {
  available <- length(library$variables$birth)
  reference_timestep <- library$reference_timestep
  if (is.null(reference_timestep)) {
    reference_timestep <- rep.int(as.integer(library$timesteps), available)
  }
  as.integer(reference_timestep)
}

stationary_human_initializer_snapshot_timesteps <- function(parameters) {
  burnin <- parameters$human_initialization_burnin_timesteps
  if (is.null(burnin)) {
    stop(
      paste(
        "stochastic_resample human initialization requires",
        "human_initialization_burnin_timesteps > 0."
      ),
      call. = FALSE
    )
  }

  stationary_initializer_snapshot_timesteps(
    burnin_timesteps = burnin,
    n_snapshots = parameters$human_initialization_n_snapshots,
    snapshot_spacing = parameters$human_initialization_snapshot_spacing,
    label = "Requested stochastic human initialization snapshots"
  )
}

stationary_human_initializer_cache_key <- function(parameters) {
  params <- stationary_human_initializer_parameters(parameters)
  params$human_population_timesteps <- 0L
  params$human_population <- stationary_human_initializer_effective_library_population(parameters)
  params <- stationary_human_initializer_parameter_context(params)
  raw <- serialize(list(
    parameters = params,
    snapshots = stationary_human_initializer_snapshot_timesteps(parameters),
    seed = parameters$human_initialization_seed
  ), NULL)

  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(raw, algo = "xxhash64"))
  }

  path <- tempfile("msimgd-stationary-human-key-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  writeBin(raw, path)
  unname(tools::md5sum(path))
}

stationary_human_initializer_decode_categorical <- function(saved, size) {
  out <- rep(NA_character_, size)
  for (nm in names(saved)) {
    idx <- as.integer(saved[[nm]])
    if (length(idx) > 0L) {
      out[idx] <- nm
    }
  }
  out
}

stationary_human_initializer_extract_library <- function(state, parameters) {
  human_names <- stationary_human_initializer_state_variable_names(
    parameters,
    state$individual$variables
  )
  saved <- state$individual$variables[seq_along(human_names)]
  names(saved) <- human_names
  size <- length(saved$birth)

  decoded <- list()
  categorical <- stationary_human_initializer_categorical_variables(parameters)
  for (nm in names(saved)) {
    if (nm %in% names(categorical)) {
      decoded[[nm]] <- stationary_human_initializer_decode_categorical(saved[[nm]], size)
    } else {
      decoded[[nm]] <- saved[[nm]]
    }
  }

  decoded$state <- as.character(decoded$state)
  decoded$zeta_group <- as.character(decoded$zeta_group)
  if (is.null(decoded$human_slot_contact_multiplier)) {
    decoded$human_slot_contact_multiplier <- resolve_human_slot_contact_multiplier(
      parameters,
      size
    )
  }

  list(
    timesteps = as.integer(state$timesteps),
    reference_timestep = rep.int(as.integer(state$timesteps), size),
    variables = decoded
  )
}

stationary_human_initializer_bind_libraries <- function(libraries) {
  if (length(libraries) < 1L) {
    stop("At least one stationary human initialization snapshot is required.", call. = FALSE)
  }

  variable_names <- names(libraries[[1L]]$variables)
  reference_timestep <- unname(do.call(c, lapply(libraries, function(x) x$reference_timestep)))
  variables <- setNames(vector("list", length(variable_names)), variable_names)

  for (nm in variable_names) {
    variables[[nm]] <- unname(do.call(c, lapply(libraries, function(x) x$variables[[nm]])))
  }

  list(
    timesteps = max(reference_timestep),
    reference_timestep = as.integer(reference_timestep),
    variables = variables
  )
}

stationary_human_initializer_extract_vector_state <- function(state, node_index = NULL) {
  if (!is.list(state) || is.null(state$timesteps) || is.null(state$malariasimulationGD)) {
    return(NULL)
  }

  saved <- state$malariasimulationGD
  if (!is.list(saved) || length(saved) < 7L) {
    return(NULL)
  }

  extract_slot <- function(slot) {
    value <- saved[[slot]]
    if (is.null(node_index)) {
      return(value)
    }
    if (!is.list(value) || length(value) < node_index) {
      return(NULL)
    }
    value[[node_index]]
  }

  list(
    timesteps = as.integer(state$timesteps),
    correlations = extract_slot(2L),
    vector_models = extract_slot(3L),
    solvers = extract_slot(4L),
    lagged_eir = extract_slot(5L),
    lagged_transmission_eir = extract_slot(6L),
    lagged_infectivity = extract_slot(7L)
  )
}

stationary_human_initializer_extract_metapop_library <- function(state, parameters) {
  if (!is.list(parameters) || length(parameters) < 1L) {
    stop("parameters must be a non-empty list of node parameter lists.", call. = FALSE)
  }

  current_block_length <- length(
    stationary_human_initializer_human_variable_names(parameters[[1L]])
  )
  legacy_block_length <- length(
    stationary_human_initializer_human_variable_names(
      parameters[[1L]],
      include_human_slot_contact_multiplier = FALSE
    )
  )
  human_block_length <- current_block_length
  n_nodes <- length(parameters)
  required <- human_block_length * n_nodes
  saved_variables <- state$individual$variables
  if (length(saved_variables) < required) {
    legacy_required <- legacy_block_length * n_nodes
    if (length(saved_variables) >= legacy_required) {
      human_block_length <- legacy_block_length
      required <- legacy_required
    } else {
      stop(
        sprintf(
          paste(
            "Metapop state contains %d saved variables, but %d are required",
            "to extract %d node-conditioned human blocks."
          ),
          length(saved_variables),
          required,
          n_nodes
        ),
        call. = FALSE
      )
    }
  }

  node_libraries <- vector("list", n_nodes)
  for (node in seq_len(n_nodes)) {
    idx <- seq.int(
      from = (node - 1L) * human_block_length + 1L,
      length.out = human_block_length
    )
    node_state <- list(
      timesteps = state$timesteps,
      individual = list(variables = saved_variables[idx])
    )
    node_libraries[[node]] <- stationary_human_initializer_extract_library(
      node_state,
      parameters[[node]]
    )
    node_libraries[[node]]$metadata <- c(
      node_libraries[[node]]$metadata,
      list(vector_state = stationary_human_initializer_extract_vector_state(state, node))
    )
  }

  names(node_libraries) <- to_char_vector(seq_len(n_nodes))
  new_stationary_human_initializer_node_conditioned_library(
    node_libraries = node_libraries,
    metadata = list(timesteps = as.integer(state$timesteps))
  )
}

stationary_human_initializer_extract_metapop_snapshot <- function(state, parameters) {
  if (!is.list(parameters) || length(parameters) < 1L) {
    stop("parameters must be a non-empty list of node parameter lists.", call. = FALSE)
  }

  current_block_length <- length(
    stationary_human_initializer_human_variable_names(parameters[[1L]])
  )
  legacy_block_length <- length(
    stationary_human_initializer_human_variable_names(
      parameters[[1L]],
      include_human_slot_contact_multiplier = FALSE
    )
  )
  human_block_length <- current_block_length
  n_nodes <- length(parameters)
  required <- human_block_length * n_nodes
  saved_variables <- state$individual$variables
  if (length(saved_variables) < required) {
    legacy_required <- legacy_block_length * n_nodes
    if (length(saved_variables) >= legacy_required) {
      human_block_length <- legacy_block_length
      required <- legacy_required
    } else {
      stop(
        sprintf(
          paste(
            "Metapop state contains %d saved variables, but %d are required",
            "to extract %d node human blocks."
          ),
          length(saved_variables),
          required,
          n_nodes
        ),
        call. = FALSE
      )
    }
  }

  node_libraries <- vector("list", n_nodes)
  for (node in seq_len(n_nodes)) {
    idx <- seq.int(
      from = (node - 1L) * human_block_length + 1L,
      length.out = human_block_length
    )
    node_state <- list(
      timesteps = state$timesteps,
      individual = list(variables = saved_variables[idx])
    )
    node_libraries[[node]] <- stationary_human_initializer_extract_library(
      node_state,
      parameters[[node]]
    )
  }

  names(node_libraries) <- to_char_vector(seq_len(n_nodes))
  stationary_human_initializer_validate_metapop_snapshot(
    list(
      node_libraries = node_libraries,
      metadata = list(
        timesteps = as.integer(state$timesteps),
        vector_state = stationary_human_initializer_extract_vector_state(state)
      )
    )
  )
}

stationary_human_initializer_bind_node_conditioned_libraries <- function(libraries) {
  if (length(libraries) < 1L) {
    stop("At least one node-conditioned stationary human library is required.", call. = FALSE)
  }

  if (!all(vapply(libraries, stationary_human_initializer_is_node_conditioned_library, logical(1)))) {
    stop("All libraries must be node-conditioned stationary human libraries.", call. = FALSE)
  }

  first_nodes <- libraries[[1L]]$node_libraries
  node_names <- names(first_nodes)
  if (is.null(node_names)) {
    node_names <- to_char_vector(seq_along(first_nodes))
  }

  node_libraries <- vector("list", length(first_nodes))
  for (i in seq_along(first_nodes)) {
    snapshot_libraries <- lapply(libraries, function(x) x$node_libraries[[i]])
    node_libraries[[i]] <- stationary_human_initializer_bind_libraries(snapshot_libraries)
    node_libraries[[i]]$metadata <- c(
      node_libraries[[i]]$metadata,
      list(snapshot_libraries = snapshot_libraries)
    )
  }
  names(node_libraries) <- node_names

  snapshot_timesteps <- vapply(
    libraries,
    function(x) {
      meta <- x$metadata
      if (!is.null(meta$snapshot_timestep)) {
        return(as.integer(meta$snapshot_timestep))
      }
      if (!is.null(meta$timesteps)) {
        return(as.integer(meta$timesteps))
      }
      NA_integer_
    },
    integer(1)
  )

  new_stationary_human_initializer_node_conditioned_library(
    node_libraries = node_libraries,
    metadata = list(snapshot_timesteps = snapshot_timesteps)
  )
}

stationary_human_initializer_bind_metapop_snapshots <- function(snapshots) {
  if (length(snapshots) < 1L) {
    stop("At least one metapop stationary snapshot is required.", call. = FALSE)
  }

  invisible(lapply(snapshots, stationary_human_initializer_validate_metapop_snapshot))
  snapshot_timesteps <- vapply(
    snapshots,
    function(x) {
      meta <- x$metadata
      if (!is.null(meta$snapshot_timestep)) {
        return(as.integer(meta$snapshot_timestep))
      }
      if (!is.null(meta$timesteps)) {
        return(as.integer(meta$timesteps))
      }
      NA_integer_
    },
    integer(1)
  )

  new_stationary_human_initializer_metapop_library(
    snapshots = snapshots,
    metadata = list(snapshot_timesteps = snapshot_timesteps)
  )
}

stationary_human_initializer_select_metapop_snapshot <- function(library,
                                                                 snapshot_index = NULL,
                                                                 timesteps = NULL) {
  library <- stationary_human_initializer_validate_library(
    stationary_human_initializer_load_library(library)
  )
  if (!stationary_human_initializer_is_metapop_library(library)) {
    stop(
      "library must be a full-metapop stationary human library.",
      call. = FALSE
    )
  }

  snapshots <- library$snapshots
  if (!is.list(snapshots) || length(snapshots) < 1L) {
    stop("full-metapop stationary human library has no snapshots.", call. = FALSE)
  }

  if (!is.null(snapshot_index) && !is.null(timesteps)) {
    stop("Specify at most one of snapshot_index or timesteps.", call. = FALSE)
  }

  idx <- NULL
  if (!is.null(snapshot_index)) {
    snapshot_index <- as.integer(snapshot_index)
    if (length(snapshot_index) != 1L || is.na(snapshot_index) || snapshot_index <= 0L) {
      stop("snapshot_index must be a single integer > 0.", call. = FALSE)
    }
    idx <- snapshot_index
  } else if (!is.null(timesteps)) {
    timesteps <- as.integer(timesteps)
    if (length(timesteps) != 1L || is.na(timesteps) || timesteps <= 0L) {
      stop("timesteps must be a single integer > 0.", call. = FALSE)
    }
    idx <- which(vapply(
      snapshots,
      function(x) {
        meta <- x$metadata
        snapshot_timestep <- if (!is.null(meta$snapshot_timestep)) {
          meta$snapshot_timestep
        } else {
          meta$timesteps
        }
        identical(as.integer(snapshot_timestep), timesteps)
      },
      logical(1)
    ))
    if (length(idx) != 1L) {
      stop(
        sprintf(
          "full-metapop stationary human library has %d matches for timesteps %d.",
          length(idx),
          timesteps
        ),
        call. = FALSE
      )
    }
  } else {
    accepted_field_present <- "accepted_snapshot_indices" %in% names(library$metadata)
    candidate_indices <- library$metadata$accepted_snapshot_indices
    candidate_indices <- unique(as.integer(candidate_indices))
    candidate_indices <- candidate_indices[
      !is.na(candidate_indices) &
        candidate_indices > 0L &
        candidate_indices <= length(snapshots)
    ]
    if (accepted_field_present && length(candidate_indices) < 1L) {
      stop(
        paste(
          "full-metapop stationary human library has no accepted snapshots;",
          "run snapshot screening again or rebuild the candidate library."
        ),
        call. = FALSE
      )
    }
    if (!accepted_field_present) {
      candidate_indices <- seq_along(snapshots)
    }
    idx <- candidate_indices[[sample.int(length(candidate_indices), size = 1L)]]
  }

  if (length(snapshots) < idx) {
    stop(sprintf("full-metapop stationary human library has no snapshot_index %d.", idx),
         call. = FALSE)
  }

  snapshot <- stationary_human_initializer_validate_metapop_snapshot(snapshots[[idx]])
  snapshot$metadata <- c(
    library$metadata,
    snapshot$metadata,
    list(selected_snapshot_index = as.integer(idx))
  )
  snapshot
}

stationary_human_initializer_metapop_node_library <- function(snapshot, node_index) {
  snapshot <- stationary_human_initializer_validate_metapop_snapshot(snapshot)
  node_index <- as.integer(node_index)
  if (length(node_index) != 1L || is.na(node_index) || node_index <= 0L) {
    stop("node_index must be a single integer > 0.", call. = FALSE)
  }

  node_name <- as.character(node_index)
  node_libraries <- snapshot$node_libraries
  node_library <- if (!is.null(names(node_libraries)) && node_name %in% names(node_libraries)) {
    node_libraries[[node_name]]
  } else if (length(node_libraries) >= node_index) {
    node_libraries[[node_index]]
  } else {
    NULL
  }

  if (is.null(node_library)) {
    stop(
      sprintf("metapop snapshot has no node library for node %d.", node_index),
      call. = FALSE
    )
  }

  snapshot_metadata <- snapshot$metadata
  node_library$metadata <- c(
    node_library$metadata,
    list(
      stationary_total_M = stationary_human_initializer_metadata_value(
        snapshot_metadata,
        field_names = c("stationary_total_M_by_node", "snapshot_total_M_by_node"),
        index = node_index
      ),
      stationary_init_foim = stationary_human_initializer_metadata_value(
        snapshot_metadata,
        field_names = c("stationary_init_foim_by_node", "snapshot_init_foim_by_node"),
        index = node_index
      ),
      preserve_human_microstate = TRUE,
      preserve_baseline_parameters = TRUE,
      selected_snapshot_index = snapshot_metadata$selected_snapshot_index,
      snapshot_timestep = .subset2(snapshot_metadata, "snapshot_timestep")
    )
  )
  node_library
}

stationary_human_initializer_prepare_metapop_parameters <- function(parameters) {
  if (!is.list(parameters) || length(parameters) < 1L) {
    return(list(parameters = parameters, context = NULL))
  }

  active <- which(vapply(
    parameters,
    function(p) {
      identical(human_initialization_mode(p), "stochastic_resample") &&
        !is.null(p$human_initialization_library)
    },
    logical(1)
  ))

  if (length(active) < 1L) {
    return(list(parameters = parameters, context = NULL))
  }

  library <- stationary_human_initializer_load_library(
    parameters[[active[[1L]]]]$human_initialization_library
  )
  if (!stationary_human_initializer_is_metapop_library(library)) {
    return(list(parameters = parameters, context = NULL))
  }

  if (length(active) != length(parameters)) {
    stop(
      paste(
        "A full-metapop stationary human library requires every metapop node",
        "to use the shared fresh-start initializer."
      ),
      call. = FALSE
    )
  }

  selected_snapshot <- stationary_human_initializer_select_metapop_snapshot(library)
  if (length(selected_snapshot$node_libraries) != length(parameters)) {
    stop(
      sprintf(
        paste(
          "Selected metapop stationary snapshot has %d node blocks, but the",
          "simulation has %d nodes."
        ),
        length(selected_snapshot$node_libraries),
        length(parameters)
      ),
      call. = FALSE
    )
  }

  prepared_parameters <- Map(
    function(p, node_index) {
      p$human_initialization_library <- stationary_human_initializer_metapop_node_library(
        selected_snapshot,
        node_index = node_index
      )
      p$human_initialization_node_index <- NULL
      p
    },
    parameters,
    seq_along(parameters)
  )

  list(
    parameters = prepared_parameters,
    context = list(
      vector_state = selected_snapshot$metadata$vector_state,
      selected_snapshot_index = selected_snapshot$metadata$selected_snapshot_index,
      snapshot_timestep = .subset2(selected_snapshot$metadata, "snapshot_timestep")
    )
  )
}

stationary_human_initializer_with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }

  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(code)
}

get_stationary_human_initialization_library <- function(parameters) {
  if (!is.null(parameters$human_initialization_library)) {
    return(stationary_human_initializer_select_library(
      parameters$human_initialization_library,
      parameters
    ))
  }

  key <- stationary_human_initializer_cache_key(parameters)
  if (exists(key, envir = .stationary_human_initializer_cache, inherits = FALSE)) {
    return(get(key, envir = .stationary_human_initializer_cache, inherits = FALSE))
  }

  aux_parameters <- stationary_human_initializer_parameters(parameters)
  snapshot_timesteps <- stationary_human_initializer_snapshot_timesteps(parameters)
  seed <- parameters$human_initialization_seed
  library <- stationary_human_initializer_with_seed(
    seed,
    {
      current_state <- NULL
      snapshots <- vector("list", length(snapshot_timesteps))

      for (i in seq_along(snapshot_timesteps)) {
        sim <- run_resumable_simulation(
          timesteps = snapshot_timesteps[[i]],
          parameters = aux_parameters,
          initial_state = current_state,
          restore_random_state = !is.null(current_state)
        )
        current_state <- sim$state
        snapshots[[i]] <- stationary_human_initializer_extract_library(
          current_state,
          aux_parameters
        )
      }

      stationary_human_initializer_bind_libraries(snapshots)
    }
  )

  assign(key, library, envir = .stationary_human_initializer_cache)
  library
}

stationary_human_initializer_resample <- function(library, size, parameters) {
  snapshot_libraries <- library$metadata$snapshot_libraries
  sampled_vector_state <- NULL
  if (is.list(snapshot_libraries) && length(snapshot_libraries) > 0L) {
    snapshot_idx <- sample.int(length(snapshot_libraries), size = 1L)
    library <- snapshot_libraries[[snapshot_idx]]
    sampled_vector_state <- library$metadata$vector_state
  }

  available <- length(library$variables$birth)
  if (available < 1L) {
    stop("Stationary human initialization library is empty.", call. = FALSE)
  }

  if (isTRUE(library$metadata$preserve_human_microstate)) {
    if (size != available) {
      stop(
        sprintf(
          paste(
            "Exact stationary human initialization requires size %d, but",
            "requested human_population is %d."
          ),
          available,
          size
        ),
        call. = FALSE
      )
    }
    index <- seq_len(available)
    sampled <- lapply(library$variables, identity)
    reference_timestep <- stationary_human_initializer_reference_timestep(library)
  } else {
    replace <- size > available
    index <- sample.int(available, size = size, replace = replace)
    sampled <- lapply(library$variables, function(x) x[index])
    reference_timestep <- stationary_human_initializer_reference_timestep(library)[index]
  }

  shift_names <- stationary_human_initializer_time_variables(parameters)
  for (nm in intersect(names(sampled), shift_names)) {
    shifted <- sampled[[nm]]
    keep <- shifted != -1
    shifted[keep] <- shifted[keep] - reference_timestep[keep]
    sampled[[nm]] <- shifted
  }
  if (is.null(sampled$human_slot_contact_multiplier)) {
    sampled$human_slot_contact_multiplier <- resolve_human_slot_contact_multiplier(
      parameters,
      size
    )
  } else {
    sampled$human_slot_contact_multiplier <- validate_human_slot_contact_multiplier(
      sampled$human_slot_contact_multiplier,
      size
    )
  }

  attr(sampled, "stationary_vector_state") <- sampled_vector_state
  unique_reference_timestep <- unique(as.integer(reference_timestep))
  unique_reference_timestep <- unique_reference_timestep[!is.na(unique_reference_timestep)]
  if (length(unique_reference_timestep) == 1L) {
    attr(sampled, "stationary_vector_control_time_offset") <- unique_reference_timestep
  }
  sampled
}

stationary_human_initializer_current_foim <- function(
    birth,
    zeta,
    human_slot_contact_multiplier = NULL,
    net_time,
    spray_time,
    infectivity,
    parameters
) {
  birth <- individual::IntegerVariable$new(as.integer(birth))
  zeta <- individual::DoubleVariable$new(as.numeric(zeta))
  human_slot_contact_multiplier <- if (is.null(human_slot_contact_multiplier)) {
    NULL
  } else {
    individual::DoubleVariable$new(
      validate_human_slot_contact_multiplier(
        human_slot_contact_multiplier,
        length(zeta$get_values())
      )
    )
  }
  net_time <- individual::IntegerVariable$new(as.integer(net_time))
  spray_time <- individual::IntegerVariable$new(as.integer(spray_time))
  age <- get_age(birth$get_values(), 0)
  psi <- unique_biting_rate(age, parameters)
  human_slot_contact_values <- if (is.null(human_slot_contact_multiplier)) {
    NULL
  } else {
    human_slot_contact_multiplier$get_values()
  }
  .pi <- human_pi(zeta$get_values(), psi, human_slot_contact_values)
  current_infectivity <- sum(.pi * as.numeric(infectivity))
  if (!is.finite(current_infectivity) || current_infectivity < 0) {
    return(NULL)
  }

  a0 <- human_blood_meal_rate(
    1L,
    list(
      birth = birth,
      zeta = zeta,
      human_slot_contact_multiplier = human_slot_contact_multiplier,
      net_time = net_time,
      spray_time = spray_time
    ),
    parameters,
    0
  )
  if (!is.finite(a0) || a0 < 0) {
    return(NULL)
  }

  a0 * current_infectivity
}

stationary_human_initializer_library_foim <- function(library, parameters) {
  stationary_human_initializer_current_foim(
    birth = library$variables$birth,
    zeta = library$variables$zeta,
    human_slot_contact_multiplier = library$variables$human_slot_contact_multiplier,
    net_time = library$variables$net_time,
    spray_time = library$variables$spray_time,
    infectivity = library$variables$infectivity,
    parameters = parameters
  )
}

stationary_human_initializer_sample_context <- function(sampled, parameters, library = NULL) {
  metadata <- if (is.list(library)) library$metadata else NULL
  preserve_baseline_parameters <- isTRUE(.subset2(metadata, "preserve_baseline_parameters"))
  vector_control_time_offset <- .subset2(metadata, "snapshot_timestep")
  if (!is.numeric(vector_control_time_offset) ||
      length(vector_control_time_offset) != 1L ||
      !is.finite(vector_control_time_offset) ||
      vector_control_time_offset < 0) {
    vector_control_time_offset <- attr(sampled, "stationary_vector_control_time_offset")
  }
  parameters_for_context <- parameters
  if (is.numeric(vector_control_time_offset) &&
      length(vector_control_time_offset) == 1L &&
      is.finite(vector_control_time_offset) &&
      vector_control_time_offset >= 0) {
    parameters_for_context$vector_control_time_offset <- as.integer(vector_control_time_offset)
  }

  list(
    init_foim = if (preserve_baseline_parameters) {
      NULL
    } else {
      stationary_human_initializer_current_foim(
        birth = sampled$birth,
        zeta = sampled$zeta,
        human_slot_contact_multiplier = sampled$human_slot_contact_multiplier,
        net_time = sampled$net_time,
        spray_time = sampled$spray_time,
        infectivity = sampled$infectivity,
        parameters = parameters_for_context
      )
    },
    total_M = if (preserve_baseline_parameters) {
      NULL
    } else {
      stationary_human_initializer_metadata_value(
        metadata,
        field_names = c(
          "stationary_total_M",
          "baseline_total_M",
          "parameter_total_M",
          "effective_total_M"
        )
      )
    },
    preserve_baseline_parameters = preserve_baseline_parameters,
    vector_control_time_offset = vector_control_time_offset,
    vector_state = attr(sampled, "stationary_vector_state")
  )
}

stationary_human_initializer_apply_context <- function(parameters, context) {
  if (is.null(context) || !is.list(context)) {
    return(parameters)
  }

  if (!isTRUE(context$preserve_baseline_parameters)) {
    total_M <- context$total_M
    if (is.numeric(total_M) && length(total_M) == 1L && is.finite(total_M) && total_M > 0) {
      parameters <- parameterise_total_M(parameters, as.numeric(total_M))
    }

    init_foim <- context$init_foim
    if (is.numeric(init_foim) && length(init_foim) == 1L && is.finite(init_foim) && init_foim >= 0) {
      parameters$init_foim <- as.numeric(init_foim)
    }
  }

  vector_control_time_offset <- context$vector_control_time_offset
  if (is.numeric(vector_control_time_offset) &&
      length(vector_control_time_offset) == 1L &&
      is.finite(vector_control_time_offset) &&
      vector_control_time_offset >= 0) {
    parameters$vector_control_time_offset <- as.integer(vector_control_time_offset)
  }

  parameters
}

stationary_human_initializer_rebase_runtime_state <- function(state, offset) {
  if (is.null(state) || !is.numeric(offset) || length(offset) != 1L ||
      !is.finite(offset) || offset == 0) {
    return(state)
  }

  if (is.data.frame(state) && "timestep" %in% names(state)) {
    state$timestep <- as.numeric(state$timestep) - offset
    return(state)
  }

  if (!is.list(state)) {
    return(state)
  }

  lapply(state, stationary_human_initializer_rebase_runtime_state, offset = offset)
}

stationary_human_initializer_rebase_vector_state <- function(vector_state) {
  if (is.null(vector_state) || !is.list(vector_state) || is.null(vector_state$timesteps)) {
    return(vector_state)
  }

  offset <- as.numeric(vector_state$timesteps)
  rebased <- vector_state
  rebased$timesteps <- 0L
  rebased$correlations <- stationary_human_initializer_rebase_runtime_state(
    rebased$correlations,
    offset
  )
  rebased$lagged_eir <- stationary_human_initializer_rebase_runtime_state(
    rebased$lagged_eir,
    offset
  )
  rebased$lagged_transmission_eir <- stationary_human_initializer_rebase_runtime_state(
    rebased$lagged_transmission_eir,
    offset
  )
  rebased$lagged_infectivity <- stationary_human_initializer_rebase_runtime_state(
    rebased$lagged_infectivity,
    offset
  )
  rebased
}

stationary_human_initializer_restore_vector_state <- function(
    context,
    correlations = NULL,
    vector_models,
    solvers,
    lagged_eir,
    lagged_transmission_eir,
    lagged_infectivity
) {
  if (is.null(context) || !is.list(context) || is.null(context$vector_state)) {
    return(invisible(NULL))
  }

  vector_state <- context$vector_state
  if (!is.list(vector_state) || is.null(vector_state$timesteps)) {
    return(invisible(NULL))
  }

  vector_state <- stationary_human_initializer_rebase_vector_state(vector_state)
  individual::restore_object_state(
    as.integer(vector_state$timesteps),
    list(
      correlations,
      vector_models,
      solvers,
      lagged_eir,
      lagged_transmission_eir,
      lagged_infectivity
    ),
    list(
      vector_state$correlations,
      vector_state$vector_models,
      vector_state$solvers,
      vector_state$lagged_eir,
      vector_state$lagged_transmission_eir,
      vector_state$lagged_infectivity
    )
  )

  invisible(NULL)
}
