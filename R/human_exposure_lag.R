HumanExposureLagBuffer <- R6::R6Class(
  "HumanExposureLagBuffer",
  private = list(
    max_rows = NULL,
    n_rows = 0L,
    next_row = 1L,
    default_exposure = NULL,
    default_weighted_exposure = NULL,
    timesteps = numeric(0),
    exposure = NULL,
    weighted_exposure = NULL,
    default_species_exposure = NULL,
    default_species_weighted_exposure = NULL,
    species_exposure = NULL,
    species_weighted_exposure = NULL,

    ordered_rows = function() {
      if (private$n_rows == 0L) {
        return(integer(0))
      }
      if (private$n_rows < private$max_rows) {
        return(seq_len(private$n_rows))
      }
      c(seq.int(private$next_row, private$max_rows), seq_len(private$next_row - 1L))
    },

    ordered_timesteps = function() {
      private$timesteps[private$ordered_rows()]
    },

    ordered_matrix = function(values) {
      if (is.null(values)) {
        return(NULL)
      }
      values[private$ordered_rows(), , drop = FALSE]
    },

    ordered_array = function(values) {
      if (is.null(values)) {
        return(NULL)
      }
      values[private$ordered_rows(), , , drop = FALSE]
    },

    check_values = function(values, label) {
      values <- as.numeric(values)
      if (length(values) != length(private$default_exposure) ||
          anyNA(values) || any(!is.finite(values))) {
        stop(sprintf("%s must be a finite numeric vector matching the human population.", label), call. = FALSE)
      }
      values
    },

    select_history = function(weighted) {
      if (isTRUE(weighted)) {
        return(private$weighted_exposure)
      }
      private$exposure
    },

    select_default = function(weighted) {
      if (isTRUE(weighted)) {
        return(private$default_weighted_exposure)
      }
      private$default_exposure
    },

    select_species_history = function(weighted) {
      if (isTRUE(weighted)) {
        return(private$species_weighted_exposure)
      }
      private$species_exposure
    },

    select_species_default = function(weighted) {
      if (isTRUE(weighted)) {
        return(private$default_species_weighted_exposure)
      }
      private$default_species_exposure
    },

    check_species_values = function(values, label) {
      values <- as.matrix(values)
      default <- private$default_species_exposure
      if (length(dim(values)) != 2L ||
          !identical(dim(values), dim(default)) ||
          anyNA(values) ||
          any(!is.finite(values))) {
        stop(sprintf("%s must be a finite numeric matrix matching the species by human exposure shape.", label), call. = FALSE)
      }
      values
    }
  ),
  public = list(
    initialize = function(
        max_lag,
        default_exposure,
        default_weighted_exposure = NULL,
        default_species_exposure = NULL,
        default_species_weighted_exposure = NULL
    ) {
      max_rows <- as.integer(ceiling(max_lag) + 2L)
      if (max_rows < 2L) {
        max_rows <- 2L
      }
      default_exposure <- as.numeric(default_exposure)
      if (length(default_exposure) < 1L || anyNA(default_exposure) || any(!is.finite(default_exposure))) {
        stop("default_exposure must be a finite numeric vector.", call. = FALSE)
      }

      private$max_rows <- max_rows
      private$n_rows <- 0L
      private$next_row <- 1L
      private$default_exposure <- default_exposure
      private$timesteps <- rep(NA_real_, max_rows)
      private$exposure <- matrix(0, nrow = max_rows, ncol = length(default_exposure))

      if (!is.null(default_species_exposure)) {
        default_species_exposure <- as.matrix(default_species_exposure)
        if (ncol(default_species_exposure) != length(default_exposure) ||
            anyNA(default_species_exposure) ||
            any(!is.finite(default_species_exposure))) {
          stop("default_species_exposure must be a finite numeric species by human matrix.", call. = FALSE)
        }
        private$default_species_exposure <- default_species_exposure
        private$species_exposure <- array(
          0,
          dim = c(max_rows, nrow(default_species_exposure), ncol(default_species_exposure))
        )
      }

      if (!is.null(default_weighted_exposure)) {
        default_weighted_exposure <- as.numeric(default_weighted_exposure)
        if (length(default_weighted_exposure) != length(default_exposure) ||
            anyNA(default_weighted_exposure) || any(!is.finite(default_weighted_exposure))) {
          stop("default_weighted_exposure must be a finite numeric vector matching default_exposure.", call. = FALSE)
        }
        private$default_weighted_exposure <- default_weighted_exposure
        private$weighted_exposure <- matrix(0, nrow = max_rows, ncol = length(default_exposure))
        if (!is.null(default_species_weighted_exposure)) {
          default_species_weighted_exposure <- as.matrix(default_species_weighted_exposure)
          if (is.null(private$default_species_exposure) ||
              !identical(dim(default_species_weighted_exposure), dim(private$default_species_exposure)) ||
              anyNA(default_species_weighted_exposure) ||
              any(!is.finite(default_species_weighted_exposure))) {
            stop("default_species_weighted_exposure must be a finite numeric matrix matching default_species_exposure.", call. = FALSE)
          }
          private$default_species_weighted_exposure <- default_species_weighted_exposure
          private$species_weighted_exposure <- array(
            0,
            dim = dim(private$species_exposure)
          )
        }
      }
    },

    save = function(
        timestep,
        exposure,
        weighted_exposure = NULL,
        species_exposure = NULL,
        species_weighted_exposure = NULL
    ) {
      exposure <- private$check_values(exposure, "exposure")
      if (!is.null(private$weighted_exposure)) {
        weighted_exposure <- private$check_values(weighted_exposure, "weighted_exposure")
      }
      if (!is.null(private$species_exposure)) {
        if (is.null(species_exposure) && nrow(private$default_species_exposure) == 1L) {
          species_exposure <- matrix(exposure, nrow = 1L)
        }
        species_exposure <- private$check_species_values(species_exposure, "species_exposure")
      }
      if (!is.null(private$species_weighted_exposure)) {
        if (is.null(species_weighted_exposure) && nrow(private$default_species_weighted_exposure) == 1L) {
          species_weighted_exposure <- matrix(weighted_exposure, nrow = 1L)
        }
        species_weighted_exposure <- private$check_species_values(species_weighted_exposure, "species_weighted_exposure")
      }

      rows <- private$ordered_rows()
      existing <- which(private$timesteps[rows] == timestep)
      if (length(existing) > 0L) {
        row <- rows[[existing[[1L]]]]
        private$exposure[row, ] <- exposure
        if (!is.null(private$weighted_exposure)) {
          private$weighted_exposure[row, ] <- weighted_exposure
        }
        if (!is.null(private$species_exposure)) {
          private$species_exposure[row, , ] <- species_exposure
        }
        if (!is.null(private$species_weighted_exposure)) {
          private$species_weighted_exposure[row, , ] <- species_weighted_exposure
        }
        return(invisible(NULL))
      }

      row <- private$next_row
      private$timesteps[[row]] <- as.numeric(timestep)
      private$exposure[row, ] <- exposure
      if (!is.null(private$weighted_exposure)) {
        private$weighted_exposure[row, ] <- weighted_exposure
      }
      if (!is.null(private$species_exposure)) {
        private$species_exposure[row, , ] <- species_exposure
      }
      if (!is.null(private$species_weighted_exposure)) {
        private$species_weighted_exposure[row, , ] <- species_weighted_exposure
      }
      if (private$n_rows < private$max_rows) {
        private$n_rows <- private$n_rows + 1L
      }
      private$next_row <- if (row == private$max_rows) 1L else row + 1L

      invisible(NULL)
    },

    get = function(timestep, weighted = FALSE, by_species = FALSE) {
      timesteps <- private$ordered_timesteps()
      if (isTRUE(by_species)) {
        history <- private$ordered_array(private$select_species_history(weighted))
        default <- private$select_species_default(weighted)
        if (is.null(history) || is.null(default)) {
          return(NULL)
        }
        if (length(timesteps) == 0L || timestep < timesteps[[1L]]) {
          return(default)
        }

        exact <- which(timesteps == timestep)
        if (length(exact) > 0L) {
          return(array_row_matrix(history, exact[[1L]]))
        }

        if (timestep > timesteps[[length(timesteps)]]) {
          warning("Exposure lag lookup is after the latest saved timestep; returning per-human defaults.", call. = FALSE)
          return(default)
        }

        after <- which(timesteps > timestep)[[1L]]
        before <- after - 1L
        weight <- (timestep - timesteps[[before]]) /
          (timesteps[[after]] - timesteps[[before]])

        before_values <- array_row_matrix(history, before)
        after_values <- array_row_matrix(history, after)
        return(before_values + weight * (after_values - before_values))
      }

      history <- private$ordered_matrix(private$select_history(weighted))
      default <- private$select_default(weighted)
      if (is.null(history) || is.null(default)) {
        return(NULL)
      }
      if (length(timesteps) == 0L || timestep < timesteps[[1L]]) {
        return(default)
      }

      exact <- which(timesteps == timestep)
      if (length(exact) > 0L) {
        return(as.numeric(history[exact[[1L]], ]))
      }

      if (timestep > timesteps[[length(timesteps)]]) {
        warning("Exposure lag lookup is after the latest saved timestep; returning per-human defaults.", call. = FALSE)
        return(default)
      }

      after <- which(timesteps > timestep)[[1L]]
      before <- after - 1L
      weight <- (timestep - timesteps[[before]]) /
        (timesteps[[after]] - timesteps[[before]])

      as.numeric(history[before, ] + weight * (history[after, ] - history[before, ]))
    },

    clear = function(target) {
      index <- if (inherits(target, "Bitset")) target$to_vector() else as.integer(target)
      if (length(index) == 0L) {
        return(invisible(NULL))
      }
      rows <- private$ordered_rows()
      if (length(rows) > 0L) {
        private$exposure[rows, index] <- matrix(
          private$default_exposure[index],
          nrow = length(rows),
          ncol = length(index),
          byrow = TRUE
        )
      }
      if (!is.null(private$weighted_exposure)) {
        if (length(rows) > 0L) {
          private$weighted_exposure[rows, index] <- matrix(
            private$default_weighted_exposure[index],
            nrow = length(rows),
            ncol = length(index),
            byrow = TRUE
          )
        }
      }
      if (!is.null(private$species_exposure) && length(rows) > 0L) {
        for (row in rows) {
          private$species_exposure[row, , index] <- private$default_species_exposure[, index, drop = FALSE]
        }
      }
      if (!is.null(private$species_weighted_exposure) && length(rows) > 0L) {
        for (row in rows) {
          private$species_weighted_exposure[row, , index] <- private$default_species_weighted_exposure[, index, drop = FALSE]
        }
      }
      invisible(NULL)
    },

    save_state = function() {
      list(
        timesteps = private$ordered_timesteps(),
        exposure = private$ordered_matrix(private$exposure),
        weighted_exposure = private$ordered_matrix(private$weighted_exposure),
        default_exposure = private$default_exposure,
        default_weighted_exposure = private$default_weighted_exposure,
        species_exposure = private$ordered_array(private$species_exposure),
        species_weighted_exposure = private$ordered_array(private$species_weighted_exposure),
        default_species_exposure = private$default_species_exposure,
        default_species_weighted_exposure = private$default_species_weighted_exposure
      )
    },

    restore_state = function(timestep, state) {
      if (is.null(state)) {
        return(invisible(NULL))
      }
      timesteps <- as.numeric(state$timesteps)
      exposure <- as.matrix(state$exposure)
      n <- min(length(timesteps), private$max_rows)
      if (n > 0L && length(timesteps) > n) {
        keep <- seq.int(length(timesteps) - n + 1L, length(timesteps))
        timesteps <- timesteps[keep]
        exposure <- exposure[keep, , drop = FALSE]
      }
      private$timesteps <- rep(NA_real_, private$max_rows)
      private$exposure <- matrix(0, nrow = private$max_rows, ncol = length(private$default_exposure))
      private$n_rows <- n
      private$next_row <- if (n == private$max_rows) 1L else n + 1L
      if (n > 0L) {
        private$timesteps[seq_len(n)] <- timesteps
        private$exposure[seq_len(n), ] <- exposure
      }
      private$default_exposure <- as.numeric(state$default_exposure)

      if (!is.null(private$weighted_exposure)) {
        private$weighted_exposure <- matrix(0, nrow = private$max_rows, ncol = length(private$default_exposure))
        if (n > 0L) {
          weighted_exposure <- as.matrix(state$weighted_exposure)
          if (nrow(weighted_exposure) > n) {
            weighted_exposure <- weighted_exposure[seq.int(nrow(weighted_exposure) - n + 1L, nrow(weighted_exposure)), , drop = FALSE]
          }
          private$weighted_exposure[seq_len(n), ] <- weighted_exposure
        }
        private$default_weighted_exposure <- as.numeric(state$default_weighted_exposure)
      }
      if (!is.null(state$species_exposure)) {
        species_exposure <- state$species_exposure
        if (dim(species_exposure)[[1L]] > n) {
          keep <- seq.int(dim(species_exposure)[[1L]] - n + 1L, dim(species_exposure)[[1L]])
          species_exposure <- species_exposure[keep, , , drop = FALSE]
        }
        private$default_species_exposure <- as.matrix(state$default_species_exposure)
        private$species_exposure <- array(
          0,
          dim = c(private$max_rows, nrow(private$default_species_exposure), ncol(private$default_species_exposure))
        )
        if (n > 0L) {
          private$species_exposure[seq_len(n), , ] <- species_exposure
        }
      } else if (!is.null(private$species_exposure) &&
          nrow(private$default_species_exposure) == 1L &&
          n > 0L) {
        private$species_exposure <- array(
          0,
          dim = c(private$max_rows, 1L, ncol(private$exposure))
        )
        private$species_exposure[seq_len(n), , ] <- private$exposure[seq_len(n), , drop = FALSE]
      }
      if (!is.null(state$species_weighted_exposure)) {
        species_weighted_exposure <- state$species_weighted_exposure
        if (dim(species_weighted_exposure)[[1L]] > n) {
          keep <- seq.int(dim(species_weighted_exposure)[[1L]] - n + 1L, dim(species_weighted_exposure)[[1L]])
          species_weighted_exposure <- species_weighted_exposure[keep, , , drop = FALSE]
        }
        private$default_species_weighted_exposure <- as.matrix(state$default_species_weighted_exposure)
        private$species_weighted_exposure <- array(
          0,
          dim = c(private$max_rows, nrow(private$default_species_weighted_exposure), ncol(private$default_species_weighted_exposure))
        )
        if (n > 0L) {
          private$species_weighted_exposure[seq_len(n), , ] <- species_weighted_exposure
        }
      } else if (!is.null(private$species_weighted_exposure) &&
          nrow(private$default_species_weighted_exposure) == 1L &&
          n > 0L) {
        private$species_weighted_exposure <- array(
          0,
          dim = c(private$max_rows, 1L, ncol(private$weighted_exposure))
        )
        private$species_weighted_exposure[seq_len(n), , ] <- private$weighted_exposure[seq_len(n), , drop = FALSE]
      }
      invisible(NULL)
    }
  )
)

abind_first_dimension <- function(array_values, matrix_values) {
  matrix_values <- as.matrix(matrix_values)
  if (dim(array_values)[[1L]] == 0L) {
    return(array(matrix_values, dim = c(1L, nrow(matrix_values), ncol(matrix_values))))
  }
  out <- array(
    numeric((dim(array_values)[[1L]] + 1L) * nrow(matrix_values) * ncol(matrix_values)),
    dim = c(dim(array_values)[[1L]] + 1L, nrow(matrix_values), ncol(matrix_values))
  )
  out[seq_len(dim(array_values)[[1L]]), , ] <- array_values
  out[dim(out)[[1L]], , ] <- matrix_values
  out
}

array_row_matrix <- function(array_values, row) {
  matrix(
    array_values[row, , ],
    nrow = dim(array_values)[[2L]],
    ncol = dim(array_values)[[3L]]
  )
}

human_exposure_lag_weighted_active <- function(parameters) {
  any(vapply(
    parameters,
    function(p) any(vapply(p$species, function(s) !is.null(human_transmission_weights_for_species(p, s)), logical(1))),
    logical(1)
  ))
}

human_exposure_lag_node_defaults <- function(
    parameters,
    variables,
    solvers,
    weighted = FALSE,
    lagged_values = NULL
) {
  vapply(
    seq_along(parameters),
    function(i) {
      if (!is.null(lagged_values)) {
        return(sum(vapply(
          lagged_values[[i]],
          function(lagged) lagged$get(-.Machine$double.xmax),
          numeric(1)
        )))
      }

      values <- vapply(
        seq_along(parameters[[i]]$species),
        function(species) {
          if (isTRUE(weighted)) {
            return(calculate_transmission_eir(species, solvers[[i]], variables[[i]], parameters[[i]], 0))
          }
          calculate_eir(species, solvers[[i]], variables[[i]], parameters[[i]], 0)
        },
        numeric(1)
      )
      sum(values)
    },
    numeric(1)
  )
}

human_exposure_lag_node_species_defaults <- function(
    parameters,
    variables,
    solvers,
    weighted = FALSE,
    lagged_values = NULL
) {
  lapply(
    seq_along(parameters),
    function(i) {
      if (!is.null(lagged_values)) {
        return(vapply(
          lagged_values[[i]],
          function(lagged) lagged$get(-.Machine$double.xmax),
          numeric(1)
        ))
      }

      vapply(
        seq_along(parameters[[i]]$species),
        function(species) {
          if (isTRUE(weighted)) {
            return(calculate_transmission_eir(species, solvers[[i]], variables[[i]], parameters[[i]], 0))
          }
          calculate_eir(species, solvers[[i]], variables[[i]], parameters[[i]], 0)
        },
        numeric(1)
      )
    }
  )
}

human_exposure_lag_validate_node_exposure <- function(values, expected_length, label) {
  values <- as.numeric(values)
  if (length(values) != expected_length ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values < 0)) {
    stop(sprintf("%s must be a nonnegative finite numeric vector matching the number of mosquito species.", label), call. = FALSE)
  }
  values
}

human_exposure_lag_validate_current_node <- function(current_node, n_nodes) {
  current_node <- as.integer(current_node)
  if (anyNA(current_node) || any(current_node < 1L) || any(current_node > n_nodes)) {
    stop("current_node must contain valid node indices for human exposure allocation.", call. = FALSE)
  }
  current_node
}

human_exposure_lag_biting_allocation_weights <- function(variables, parameters, timestep, species) {
  age <- get_age(variables$birth$get_values(), timestep)
  psi <- unique_biting_rate(age, parameters)
  raw_weights <- human_biting_weights(
    variables$zeta$get_values(),
    psi,
    human_slot_contact_multiplier_values(variables)
  )
  bite_probability <- prob_bitten(timestep, variables, species, parameters)$prob_bitten
  weights <- raw_weights * bite_probability
  if (length(weights) != parameters$human_population ||
      anyNA(weights) ||
      any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("Human biting allocation weights must be nonnegative and finite.", call. = FALSE)
  }
  weights
}

human_exposure_lag_present_psi <- function(context, destination_node, timestep) {
  present_psi <- numeric(0)
  for (home_node in seq_len(context$n_nodes)) {
    current_node <- human_exposure_lag_validate_current_node(
      context$variables[[home_node]]$current_node$get_values(),
      context$n_nodes
    )
    present <- current_node == destination_node
    if (!any(present)) {
      next
    }
    age <- get_age(context$variables[[home_node]]$birth$get_values(), timestep)
    present_psi <- c(present_psi, unique_biting_rate(age, context$parameters[[home_node]])[present])
  }
  present_psi
}

human_exposure_lag_expected_bite_count <- function(context, destination_node, timestep, eir) {
  present_psi <- human_exposure_lag_present_psi(context, destination_node, timestep)
  if (length(present_psi) == 0L) {
    return(0)
  }
  legacy_expected_infectious_bites(
    species_eir = eir,
    psi = present_psi,
    parameters = context$parameters[[destination_node]]
  )
}

human_exposure_lag_destination_shares <- function(context, destination_node, species, timestep) {
  shares <- vector("list", context$n_nodes)
  total_weight <- 0

  for (home_node in seq_len(context$n_nodes)) {
    current_node <- human_exposure_lag_validate_current_node(
      context$variables[[home_node]]$current_node$get_values(),
      context$n_nodes
    )
    weights <- human_exposure_lag_biting_allocation_weights(
      context$variables[[home_node]],
      context$parameters[[home_node]],
      timestep,
      species
    )
    weights[current_node != destination_node] <- 0
    shares[[home_node]] <- weights
    total_weight <- total_weight + sum(weights)
  }

  if (!is.finite(total_weight) || total_weight <= 0) {
    return(lapply(shares, function(x) rep(0, length(x))))
  }

  lapply(shares, function(x) x / total_weight)
}

human_exposure_lag_allocate_exposure <- function(context, timestep, node_exposure, node_weighted_exposure = NULL) {
  exposure_by_home <- lapply(
    context$parameters,
    function(parameters) rep(0, parameters$human_population)
  )
  species_by_home <- lapply(
    context$parameters,
    function(parameters) matrix(0, nrow = length(parameters$species), ncol = parameters$human_population)
  )
  weighted_by_home <- if (isTRUE(context$weighted_active)) {
    lapply(context$parameters, function(parameters) rep(0, parameters$human_population))
  } else {
    NULL
  }
  species_weighted_by_home <- if (isTRUE(context$weighted_active)) {
    lapply(
      context$parameters,
      function(parameters) matrix(0, nrow = length(parameters$species), ncol = parameters$human_population)
    )
  } else {
    NULL
  }

  for (destination_node in seq_len(context$n_nodes)) {
    species_exposure <- node_exposure[[destination_node]]
    species_weighted_exposure <- if (isTRUE(context$weighted_active)) {
      node_weighted_exposure[[destination_node]]
    } else {
      NULL
    }

    for (species in seq_along(species_exposure)) {
      expected_bites <- human_exposure_lag_expected_bite_count(
        context,
        destination_node,
        timestep,
        species_exposure[[species]]
      )
      weighted_expected_bites <- if (isTRUE(context$weighted_active)) {
        human_exposure_lag_expected_bite_count(
          context,
          destination_node,
          timestep,
          species_weighted_exposure[[species]]
        )
      } else {
        NULL
      }
      shares <- human_exposure_lag_destination_shares(
        context,
        destination_node,
        species,
        timestep
      )
      for (home_node in seq_len(context$n_nodes)) {
        exposure_by_home[[home_node]] <- exposure_by_home[[home_node]] +
          expected_bites * shares[[home_node]]
        species_by_home[[home_node]][species, ] <- species_by_home[[home_node]][species, ] +
          expected_bites * shares[[home_node]]
        if (isTRUE(context$weighted_active)) {
          weighted_by_home[[home_node]] <- weighted_by_home[[home_node]] +
            weighted_expected_bites * shares[[home_node]]
          species_weighted_by_home[[home_node]][species, ] <- species_weighted_by_home[[home_node]][species, ] +
            weighted_expected_bites * shares[[home_node]]
        }
      }
    }
  }

  list(
    exposure = exposure_by_home,
    weighted_exposure = weighted_by_home,
    species_exposure = species_by_home,
    species_weighted_exposure = species_weighted_by_home
  )
}

create_human_exposure_lag_context <- function(
    parameters,
    variables,
    solvers,
    lagged_eir = NULL,
    lagged_transmission_eir = NULL
) {
  if (!human_mobility_enabled_any(parameters)) {
    return(NULL)
  }

  n_nodes <- length(parameters)
  weighted_active <- human_exposure_lag_weighted_active(parameters)
  exposure_defaults <- human_exposure_lag_node_species_defaults(
    parameters,
    variables,
    solvers,
    lagged_values = lagged_eir
  )
  weighted_defaults <- if (isTRUE(weighted_active)) {
    human_exposure_lag_node_species_defaults(
      parameters,
      variables,
      solvers,
      weighted = TRUE,
      lagged_values = lagged_transmission_eir
    )
  } else {
    NULL
  }

  context <- new.env(parent = emptyenv())
  context$n_nodes <- n_nodes
  context$parameters <- parameters
  context$variables <- variables
  context$weighted_active <- weighted_active
  context$reported_timestep <- NA_real_
  context$reported <- rep(FALSE, n_nodes)
  context$node_exposure <- vector("list", n_nodes)
  context$node_weighted_exposure <- if (isTRUE(weighted_active)) vector("list", n_nodes) else NULL

  default_context <- new.env(parent = emptyenv())
  default_context$n_nodes <- n_nodes
  default_context$parameters <- parameters
  default_context$variables <- variables
  default_context$weighted_active <- weighted_active
  default_context$reported_timestep <- NA_real_
  allocated_defaults <- human_exposure_lag_allocate_exposure(
    default_context,
    timestep = 0,
    node_exposure = exposure_defaults,
    node_weighted_exposure = weighted_defaults
  )

  context$buffers <- lapply(
    seq_along(parameters),
    function(i) {
      default_exposure <- allocated_defaults$exposure[[i]]
      default_weighted <- if (isTRUE(weighted_active)) {
        allocated_defaults$weighted_exposure[[i]]
      } else {
        NULL
      }
      HumanExposureLagBuffer$new(
        max_lag = parameters[[i]]$de + 1,
        default_exposure = default_exposure,
        default_weighted_exposure = default_weighted,
        default_species_exposure = allocated_defaults$species_exposure[[i]],
        default_species_weighted_exposure = if (isTRUE(weighted_active)) {
          allocated_defaults$species_weighted_exposure[[i]]
        } else {
          NULL
        }
      )
    }
  )
  context
}

human_exposure_lag_record_node <- function(context, node_index, timestep, exposure, weighted_exposure = NULL) {
  if (is.null(context)) {
    return(invisible(NULL))
  }

  if (is.na(context$reported_timestep) || context$reported_timestep != timestep) {
    context$reported_timestep <- timestep
    context$reported[] <- FALSE
    context$node_exposure[] <- vector("list", context$n_nodes)
    if (isTRUE(context$weighted_active)) {
      context$node_weighted_exposure[] <- vector("list", context$n_nodes)
    }
  }

  expected_species <- length(context$parameters[[node_index]]$species)
  context$reported[[node_index]] <- TRUE
  context$node_exposure[[node_index]] <- human_exposure_lag_validate_node_exposure(
    exposure,
    expected_species,
    "exposure"
  )
  if (isTRUE(context$weighted_active)) {
    context$node_weighted_exposure[[node_index]] <- human_exposure_lag_validate_node_exposure(
      weighted_exposure,
      expected_species,
      "weighted_exposure"
    )
  }

  if (!all(context$reported)) {
    return(invisible(NULL))
  }

  allocated <- human_exposure_lag_allocate_exposure(
    context,
    timestep,
    context$node_exposure,
    context$node_weighted_exposure
  )

  for (home_node in seq_len(context$n_nodes)) {
    context$buffers[[home_node]]$save(
      timestep,
      allocated$exposure[[home_node]],
      if (isTRUE(context$weighted_active)) allocated$weighted_exposure[[home_node]] else NULL,
      species_exposure = allocated$species_exposure[[home_node]],
      species_weighted_exposure = if (isTRUE(context$weighted_active)) {
        allocated$species_weighted_exposure[[home_node]]
      } else {
        NULL
      }
    )
  }

  invisible(NULL)
}

human_exposure_lag_validate_infection_vector <- function(values, expected_length, label) {
  values <- as.numeric(values)
  if (length(values) != expected_length ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values < 0)) {
    stop(sprintf("%s must be a nonnegative finite numeric vector matching the human population.", label), call. = FALSE)
  }
  values
}

human_exposure_lag_get_infection_input <- function(context, node_index, timestep, parameters) {
  if (is.null(context)) {
    stop("Explicit human mobility requires a human exposure lag context.", call. = FALSE)
  }

  node_index <- as.integer(node_index)
  if (length(node_index) != 1L ||
      is.na(node_index) ||
      node_index < 1L ||
      node_index > length(context$buffers)) {
    stop("node_index must select one human exposure lag buffer.", call. = FALSE)
  }

  lookup_timestep <- timestep - parameters$de
  exposure <- context$buffers[[node_index]]$get(lookup_timestep)
  exposure <- human_exposure_lag_validate_infection_vector(
    exposure,
    parameters$human_population,
    "Delayed human bite-sampling exposure"
  )
  species_exposure <- context$buffers[[node_index]]$get(lookup_timestep, by_species = TRUE)
  if (is.null(species_exposure)) {
    species_exposure <- matrix(exposure, nrow = 1L)
  }
  if (nrow(species_exposure) != length(parameters$species) ||
      ncol(species_exposure) != parameters$human_population ||
      anyNA(species_exposure) ||
      any(!is.finite(species_exposure)) ||
      any(species_exposure < 0)) {
    stop("Delayed human bite-sampling exposure by species must be a nonnegative finite matrix.", call. = FALSE)
  }

  weighted_exposure <- NULL
  species_weighted_exposure <- NULL
  transmission_multiplier <- 1
  if (isTRUE(context$weighted_active)) {
    weighted_exposure <- context$buffers[[node_index]]$get(lookup_timestep, weighted = TRUE)
    weighted_exposure <- human_exposure_lag_validate_infection_vector(
      weighted_exposure,
      parameters$human_population,
      "Delayed weighted human bite-sampling exposure"
    )
    species_weighted_exposure <- context$buffers[[node_index]]$get(
      lookup_timestep,
      weighted = TRUE,
      by_species = TRUE
    )

    transmission_multiplier <- rep.int(1, length(exposure))
    positive_exposure <- exposure > 0
    transmission_multiplier[positive_exposure] <- weighted_exposure[positive_exposure] /
      exposure[positive_exposure]
  }

  list(
    infection_exposure = exposure,
    bite_rates_by_species = species_exposure,
    weighted_exposure = weighted_exposure,
    weighted_bite_rates_by_species = species_weighted_exposure,
    transmission_multiplier = transmission_multiplier
  )
}

human_exposure_lag_clear <- function(context, node_index, target) {
  if (is.null(context)) {
    return(invisible(NULL))
  }
  context$buffers[[node_index]]$clear(target)
  invisible(NULL)
}

human_exposure_lag_save_state <- function(context) {
  if (is.null(context)) {
    return(NULL)
  }
  lapply(context$buffers, function(buffer) buffer$save_state())
}

human_exposure_lag_restore_state <- function(context, timestep, state) {
  if (is.null(context) || is.null(state)) {
    return(invisible(NULL))
  }
  for (i in seq_along(context$buffers)) {
    if (length(state) >= i) {
      context$buffers[[i]]$restore_state(timestep, state[[i]])
    }
  }
  invisible(NULL)
}
