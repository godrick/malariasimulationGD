HumanInfectivityLagBuffer <- R6::R6Class(
  "HumanInfectivityLagBuffer",
  private = list(
    max_rows = NULL,
    n_rows = 0L,
    next_row = 1L,
    default_infectivity = NULL,
    timesteps = numeric(0),
    infectivity = NULL,

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

    ordered_infectivity = function() {
      rows <- private$ordered_rows()
      private$infectivity[rows, , drop = FALSE]
    },

    check_values = function(values, label) {
      values <- as.numeric(values)
      if (length(values) != length(private$default_infectivity) ||
          anyNA(values) ||
          any(!is.finite(values)) ||
          any(values < 0)) {
        stop(sprintf("%s must be a nonnegative finite numeric vector matching the human population.", label), call. = FALSE)
      }
      values
    }
  ),
  public = list(
    initialize = function(max_lag, default_infectivity) {
      max_rows <- as.integer(ceiling(max_lag) + 2L)
      if (max_rows < 2L) {
        max_rows <- 2L
      }

      default_infectivity <- as.numeric(default_infectivity)
      if (length(default_infectivity) < 1L ||
          anyNA(default_infectivity) ||
          any(!is.finite(default_infectivity)) ||
          any(default_infectivity < 0)) {
        stop("default_infectivity must be a nonnegative finite numeric vector.", call. = FALSE)
      }

      private$max_rows <- max_rows
      private$n_rows <- 0L
      private$next_row <- 1L
      private$default_infectivity <- default_infectivity
      private$timesteps <- rep(NA_real_, max_rows)
      private$infectivity <- matrix(0, nrow = max_rows, ncol = length(default_infectivity))
    },

    save = function(timestep, infectivity) {
      infectivity <- private$check_values(infectivity, "infectivity")

      rows <- private$ordered_rows()
      existing <- which(private$timesteps[rows] == timestep)
      if (length(existing) > 0L) {
        private$infectivity[rows[[existing[[1L]]]], ] <- infectivity
        return(invisible(NULL))
      }

      row <- private$next_row
      private$timesteps[[row]] <- as.numeric(timestep)
      private$infectivity[row, ] <- infectivity
      if (private$n_rows < private$max_rows) {
        private$n_rows <- private$n_rows + 1L
      }
      private$next_row <- if (row == private$max_rows) 1L else row + 1L

      invisible(NULL)
    },

    get = function(timestep) {
      timesteps <- private$ordered_timesteps()
      if (length(timesteps) == 0L || timestep < timesteps[[1L]]) {
        return(private$default_infectivity)
      }

      history <- private$ordered_infectivity()
      exact <- which(timesteps == timestep)
      if (length(exact) > 0L) {
        return(as.numeric(history[exact[[1L]], ]))
      }

      if (timestep > timesteps[[length(timesteps)]]) {
        warning("Infectivity lag lookup is after the latest saved timestep; returning per-human defaults.", call. = FALSE)
        return(private$default_infectivity)
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

      private$default_infectivity[index] <- 0
      rows <- private$ordered_rows()
      if (length(rows) > 0L) {
        private$infectivity[rows, index] <- 0
      }
      invisible(NULL)
    },

    save_state = function() {
      list(
        timesteps = private$ordered_timesteps(),
        infectivity = private$ordered_infectivity(),
        default_infectivity = private$default_infectivity
      )
    },

    restore_state = function(timestep, state) {
      if (is.null(state)) {
        return(invisible(NULL))
      }
      timesteps <- as.numeric(state$timesteps)
      infectivity <- as.matrix(state$infectivity)
      n <- min(length(timesteps), private$max_rows)
      if (n > 0L && length(timesteps) > n) {
        keep <- seq.int(length(timesteps) - n + 1L, length(timesteps))
        timesteps <- timesteps[keep]
        infectivity <- infectivity[keep, , drop = FALSE]
      }
      private$timesteps <- rep(NA_real_, private$max_rows)
      private$infectivity <- matrix(0, nrow = private$max_rows, ncol = length(private$default_infectivity))
      private$n_rows <- n
      private$next_row <- if (n == private$max_rows) 1L else n + 1L
      if (n > 0L) {
        private$timesteps[seq_len(n)] <- timesteps
        private$infectivity[seq_len(n), ] <- infectivity
      }
      private$default_infectivity <- as.numeric(state$default_infectivity)
      invisible(NULL)
    }
  )
)

human_infectivity_lag_current_values <- function(timestep, variables, parameters) {
  human_infectivity <- variables$infectivity$get_values()
  if (parameters$tbv) {
    human_infectivity <- account_for_tbv(
      timestep,
      human_infectivity,
      variables,
      parameters
    )
  }
  as.numeric(human_infectivity)
}

create_human_infectivity_lag_context <- function(parameters, variables) {
  if (!human_mobility_enabled_any(parameters)) {
    return(NULL)
  }

  context <- new.env(parent = emptyenv())
  context$n_nodes <- length(parameters)
  context$parameters <- parameters
  context$variables <- variables
  context$reservoir_timestep <- NA_real_
  context$reservoir <- rep(NA_real_, context$n_nodes)
  context$buffers <- lapply(
    seq_along(parameters),
    function(i) {
      HumanInfectivityLagBuffer$new(
        max_lag = parameters[[i]]$delay_gam + 2,
        default_infectivity = variables[[i]]$infectivity$get_values()
      )
    }
  )
  context
}

human_infectivity_lag_record_node <- function(context, node_index, timestep) {
  if (is.null(context)) {
    return(invisible(NULL))
  }

  node_index <- as.integer(node_index)
  infectivity <- human_infectivity_lag_current_values(
    timestep,
    context$variables[[node_index]],
    context$parameters[[node_index]]
  )
  context$buffers[[node_index]]$save(timestep, infectivity)
  context$reservoir_timestep <- NA_real_
  invisible(NULL)
}

create_human_infectivity_lag_process <- function(context, node_index) {
  if (is.null(context)) {
    return(NULL)
  }
  function(timestep) {
    human_infectivity_lag_record_node(context, node_index, timestep)
  }
}

human_infectivity_lag_biting_weights <- function(variables, parameters, timestep) {
  age <- get_age(variables$birth$get_values(), timestep)
  psi <- unique_biting_rate(age, parameters)
  human_biting_weights(
    variables$zeta$get_values(),
    psi,
    human_slot_contact_multiplier_values(variables)
  )
}

human_infectivity_lag_validate_current_node <- function(current_node, n_nodes) {
  current_node <- as.integer(current_node)
  if (anyNA(current_node) || any(current_node < 1L) || any(current_node > n_nodes)) {
    stop("current_node must contain valid node indices for human infectivity aggregation.", call. = FALSE)
  }
  current_node
}

human_infectivity_lag_sum_by_node <- function(values, current_node, n_nodes) {
  out <- rep(0, n_nodes)
  if (length(values) == 0L) {
    return(out)
  }

  grouped <- rowsum(
    matrix(as.numeric(values), ncol = 1L),
    group = current_node,
    reorder = FALSE
  )
  out[as.integer(rownames(grouped))] <- as.numeric(grouped[, 1L])
  out
}

human_infectivity_lag_calculate_reservoir <- function(context, timestep) {
  n_nodes <- context$n_nodes
  numerator <- rep(0, n_nodes)
  denominator <- rep(0, n_nodes)

  for (home_node in seq_len(n_nodes)) {
    params <- context$parameters[[home_node]]
    vars <- context$variables[[home_node]]
    current_node <- human_infectivity_lag_validate_current_node(
      vars$current_node$get_values(),
      n_nodes
    )
    weights <- human_infectivity_lag_biting_weights(vars, params, timestep)
    if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
      stop("Human biting weights must be nonnegative and finite.", call. = FALSE)
    }

    lagged_infectivity <- context$buffers[[home_node]]$get(timestep - params$delay_gam)
    numerator <- numerator + human_infectivity_lag_sum_by_node(
      weights * lagged_infectivity,
      current_node,
      n_nodes
    )
    denominator <- denominator + human_infectivity_lag_sum_by_node(
      weights,
      current_node,
      n_nodes
    )
  }

  reservoir <- rep(0, n_nodes)
  has_humans <- denominator > 0
  reservoir[has_humans] <- numerator[has_humans] / denominator[has_humans]
  reservoir
}

human_infectivity_lag_get_reservoir <- function(context, timestep) {
  if (is.null(context)) {
    stop("Explicit human mobility requires a human infectivity lag context.", call. = FALSE)
  }

  if (is.na(context$reservoir_timestep) || context$reservoir_timestep != timestep) {
    context$reservoir <- human_infectivity_lag_calculate_reservoir(context, timestep)
    context$reservoir_timestep <- timestep
  }
  context$reservoir
}

human_infectivity_lag_get_node_reservoir <- function(context, node_index, timestep) {
  reservoir <- human_infectivity_lag_get_reservoir(context, timestep)
  reservoir[[as.integer(node_index)]]
}

human_infectivity_lag_clear <- function(context, node_index, target) {
  if (is.null(context)) {
    return(invisible(NULL))
  }
  context$buffers[[node_index]]$clear(target)
  context$reservoir_timestep <- NA_real_
  invisible(NULL)
}

human_infectivity_lag_save_state <- function(context) {
  if (is.null(context)) {
    return(NULL)
  }
  lapply(context$buffers, function(buffer) buffer$save_state())
}

human_infectivity_lag_restore_state <- function(context, timestep, state) {
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
