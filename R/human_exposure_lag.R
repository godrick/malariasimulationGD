HumanExposureLagBuffer <- R6::R6Class(
  "HumanExposureLagBuffer",
  private = list(
    max_rows = NULL,
    default_exposure = NULL,
    default_weighted_exposure = NULL,
    timesteps = numeric(0),
    exposure = NULL,
    weighted_exposure = NULL,

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
    }
  ),
  public = list(
    initialize = function(max_lag, default_exposure, default_weighted_exposure = NULL) {
      max_rows <- as.integer(ceiling(max_lag) + 2L)
      if (max_rows < 2L) {
        max_rows <- 2L
      }
      default_exposure <- as.numeric(default_exposure)
      if (length(default_exposure) < 1L || anyNA(default_exposure) || any(!is.finite(default_exposure))) {
        stop("default_exposure must be a finite numeric vector.", call. = FALSE)
      }

      private$max_rows <- max_rows
      private$default_exposure <- default_exposure
      private$exposure <- matrix(numeric(0), nrow = 0L, ncol = length(default_exposure))

      if (!is.null(default_weighted_exposure)) {
        default_weighted_exposure <- as.numeric(default_weighted_exposure)
        if (length(default_weighted_exposure) != length(default_exposure) ||
            anyNA(default_weighted_exposure) || any(!is.finite(default_weighted_exposure))) {
          stop("default_weighted_exposure must be a finite numeric vector matching default_exposure.", call. = FALSE)
        }
        private$default_weighted_exposure <- default_weighted_exposure
        private$weighted_exposure <- matrix(numeric(0), nrow = 0L, ncol = length(default_exposure))
      }
    },

    save = function(timestep, exposure, weighted_exposure = NULL) {
      exposure <- private$check_values(exposure, "exposure")
      if (!is.null(private$weighted_exposure)) {
        weighted_exposure <- private$check_values(weighted_exposure, "weighted_exposure")
      }

      existing <- which(private$timesteps == timestep)
      if (length(existing) > 0L) {
        row <- existing[[1L]]
        private$exposure[row, ] <- exposure
        if (!is.null(private$weighted_exposure)) {
          private$weighted_exposure[row, ] <- weighted_exposure
        }
        return(invisible(NULL))
      }

      private$timesteps <- c(private$timesteps, as.numeric(timestep))
      private$exposure <- rbind(private$exposure, exposure)
      if (!is.null(private$weighted_exposure)) {
        private$weighted_exposure <- rbind(private$weighted_exposure, weighted_exposure)
      }

      if (length(private$timesteps) > private$max_rows) {
        keep <- seq.int(length(private$timesteps) - private$max_rows + 1L, length(private$timesteps))
        private$timesteps <- private$timesteps[keep]
        private$exposure <- private$exposure[keep, , drop = FALSE]
        if (!is.null(private$weighted_exposure)) {
          private$weighted_exposure <- private$weighted_exposure[keep, , drop = FALSE]
        }
      }

      invisible(NULL)
    },

    get = function(timestep, weighted = FALSE) {
      history <- private$select_history(weighted)
      default <- private$select_default(weighted)
      if (is.null(history) || is.null(default)) {
        return(NULL)
      }
      if (length(private$timesteps) == 0L || timestep < private$timesteps[[1L]]) {
        return(default)
      }

      exact <- which(private$timesteps == timestep)
      if (length(exact) > 0L) {
        return(as.numeric(history[exact[[1L]], ]))
      }

      if (timestep > private$timesteps[[length(private$timesteps)]]) {
        warning("Exposure lag lookup is after the latest saved timestep; returning per-human defaults.", call. = FALSE)
        return(default)
      }

      after <- which(private$timesteps > timestep)[[1L]]
      before <- after - 1L
      weight <- (timestep - private$timesteps[[before]]) /
        (private$timesteps[[after]] - private$timesteps[[before]])

      as.numeric(history[before, ] + weight * (history[after, ] - history[before, ]))
    },

    clear = function(target) {
      index <- if (inherits(target, "Bitset")) target$to_vector() else as.integer(target)
      if (length(index) == 0L) {
        return(invisible(NULL))
      }
      if (nrow(private$exposure) > 0L) {
        private$exposure[, index] <- matrix(
          private$default_exposure[index],
          nrow = nrow(private$exposure),
          ncol = length(index),
          byrow = TRUE
        )
      }
      if (!is.null(private$weighted_exposure)) {
        if (nrow(private$weighted_exposure) > 0L) {
          private$weighted_exposure[, index] <- matrix(
            private$default_weighted_exposure[index],
            nrow = nrow(private$weighted_exposure),
            ncol = length(index),
            byrow = TRUE
          )
        }
      }
      invisible(NULL)
    },

    save_state = function() {
      list(
        timesteps = private$timesteps,
        exposure = private$exposure,
        weighted_exposure = private$weighted_exposure,
        default_exposure = private$default_exposure,
        default_weighted_exposure = private$default_weighted_exposure
      )
    },

    restore_state = function(timestep, state) {
      if (is.null(state)) {
        return(invisible(NULL))
      }
      private$timesteps <- as.numeric(state$timesteps)
      private$exposure <- as.matrix(state$exposure)
      private$default_exposure <- as.numeric(state$default_exposure)

      if (!is.null(private$weighted_exposure)) {
        private$weighted_exposure <- as.matrix(state$weighted_exposure)
        private$default_weighted_exposure <- as.numeric(state$default_weighted_exposure)
      }
      invisible(NULL)
    }
  )
)

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
  exposure_defaults <- human_exposure_lag_node_defaults(
    parameters,
    variables,
    solvers,
    lagged_values = lagged_eir
  )
  weighted_defaults <- if (isTRUE(weighted_active)) {
    human_exposure_lag_node_defaults(
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
  context$variables <- variables
  context$weighted_active <- weighted_active
  context$reported_timestep <- NA_real_
  context$reported <- rep(FALSE, n_nodes)
  context$node_exposure <- rep(NA_real_, n_nodes)
  context$node_weighted_exposure <- if (isTRUE(weighted_active)) rep(NA_real_, n_nodes) else NULL
  context$buffers <- lapply(
    seq_along(parameters),
    function(i) {
      default_exposure <- rep(exposure_defaults[[i]], parameters[[i]]$human_population)
      default_weighted <- if (isTRUE(weighted_active)) {
        rep(weighted_defaults[[i]], parameters[[i]]$human_population)
      } else {
        NULL
      }
      HumanExposureLagBuffer$new(
        max_lag = parameters[[i]]$de + 1,
        default_exposure = default_exposure,
        default_weighted_exposure = default_weighted
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
    context$node_exposure[] <- NA_real_
    if (isTRUE(context$weighted_active)) {
      context$node_weighted_exposure[] <- NA_real_
    }
  }

  context$reported[[node_index]] <- TRUE
  context$node_exposure[[node_index]] <- as.numeric(exposure)
  if (isTRUE(context$weighted_active)) {
    context$node_weighted_exposure[[node_index]] <- as.numeric(weighted_exposure)
  }

  if (!all(context$reported)) {
    return(invisible(NULL))
  }

  for (home_node in seq_len(context$n_nodes)) {
    current_node <- context$variables[[home_node]]$current_node$get_values()
    exposure_by_human <- context$node_exposure[current_node]
    weighted_by_human <- if (isTRUE(context$weighted_active)) {
      context$node_weighted_exposure[current_node]
    } else {
      NULL
    }
    context$buffers[[home_node]]$save(timestep, exposure_by_human, weighted_by_human)
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
    "Delayed human infection exposure"
  )

  weighted_exposure <- NULL
  transmission_multiplier <- 1
  if (isTRUE(context$weighted_active)) {
    weighted_exposure <- context$buffers[[node_index]]$get(lookup_timestep, weighted = TRUE)
    weighted_exposure <- human_exposure_lag_validate_infection_vector(
      weighted_exposure,
      parameters$human_population,
      "Delayed weighted human infection exposure"
    )

    transmission_multiplier <- rep.int(1, length(exposure))
    positive_exposure <- exposure > 0
    transmission_multiplier[positive_exposure] <- weighted_exposure[positive_exposure] /
      exposure[positive_exposure]
  }

  list(
    infection_exposure = exposure,
    weighted_exposure = weighted_exposure,
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
