human_mobility_enabled <- function(parameters) {
  isTRUE(parameters$human_mobility_enabled)
}

create_human_mobility_context <- function(parameters, variables, timesteps, render_output = TRUE) {
  if (!human_mobility_enabled_any(parameters)) {
    return(NULL)
  }

  n_nodes <- length(parameters)
  store_diagnostics <- any(vapply(
    parameters,
    function(p) isTRUE(p$human_mobility_store_diagnostics),
    logical(1)
  ))

  context <- new.env(parent = emptyenv())
  context$parameters <- parameters
  context$variables <- variables
  context$n_nodes <- n_nodes
  context$move_probs <- parameters[[1L]]$human_move_probs
  context$last_updated_timestep <- NA_integer_
  context$active_od <- matrix(0L, nrow = n_nodes, ncol = n_nodes)
  context$started_od <- matrix(0L, nrow = n_nodes, ncol = n_nodes)
  context$mean_remaining_trip_duration <- rep(0, n_nodes)
  context$store_diagnostics <- isTRUE(store_diagnostics) && isTRUE(render_output)
  if (context$store_diagnostics) {
    context$OD_started_trips <- array(0L, dim = c(timesteps, n_nodes, n_nodes))
    context$OD_active_overnight_stays <- array(0L, dim = c(timesteps, n_nodes, n_nodes))
    context$mean_remaining_trip_duration_history <- matrix(0, nrow = timesteps, ncol = n_nodes)
  }
  context
}

human_mobility_draw_trip_duration <- function(n, parameters) {
  if (n <= 0L) {
    return(integer(0))
  }
  if (identical(parameters$human_trip_duration_type, "fixed")) {
    return(rep.int(as.integer(parameters$human_trip_duration_mean), n))
  }
  as.integer(stats::rgeom(n, prob = 1 / parameters$human_trip_duration_mean) + 1L)
}

human_mobility_update_var <- function(variable, values, index) {
  if (length(index) > 0L) {
    variable$queue_update(values, index)
  }
}

human_mobility_update_current_node <- function(variables) {
  variables$current_node$.update()
}

human_mobility_update_node <- function(context, node_index) {
  vars <- context$variables[[node_index]]
  params <- context$parameters[[node_index]]
  n_nodes <- context$n_nodes

  home_node <- vars$home_node$get_values()
  current_node <- vars$current_node$get_values()
  remaining <- vars$travel_remaining_nights$get_values()
  is_travelling <- vars$is_travelling$get_values()
  exposure_remaining <- remaining
  exposure_is_travelling <- is_travelling
  started_destination <- integer(0)

  cooldown <- which(remaining == -1L)

  active_travel <- which(is_travelling == 1L & remaining > 0L)
  stay_away <- integer(0)
  returning_after_exposure <- integer(0)
  if (length(active_travel) > 0L) {
    stay_away <- active_travel[remaining[active_travel] > 1L]
    returning_after_exposure <- active_travel[remaining[active_travel] == 1L]
  }

  eligible <- which(
    current_node == home_node &
      is_travelling == 0L &
      remaining == 0L
  )
  if (length(eligible) > 0L) {
    draw <- sample.int(
      n_nodes,
      size = length(eligible),
      replace = TRUE,
      prob = context$move_probs[node_index, ]
    )
    travellers <- eligible[draw != node_index]
    if (length(travellers) > 0L) {
      started_destination <- draw[draw != node_index]
      duration <- human_mobility_draw_trip_duration(length(travellers), params)
      vars$current_node$queue_update(started_destination, travellers)
      human_mobility_update_current_node(vars)

      exposure_remaining[travellers] <- duration
      exposure_is_travelling[travellers] <- 1L

      post_exposure_remaining <- ifelse(duration == 1L, -1L, duration - 1L)
      post_exposure_destination <- started_destination
      one_night_travellers <- travellers[duration == 1L]
      if (length(one_night_travellers) > 0L) {
        vars$current_node$queue_update(home_node[one_night_travellers], one_night_travellers)
        post_exposure_destination[duration == 1L] <- home_node[one_night_travellers]
      }
      vars$travel_destination$queue_update(post_exposure_destination, travellers)
      vars$travel_remaining_nights$queue_update(post_exposure_remaining, travellers)
      vars$is_travelling$queue_update(ifelse(duration == 1L, 0L, 1L), travellers)
    }
  }

  current_node <- vars$current_node$get_values()

  if (length(cooldown) > 0L) {
    vars$travel_remaining_nights$queue_update(0L, cooldown)
  }
  if (length(stay_away) > 0L) {
    vars$travel_remaining_nights$queue_update(remaining[stay_away] - 1L, stay_away)
  }
  if (length(returning_after_exposure) > 0L) {
    human_mobility_update_var(vars$current_node, home_node[returning_after_exposure], returning_after_exposure)
    human_mobility_update_var(vars$travel_destination, home_node[returning_after_exposure], returning_after_exposure)
    vars$travel_remaining_nights$queue_update(-1L, returning_after_exposure)
    vars$is_travelling$queue_update(0L, returning_after_exposure)
  }

  context$active_od[node_index, ] <- tabulate(current_node, nbins = n_nodes)
  if (length(started_destination) > 0L) {
    context$started_od[node_index, ] <- tabulate(started_destination, nbins = n_nodes)
  }

  travelling_now <- exposure_is_travelling == 1L & exposure_remaining > 0L
  if (any(travelling_now)) {
    context$mean_remaining_trip_duration[[node_index]] <- mean(exposure_remaining[travelling_now])
  } else {
    context$mean_remaining_trip_duration[[node_index]] <- 0
  }
}

human_mobility_update_context <- function(context, timestep) {
  if (!is.na(context$last_updated_timestep) && context$last_updated_timestep == timestep) {
    return(invisible(NULL))
  }

  context$active_od[,] <- 0L
  context$started_od[,] <- 0L
  context$mean_remaining_trip_duration[] <- 0
  for (node_index in seq_len(context$n_nodes)) {
    human_mobility_update_node(context, node_index)
  }

  if (context$store_diagnostics && timestep <= dim(context$OD_started_trips)[[1L]]) {
    context$OD_started_trips[timestep, , ] <- context$started_od
    context$OD_active_overnight_stays[timestep, , ] <- context$active_od
    context$mean_remaining_trip_duration_history[timestep, ] <- context$mean_remaining_trip_duration
  }

  context$last_updated_timestep <- timestep
  invisible(NULL)
}

human_mobility_render_node <- function(context, node_index, renderer, timestep) {
  active_od <- context$active_od
  started_od <- context$started_od
  humans_present <- colSums(active_od)
  visitors_present <- humans_present - diag(active_od)
  residents_away <- rowSums(active_od) - diag(active_od)
  trips_started <- rowSums(started_od) - diag(started_od)

  renderer$render("humans_present", humans_present[[node_index]], timestep)
  renderer$render("visitors_present", visitors_present[[node_index]], timestep)
  renderer$render("residents_away", residents_away[[node_index]], timestep)
  renderer$render("trips_started", trips_started[[node_index]], timestep)
}

create_human_mobility_process <- function(context, node_index, renderer) {
  if (is.null(context)) {
    return(NULL)
  }
  function(timestep) {
    human_mobility_update_context(context, timestep)
    human_mobility_render_node(context, node_index, renderer, timestep)
  }
}

attach_human_mobility_diagnostics <- function(outputs, context, initial_timesteps = NULL) {
  if (is.null(outputs) || is.null(context) || !isTRUE(context$store_diagnostics)) {
    return(outputs)
  }

  trim <- function(x) {
    if (is.null(initial_timesteps) || initial_timesteps <= 0L) {
      return(x)
    }
    keep <- seq.int(as.integer(initial_timesteps) + 1L, dim(x)[[1L]])
    x[keep, , , drop = FALSE]
  }
  trim_matrix <- function(x) {
    if (is.null(initial_timesteps) || initial_timesteps <= 0L) {
      return(x)
    }
    keep <- seq.int(as.integer(initial_timesteps) + 1L, nrow(x))
    x[keep, , drop = FALSE]
  }

  attr(outputs, "OD_started_trips") <- trim(context$OD_started_trips)
  attr(outputs, "OD_active_overnight_stays") <- trim(context$OD_active_overnight_stays)
  attr(outputs, "mean_remaining_trip_duration") <- trim_matrix(context$mean_remaining_trip_duration_history)
  outputs
}
