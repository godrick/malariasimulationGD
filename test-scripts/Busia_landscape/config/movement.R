# ------------------------------------------------------------------------------
# movement.R
# ------------------------------------------------------------------------------
# Movement parameters and the function that turns them, together with the
# landscape's distance matrix, into the (move_probs, move_rates) pair that
# msimGD's metapop backend consumes.
#
# Single source of truth for the example's movement assumptions. To change
# how mobile the mosquitoes are or how concentrated their destination
# preference is, edit `seven_node_movement_settings()` here.
#
# The quick workflow consumes the returned `move_probs` matrix and the
# per-origin `move_rates` vector directly. The careful workflow consumes
# just `mu` and `p_move` from `seven_node_movement_settings()` and lets the
# production helpers (`msimGD_build_baseline_checkpoint`, `msimGD_run_truth`)
# build the matrix internally.
# ------------------------------------------------------------------------------

seven_node_movement_settings <- function() {
  # For the default 14 km, min_dist=1.5 km landscape the feasibility band
  # `mu_feasible_range()` returns is roughly [3.4, 7.6] km. mu = 5 km sits
  # comfortably in the middle. If you change `seven_node_landscape.R` knobs
  # the feasibility band shifts; pick mu inside the new band, or the
  # workflow scripts will error out before calling msimGD.
  list(
    mu = 5.0,        # target mean realised move distance (km)
    p_move = 0.05    # per-origin probability of a move-out event per timestep
  )
}


#' Build the movement matrix and move-rate vector for the 7-node example.
#'
#' Wraps the production helper `mosquito_movement_from_mu()` from
#' `lib/movement_mu.R` and adapts its output to the field names this
#' example's `quick/` workflow expects (`move_probs`, `move_rates`,
#' `mu_achieved`, `mu_min`, `mu_beta0`, `beta`).
#'
#' Conversion from `p_move` (per-origin probability of moving in one
#' timestep) to the per-origin RATE expected by the kernel:
#'
#'   rate_v = (p_move * muF) / (1 - p_move)
#'
#' where `muF` is the adult-female mortality rate (here 0.132 per day,
#' matching the example's `theta$muF`). This is the same conversion the
#' production audit-cell pipeline uses
#' (`customMGDrive2::calc_move_rate(mu = muF, P = p_move)` inlined).
#'
#' @param D 7x7 distance matrix in km (from build_seven_node_landscape())
#' @param settings result of seven_node_movement_settings()
#' @param muF adult-female daily mortality rate (default 0.132 = An. gambiae)
#' @return list with move_probs (7x7), move_rates (length 7), mu_achieved,
#'   mu_min, mu_beta0, beta, p_move (passed through for diagnostics).
build_seven_node_movement <- function(D, settings = seven_node_movement_settings(),
                                      muF = 0.132) {
  if (!exists("mosquito_movement_from_mu", mode = "function")) {
    stop("mosquito_movement_from_mu() not in scope. ",
         "Source seven_node_example/lib/movement_mu.R first.",
         call. = FALSE)
  }
  n <- nrow(as.matrix(D))
  # Per-origin RATE (not probability). Same formula production uses.
  origin_rate <- (settings$p_move * muF) / (1 - settings$p_move)
  mov <- mosquito_movement_from_mu(
    D              = D,
    mu             = settings$mu,
    move_rates     = rep(origin_rate, n),
    attractiveness = rep(1, n),
    return_diagnostics = TRUE,
    verbose        = FALSE
  )
  # Adapt to slim-format field names the quick workflow uses.
  list(
    move_probs   = mov$mosquito_move_probs,
    move_rates   = mov$mosquito_move_rates,
    mu_achieved  = mov$diagnostics$mu_achieved %||%
                   mov$diagnostics$mu_target %||% settings$mu,
    mu_min       = mov$diagnostics$mu_min,
    mu_beta0     = mov$diagnostics$mu_beta0,
    beta         = mov$beta,
    p_move       = settings$p_move,
    origin_rate  = origin_rate
  )
}

# Provide `%||%` if not already defined by an earlier source.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
