# ------------------------------------------------------------------------------
# covariate.R
# ------------------------------------------------------------------------------
# Generate one synthetic per-village covariate z, weakly correlated to an
# underlying Gaussian-process field on the landscape. The covariate is then
# converted into a per-village contact multiplier `c_v` via a log-linear map:
#
#   log(c_v) = delta * z_v,   so   c_v = exp(delta * z_v)
#
# This example keeps the historical `contact_multiplier` / contact-surface
# field names for compatibility with saved metadata and old scripts. In the
# current package, those scalar node values are folded into uniform
# `human_slot_contact_multiplier` values when human variables are initialized.
# New truth-generation code should write household/slot contact values to
# `human_slot_contact_multiplier` directly instead of adding new uses of
# `contact_multiplier`.
#
# Villages with z above the area mean have higher human blood-meal contact
# intensity than villages with z below. By construction z is approximately
# mean-zero with sd 1 across the 7 villages, so the geometric mean of c_v across
# the population is close to 1.
# ------------------------------------------------------------------------------

seven_node_covariate_settings <- function() {
  list(
    rho = 0.6,         # correlation between z and underlying GP field S
    tau = 1.0,         # GP marginal SD
    theta_km = 5.0,    # GP length scale (km)
    delta = 0.4,       # log-linear coefficient mapping z to contact multiplier
    seed = 9701L
  )
}


#' Build the per-village covariate and contact-multiplier vector.
#'
#' @param D 7x7 distance matrix in km (for the underlying GP)
#' @param settings result of seven_node_covariate_settings()
#' @return list(z, contact_multiplier, settings, underlying_S)
build_seven_node_covariate <- function(D, settings = seven_node_covariate_settings()) {
  # `lib/synthetic_covariate.R` must be sourced before calling this.
  if (!exists("draw_synthetic_covariate", mode = "function")) {
    stop("draw_synthetic_covariate() not in scope. ",
         "Source seven_node_example/lib/synthetic_covariate.R first.",
         call. = FALSE)
  }
  draw <- draw_synthetic_covariate(
    D = D,
    rho = settings$rho,
    tau = settings$tau,
    theta_km = settings$theta_km,
    seed = settings$seed
  )
  contact_multiplier <- exp(settings$delta * draw$z)
  list(
    z = draw$z,
    contact_multiplier = contact_multiplier,
    underlying_S = draw$S,
    settings = settings
  )
}
