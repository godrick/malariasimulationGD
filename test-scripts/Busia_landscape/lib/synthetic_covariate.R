# ------------------------------------------------------------------------------
# synthetic_covariate.R
# ------------------------------------------------------------------------------
# Draw a per-village synthetic covariate z (one number per village) that is
# weakly correlated with an underlying Gaussian-process "field" on the
# landscape. We use this in the example to demonstrate that the contact-rate
# multiplier can vary across villages in a spatially-structured way.
#
# The construction:
#   1. Build a kernel matrix K_ij = exp(-D_ij / theta_km).
#   2. Draw S ~ MVN(0, tau^2 * K). S is a single GP draw across villages.
#   3. Draw eps ~ N(0, 1) per village independently.
#   4. Set z_v = rho * S_v / sd(S) + sqrt(1 - rho^2) * eps_v, so z is
#      approximately standard normal across villages with correlation rho to S.
#
# Args:
#   D : n x n distance matrix in km
#   rho : target correlation between z and S; 0 = independent, 1 = identical
#   tau : GP marginal SD (for the underlying field)
#   theta_km : GP length scale in km
#   seed : RNG seed (use a fresh integer per call to get reproducible draws)
# ------------------------------------------------------------------------------

draw_synthetic_covariate <- function(D, rho = 0.6, tau = 1.0, theta_km = 5.0,
                                     seed = NULL) {
  D <- as.matrix(D)
  storage.mode(D) <- "double"
  n <- nrow(D)
  stopifnot(n == ncol(D), n >= 2)
  stopifnot(is.numeric(rho), length(rho) == 1L, rho >= 0, rho <= 1)
  stopifnot(is.numeric(theta_km), length(theta_km) == 1L, theta_km > 0)
  stopifnot(is.numeric(tau), length(tau) == 1L, tau > 0)

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # Build the exponential kernel with a small ridge for numerical safety.
  K <- tau^2 * exp(-D / theta_km)
  diag(K) <- diag(K) + 1e-8

  # Cholesky factorisation, then S = L %*% z
  L <- chol(K)
  z_S <- stats::rnorm(n)
  S <- as.numeric(crossprod(L, z_S))   # length n

  # Standardise S to unit SD before mixing
  S_std <- (S - mean(S)) / stats::sd(S)

  # Independent noise component
  eps <- stats::rnorm(n)

  # Mix to target correlation
  z <- rho * S_std + sqrt(1 - rho^2) * eps

  list(
    z = as.numeric(z),
    S = S,                    # underlying GP draw (for reference)
    rho = rho,
    tau = tau,
    theta_km = theta_km
  )
}
