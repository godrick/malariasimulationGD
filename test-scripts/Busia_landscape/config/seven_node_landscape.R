# ------------------------------------------------------------------------------
# seven_node_landscape.R
# ------------------------------------------------------------------------------
# Specification of the 7-village landscape: village positions in a 14x14 km
# window, drawn by inhibitory rejection sampling with a 1.5 km minimum
# pairwise distance, plus per-village human-population NH drawn from a clipped
# normal so there is village-to-village heterogeneity. Reproducible via a
# fixed seed.
#
# Exposes the function `build_seven_node_landscape()` returning a list:
#   nodes : data.frame with columns (node, x, y, NH_per_node)
#   D     : 7x7 distance matrix (km, symmetric, zero diagonal)
#   n_nodes : 7
#
# The example uses a tiny landscape (14 km side) on purpose. Distances are
# well below the package's default flight-range assumptions, so the movement
# kernel has a meaningful between-village mixing range.
# ------------------------------------------------------------------------------

#' Draw 7-node positions with a min-distance constraint (inhibitory sampling).
#'
#' @param side side length of the square window (km)
#' @param min_dist minimum pairwise distance (km)
#' @param n_nodes number of nodes
#' @param max_tries safety cap on rejection-sampling attempts
#' @param seed RNG seed
#' @return data.frame with columns (node, x, y)
.inhibitory_positions <- function(side = 14, min_dist = 1.5, n_nodes = 7L,
                                  max_tries = 10000L, seed = 20260514L) {
  set.seed(as.integer(seed))
  half <- side / 2
  pts <- matrix(NA_real_, n_nodes, 2)
  filled <- 0L
  tries <- 0L
  while (filled < n_nodes && tries < max_tries) {
    tries <- tries + 1L
    cand <- c(stats::runif(1, -half, half), stats::runif(1, -half, half))
    ok <- TRUE
    if (filled > 0L) {
      d <- sqrt(rowSums((pts[seq_len(filled), , drop = FALSE] -
                          matrix(cand, filled, 2, byrow = TRUE))^2))
      if (any(d < min_dist)) ok <- FALSE
    }
    if (ok) {
      filled <- filled + 1L
      pts[filled, ] <- cand
    }
  }
  if (filled < n_nodes) {
    stop(sprintf(
      "Could not place %d points with min_dist=%.2f km in a %.1f km window after %d tries; try relaxing min_dist or shrinking n_nodes.",
      n_nodes, min_dist, side, max_tries
    ), call. = FALSE)
  }
  data.frame(node = seq_len(n_nodes), x = pts[, 1], y = pts[, 2])
}


#' Draw per-village NH from a clipped normal.
.draw_NH_per_node <- function(n_nodes = 7L, mu_NH = 1500, sd_NH = 500,
                              lower = 800, upper = 2500, seed = 20260515L) {
  set.seed(as.integer(seed))
  NH <- stats::rnorm(n_nodes, mean = mu_NH, sd = sd_NH)
  NH <- pmin(pmax(round(NH), lower), upper)
  as.integer(NH)
}


#' Build the complete 7-node landscape (positions + NH + distance matrix).
build_seven_node_landscape <- function(side = 14, min_dist = 1.5,
                                       n_nodes = 7L,
                                       NH_mu = 1500, NH_sd = 500,
                                       NH_lower = 800, NH_upper = 2500,
                                       seed_positions = 20260514L,
                                       seed_NH = 20260515L) {
  pos <- .inhibitory_positions(side = side, min_dist = min_dist,
                               n_nodes = n_nodes, seed = seed_positions)
  NH_per_node <- .draw_NH_per_node(n_nodes = n_nodes, mu_NH = NH_mu,
                                   sd_NH = NH_sd, lower = NH_lower,
                                   upper = NH_upper, seed = seed_NH)
  nodes <- data.frame(
    node = pos$node,
    x = pos$x,
    y = pos$y,
    NH_per_node = NH_per_node
  )
  D <- as.matrix(stats::dist(nodes[, c("x", "y")]))
  list(nodes = nodes, D = D, n_nodes = nrow(nodes))
}
