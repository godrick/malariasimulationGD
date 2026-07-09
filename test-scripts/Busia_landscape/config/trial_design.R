# ------------------------------------------------------------------------------
# trial_design.R
# ------------------------------------------------------------------------------
# Trial-design knobs for the 7-node example: which villages receive the
# release, how many drive males per release event, when the release fires,
# and the trial horizon. `horizon_day` is measured from release; the absolute
# readout day in a release simulation is `release_day + horizon_day`.
# ------------------------------------------------------------------------------

seven_node_trial_design <- function() {
  list(
    # Three release villages selected by spatial spread from the seven.
    # `pick_release_nodes_by_spread()` below makes this deterministic given
    # the landscape's `nodes` table.
    n_release_nodes = 3L,

    # Release size: number of HH males released per release event per village.
    # In the production audit-cell design this is sized to overwhelm natural
    # mortality so the drive establishes reliably:
    #   release_size = males_per_person * reference_population (26 * 770 = 20020).
    # We use that same value so the example matches what the production
    # 7-node and 20-node audit cells actually release.
    release_size = 20020L,
    release_day  = 90L,     # day of the release run (post-warmup for careful)
    horizon_day  = 365L     # days after release
  )
}


#' Pick `n_release` villages with maximum spatial spread.
#'
#' Greedy maximin: start from the village closest to the area centroid, then
#' at each step add the village whose minimum distance to the already-chosen
#' set is largest. Reproducible (no RNG).
#'
#' @param nodes data.frame with columns (node, x, y)
#' @param n_release number of release villages to pick
#' @return integer vector of node indices, length n_release
pick_release_nodes_by_spread <- function(nodes, n_release = 3L) {
  stopifnot(is.data.frame(nodes), all(c("node", "x", "y") %in% names(nodes)))
  n_release <- as.integer(n_release)
  stopifnot(n_release >= 1L, n_release <= nrow(nodes))

  centroid <- c(mean(nodes$x), mean(nodes$y))
  d_to_centroid <- sqrt((nodes$x - centroid[1])^2 + (nodes$y - centroid[2])^2)
  chosen <- as.integer(nodes$node[which.min(d_to_centroid)])

  D <- as.matrix(stats::dist(nodes[, c("x", "y")]))

  while (length(chosen) < n_release) {
    remaining <- setdiff(nodes$node, chosen)
    chosen_idx <- match(chosen, nodes$node)
    remain_idx <- match(remaining, nodes$node)
    # For each candidate, the min distance to any already-chosen village
    d_min <- vapply(remain_idx,
                    function(i) min(D[i, chosen_idx]),
                    numeric(1))
    chosen <- c(chosen, as.integer(remaining[which.max(d_min)]))
  }
  sort(chosen)
}
