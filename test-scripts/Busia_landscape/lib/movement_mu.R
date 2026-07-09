# ------------------------------------------------------------------------------
# movement_mu.R
# ------------------------------------------------------------------------------
# Mosquito movement kernel parameterized by an interpretable mean move distance `mu`.
#
# Exports:
#   - mosquito_movement_from_mu()
#   - mu_feasible_range()
#   - pick_mu_from_fraction()
#
# Notes:
#   - `allowed` is an adjacency mask: TRUE means i->j moves are allowed to exist at all.
#   - Returned `m_move` should be passed into customMGDrive2::spn_T_*_network(..., m_move = m_move).
#   - Returned `mosquito_move_rates` and `mosquito_move_probs` go into the params list
#     that you pass into customMGDrive2::spn_hazards(...).
# ------------------------------------------------------------------------------

mosquito_movement_from_mu <- function(
    D,
    mu,
    move_rates = NULL,           # length n vector, optional
    move_rate  = 1,              # scalar default if move_rates is NULL
    origin_event_weights = NULL, # length n vector; weights for "random move event" origins (will be normalized)
    attractiveness = NULL,       # length n vector A_j >= 0 (destination weights); default all 1
    allowed = NULL,              # optional n x n logical matrix of allowed moves; diag will be forced FALSE
    method = c("auto", "dense", "block"),
    block_size = 256L,           # for block method
    mu_tol = 1e-10,              # tolerance for mu boundary cases
    root_tol = 1e-10,            # tolerance passed to uniroot
    max_beta = 1e6,              # maximum beta used for bracketing
    max_bracket_iter = 60L,
    verbose = FALSE,
    return_m_move = TRUE,
    return_pairwise_rates = FALSE,
    return_diagnostics = TRUE
) {
  # ----------------------------
  # Helpers (fast rowMax/rowMin if matrixStats is available)
  # ----------------------------
  .rowMaxs <- function(x) {
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      matrixStats::rowMaxs(x)
    } else {
      apply(x, 1L, max)
    }
  }
  .rowMins <- function(x) {
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      matrixStats::rowMins(x)
    } else {
      apply(x, 1L, min)
    }
  }

  # ----------------------------
  # Basic validation + coercions
  # ----------------------------
  if (inherits(D, "dist")) D <- as.matrix(D)
  if (inherits(D, "Matrix")) D <- as.matrix(D)
  if (is.data.frame(D))   D <- as.matrix(D)

  if (!is.matrix(D) || !is.numeric(D)) {
    stop("`D` must be a numeric matrix (or coercible to one).")
  }
  n <- nrow(D)
  if (n != ncol(D)) stop("`D` must be square (n x n).")
  if (n < 1) stop("`D` must have n >= 1.")

  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) || mu < 0) {
    stop("`mu` must be a single finite nonnegative number in the same distance units as `D`.")
  }

  # Ensure double storage (faster math)
  storage.mode(D) <- "double"

  # Negative finite distances are not allowed
  if (any(is.finite(D) & (D < 0), na.rm = TRUE)) {
    stop("`D` contains negative finite distances. Distances must be >= 0.")
  }

  # ----------------------------
  # Attractiveness (destination weights A_j)
  # ----------------------------
  if (is.null(attractiveness)) {
    A <- rep(1, n)
  } else {
    A <- as.numeric(attractiveness)
    if (length(A) != n) stop("`attractiveness` must have length n.")
    if (any(!is.finite(A))) stop("`attractiveness` must be finite (use 0 to exclude a destination).")
    if (any(A < 0)) stop("`attractiveness` must be >= 0.")
    if (all(A == 0)) stop("All entries of `attractiveness` are 0; no destinations can ever be chosen.")
  }
  logA <- log(A)  # log(0) = -Inf is OK

  # ----------------------------
  # Allowed moves mask
  # Default: allowed if finite, not NA, and i != j
  # ----------------------------
  base_allowed <- is.finite(D) & !is.na(D)
  diag(base_allowed) <- FALSE

  if (!is.null(allowed)) {
    if (inherits(allowed, "Matrix")) allowed <- as.matrix(allowed)
    if (!is.matrix(allowed) || any(dim(allowed) != c(n, n))) {
      stop("`allowed` must be an n x n logical matrix.")
    }
    if (!is.logical(allowed)) allowed <- as.logical(allowed)
    diag(allowed) <- FALSE
    base_allowed <- base_allowed & allowed
  }

  # Exclude destinations with A_j == 0 (never choose those columns)
  if (any(A == 0)) {
    base_allowed[, A == 0] <- FALSE
  }

  # This adjacency is what you should pass as `m_move` to customMGDrive2::spn_T_*_network().
  # It is deliberately binary; rates/probabilities are handled separately in hazards.
  m_move <- base_allowed

  row_has_dest <- (rowSums(base_allowed) > 0)

  # ----------------------------
  # Move rates (mosquito_move_rates)
  # ----------------------------
  if (!is.null(move_rates)) {
    move_rates <- as.numeric(move_rates)
    if (length(move_rates) != n) stop("`move_rates` must have length n.")
    if (any(!is.finite(move_rates))) stop("`move_rates` must be finite.")
    if (any(move_rates < 0)) stop("`move_rates` must be >= 0.")

    # If move_rate > 0 but no destinations exist, that's invalid (cannot define conditional probs)
    bad <- which(move_rates > 0 & !row_has_dest)
    if (length(bad) > 0) {
      stop("Some origins have move_rates > 0 but have no allowed destinations. Origins: ",
           paste(bad, collapse = ", "))
    }
  } else {
    if (!is.numeric(move_rate) || length(move_rate) != 1L || !is.finite(move_rate) || move_rate < 0) {
      stop("`move_rate` must be a single finite nonnegative number (used when move_rates is NULL).")
    }
    move_rates <- rep(move_rate, n)
    # If an origin has no destinations, force its move_rate to 0 (so row-sum requirement doesn't apply)
    move_rates[!row_has_dest] <- 0
  }

  # Active origins = those that can generate move events
  active_origins <- which(move_rates > 0 & row_has_dest)

  if (length(active_origins) == 0L) {
    # No moves ever happen
    if (mu > 0) {
      stop("All move_rates are 0 (or no destinations exist). You requested mu > 0, which is incompatible.")
    }
    P0 <- matrix(0, n, n)
    out <- list(
      mosquito_move_rates = move_rates,
      mosquito_move_probs = P0,
      beta = NA_real_
    )
    if (return_m_move) out$m_move <- m_move
    if (return_pairwise_rates) out$pairwise_move_rates <- P0
    if (return_diagnostics) {
      out$diagnostics <- list(
        mu_target = mu,
        mu_achieved = 0,
        mu_min = NA_real_,
        mu_beta0 = NA_real_,
        origin_event_weights = rep(0, n),
        method = "none (no movement)"
      )
    }
    return(out)
  }

  # ----------------------------
  # Origin weights for "random move event"
  # ----------------------------
  if (is.null(origin_event_weights)) {
    w <- rep(0, n)
    w[active_origins] <- 1
  } else {
    w <- as.numeric(origin_event_weights)
    if (length(w) != n) stop("`origin_event_weights` must have length n.")
    if (any(!is.finite(w))) stop("`origin_event_weights` must be finite.")
    if (any(w < 0)) stop("`origin_event_weights` must be >= 0.")
    if (sum(w) <= 0) stop("`origin_event_weights` must have positive sum.")

    # Enforce: no weight on origins that cannot generate move events
    badw <- which(w > 0 & !(move_rates > 0 & row_has_dest))
    if (length(badw) > 0) {
      stop("`origin_event_weights` assigns positive weight to origins with move_rate==0 or no destinations. Origins: ",
           paste(badw, collapse = ", "))
    }
  }
  # Normalize over active origins only (so weights sum to 1 over move events)
  w_sum <- sum(w[active_origins])
  if (w_sum <= 0) stop("After restricting to active origins, origin_event_weights sum to 0.")
  w <- w / w_sum

  idx_event <- active_origins  # where w>0 by construction

  # ----------------------------
  # Precompute safe distance matrices:
  # D_exp: used inside exp(-beta*D); disallowed -> Inf, diag -> Inf
  # D_num: used in numerator sums; disallowed -> 0, diag -> 0
  # ----------------------------
  D_exp <- D
  D_exp[!is.finite(D_exp)] <- Inf
  D_num <- D
  D_num[!is.finite(D_num)] <- 0

  # Apply allowed mask to both
  D_exp[!base_allowed] <- Inf
  D_num[!base_allowed] <- 0
  diag(D_exp) <- Inf
  diag(D_num) <- 0

  # ----------------------------
  # Compute boundary means:
  # mu_min = sum_i w_i * min_j D_ij (over allowed destinations)
  # mu_beta0 = mean at beta=0  (p_ij ∝ A_j among allowed)
  # ----------------------------
  dmin <- .rowMins(D_exp)  # Inf for rows with no destinations
  mu_min <- sum(w[idx_event] * dmin[idx_event])

  # mu(beta=0): p_ij ∝ A_j among allowed
  weights0 <- sweep(base_allowed * 1.0, 2L, A, `*`) # n x n
  denom0 <- rowSums(weights0)
  denom0_safe <- denom0
  denom0_safe[denom0_safe == 0] <- 1
  P_beta0 <- weights0 / denom0_safe
  if (any(denom0 == 0)) P_beta0[denom0 == 0, ] <- 0

  mu_beta0 <- sum(w[idx_event] * rowSums(P_beta0[idx_event, , drop = FALSE] * D_num[idx_event, , drop = FALSE]))

  if (verbose) {
    message(sprintf("Feasible mean range under distance-decay kernel: [mu_min=%.6g, mu_beta0=%.6g]",
                    mu_min, mu_beta0))
  }

  # Feasibility checks for beta >= 0 distance-decay
  if (mu < mu_min - mu_tol) {
    stop(sprintf("Requested mu=%.6g is smaller than the minimum achievable mu_min=%.6g given allowed moves & distances.",
                 mu, mu_min))
  }
  if (mu > mu_beta0 + mu_tol) {
    stop(sprintf(paste0(
      "Requested mu=%.6g is larger than mu(beta=0)=%.6g under the distance-decay family p_ij ∝ A_j exp(-beta D_ij) with beta>=0.\n",
      "This family cannot produce such a large mean distance without allowing beta<0 (distance preference) or changing the model."
    ), mu, mu_beta0))
  }

  # ----------------------------
  # Choose computation method
  # ----------------------------
  method <- match.arg(method)
  if (method == "auto") {
    # Simple heuristic: dense is usually faster up to ~1000-1500
    method <- if (n <= 1200L) "dense" else "block"
  }
  block_size <- as.integer(block_size)
  if (block_size < 1L) stop("`block_size` must be >= 1.")

  # ----------------------------
  # Core: compute mu(beta)
  # ----------------------------
  mu_of_beta_dense <- function(beta) {
    if (beta <= 0) return(mu_beta0)

    L <- -beta * D_exp
    L <- sweep(L, 2L, logA, `+`)
    m <- .rowMaxs(L)
    m[!is.finite(m)] <- 0
    W <- exp(L - m)
    denom <- rowSums(W)
    num <- rowSums(W * D_num)
    mu_i <- num / denom
    sum(w[idx_event] * mu_i[idx_event])
  }

  mu_of_beta_block <- function(beta) {
    if (beta <= 0) return(mu_beta0)

    m <- rep(-Inf, n)

    # Pass 1: row maxima
    for (start in seq.int(1L, n, by = block_size)) {
      end <- min(n, start + block_size - 1L)
      idx <- start:end
      Lb <- -beta * D_exp[, idx, drop = FALSE]
      Lb <- sweep(Lb, 2L, logA[idx], `+`)
      mb <- .rowMaxs(Lb)
      m <- pmax(m, mb)
    }
    m[!is.finite(m)] <- 0

    # Pass 2: denom and numerator
    denom <- numeric(n)
    num   <- numeric(n)

    for (start in seq.int(1L, n, by = block_size)) {
      end <- min(n, start + block_size - 1L)
      idx <- start:end
      Lb <- -beta * D_exp[, idx, drop = FALSE]
      Lb <- sweep(Lb, 2L, logA[idx], `+`)
      Wb <- exp(Lb - m)
      denom <- denom + rowSums(Wb)
      num   <- num   + rowSums(Wb * D_num[, idx, drop = FALSE])
    }

    mu_i <- num / denom
    sum(w[idx_event] * mu_i[idx_event])
  }

  mu_of_beta <- if (method == "dense") mu_of_beta_dense else mu_of_beta_block

  # ----------------------------
  # Solve for beta
  # ----------------------------
  if (abs(mu - mu_beta0) <= mu_tol) {
    beta_hat <- 0
  } else if (abs(mu - mu_min) <= mu_tol) {
    beta_hat <- Inf
  } else {
    beta_hi <- 1
    mu_hi <- mu_of_beta(beta_hi)
    iter <- 0L
    while (mu_hi > mu && beta_hi < max_beta && iter < max_bracket_iter) {
      beta_hi <- beta_hi * 2
      if (beta_hi > max_beta) beta_hi <- max_beta
      mu_hi <- mu_of_beta(beta_hi)
      iter <- iter + 1L
      if (verbose) message(sprintf("Bracketing: beta_hi=%.6g, mu(beta_hi)=%.6g", beta_hi, mu_hi))
      if (beta_hi >= max_beta) break
    }
    if (mu_hi > mu) {
      stop(sprintf(
        paste0("Failed to bracket beta up to max_beta=%.6g.\n",
               "mu(beta) at max_beta is %.6g but target mu is %.6g.\n",
               "Try increasing `max_beta`, rescaling distances (e.g., km not m), or check mu feasibility."),
        max_beta, mu_hi, mu
      ))
    }

    f <- function(b) mu_of_beta(b) - mu
    beta_hat <- uniroot(f, interval = c(0, beta_hi), tol = root_tol)$root
  }

  # ----------------------------
  # Build the final move probability matrix P
  # ----------------------------
  make_P_beta0 <- function() {
    P <- P_beta0
    if (any(move_rates == 0)) P[move_rates == 0, ] <- 0
    diag(P) <- 0
    P
  }

  make_P_infty <- function() {
    P <- matrix(0, n, n)
    for (i in seq_len(n)) {
      if (move_rates[i] <= 0) next
      di <- D_exp[i, ]
      dmin_i <- min(di)
      if (!is.finite(dmin_i)) next
      idx <- which(di == dmin_i)
      if (length(idx) == 0L) next
      wi <- A[idx]
      s  <- sum(wi)
      if (s <= 0) next
      P[i, idx] <- wi / s
    }
    diag(P) <- 0
    P
  }

  make_P_finite_beta_dense <- function(beta) {
    if (beta <= 0) return(make_P_beta0())

    L <- -beta * D_exp
    L <- sweep(L, 2L, logA, `+`)
    m <- .rowMaxs(L)
    m[!is.finite(m)] <- 0
    W <- exp(L - m)
    denom <- rowSums(W)
    denom_safe <- denom
    denom_safe[denom_safe == 0] <- 1
    P <- W / denom_safe
    if (any(denom == 0)) P[denom == 0, ] <- 0
    if (any(move_rates == 0)) P[move_rates == 0, ] <- 0
    diag(P) <- 0

    rs <- rowSums(P)
    fix <- which(move_rates > 0 & rs > 0)
    if (length(fix) > 0L) {
      P[fix, ] <- P[fix, , drop = FALSE] / rs[fix]
    }
    P
  }

  make_P_finite_beta_block <- function(beta) {
    if (beta <= 0) return(make_P_beta0())

    m <- rep(-Inf, n)
    for (start in seq.int(1L, n, by = block_size)) {
      end <- min(n, start + block_size - 1L)
      idx <- start:end
      Lb <- -beta * D_exp[, idx, drop = FALSE]
      Lb <- sweep(Lb, 2L, logA[idx], `+`)
      mb <- .rowMaxs(Lb)
      m <- pmax(m, mb)
    }
    m[!is.finite(m)] <- 0

    denom <- numeric(n)
    for (start in seq.int(1L, n, by = block_size)) {
      end <- min(n, start + block_size - 1L)
      idx <- start:end
      Lb <- -beta * D_exp[, idx, drop = FALSE]
      Lb <- sweep(Lb, 2L, logA[idx], `+`)
      Wb <- exp(Lb - m)
      denom <- denom + rowSums(Wb)
    }
    denom_safe <- denom
    denom_safe[denom_safe == 0] <- 1

    P <- matrix(0, n, n)
    for (start in seq.int(1L, n, by = block_size)) {
      end <- min(n, start + block_size - 1L)
      idx <- start:end
      Lb <- -beta * D_exp[, idx, drop = FALSE]
      Lb <- sweep(Lb, 2L, logA[idx], `+`)
      Wb <- exp(Lb - m)
      Pb <- Wb / denom_safe
      P[, idx] <- Pb
    }

    if (any(denom == 0)) P[denom == 0, ] <- 0
    if (any(move_rates == 0)) P[move_rates == 0, ] <- 0
    diag(P) <- 0

    rs <- rowSums(P)
    fix <- which(move_rates > 0 & rs > 0)
    if (length(fix) > 0L) {
      P[fix, ] <- P[fix, , drop = FALSE] / rs[fix]
    }
    P
  }

  if (is.infinite(beta_hat)) {
    P <- make_P_infty()
  } else {
    P <- if (method == "dense") make_P_finite_beta_dense(beta_hat) else make_P_finite_beta_block(beta_hat)
  }

  # Final sanity checks
  if (any(P < -1e-12, na.rm = TRUE)) stop("Internal error: negative probabilities produced.")
  P[P < 0] <- 0

  rs <- rowSums(P)
  bad_rows <- which(move_rates > 0 & abs(rs - 1) > 1e-8)
  if (length(bad_rows) > 0L) {
    stop("Row-sum check failed for some rows (move_rates>0). Rows: ", paste(bad_rows, collapse = ", "))
  }
  if (any(diag(P) != 0)) stop("Internal error: diag(P) is not zero.")

  mu_achieved <- sum(w[idx_event] * rowSums(P[idx_event, , drop = FALSE] * D_num[idx_event, , drop = FALSE]))

  out <- list(
    mosquito_move_rates = move_rates,
    mosquito_move_probs = P,
    beta = beta_hat
  )
  if (return_m_move) out$m_move <- m_move
  if (return_pairwise_rates) {
    out$pairwise_move_rates <- P * move_rates
  }
  if (return_diagnostics) {
    out$diagnostics <- list(
      mu_target = mu,
      mu_achieved = mu_achieved,
      mu_min = mu_min,
      mu_beta0 = mu_beta0,
      origin_event_weights = w,
      active_origins = idx_event,
      method = method,
      block_size = if (method == "block") block_size else NA_integer_
    )
  }
  out
}

mu_feasible_range <- function(
    D,
    allowed = NULL,
    attractiveness = NULL,
    origin_event_weights = NULL,
    move_rates = NULL,
    move_rate = 1
) {
  if (inherits(D, "dist")) D <- as.matrix(D)
  if (!is.matrix(D) || !is.numeric(D)) stop("D must be a numeric matrix.", call. = FALSE)
  n <- nrow(D)
  if (n != ncol(D)) stop("D must be square.", call. = FALSE)

  # Attractiveness A_j
  if (is.null(attractiveness)) {
    A <- rep(1, n)
  } else {
    A <- as.numeric(attractiveness)
    if (length(A) != n) stop("attractiveness must have length n.", call. = FALSE)
    if (any(!is.finite(A)) || any(A < 0)) stop("attractiveness must be finite and >= 0.", call. = FALSE)
    if (all(A == 0)) stop("All attractiveness are 0; no destinations.", call. = FALSE)
  }

  # Allowed mask base: finite + not NA + not diagonal
  base_allowed <- is.finite(D) & !is.na(D)
  diag(base_allowed) <- FALSE

  if (!is.null(allowed)) {
    if (!is.matrix(allowed) || any(dim(allowed) != c(n, n))) {
      stop("allowed must be n x n logical matrix.", call. = FALSE)
    }
    if (!is.logical(allowed)) allowed <- as.logical(allowed)
    diag(allowed) <- FALSE
    base_allowed <- base_allowed & allowed
  }

  if (any(A == 0)) base_allowed[, A == 0] <- FALSE

  row_has_dest <- rowSums(base_allowed) > 0

  # Move rates -> active origins
  if (!is.null(move_rates)) {
    lam <- as.numeric(move_rates)
    if (length(lam) != n) stop("move_rates must have length n.", call. = FALSE)
    if (any(!is.finite(lam)) || any(lam < 0)) stop("move_rates must be finite and >=0.", call. = FALSE)
    bad <- which(lam > 0 & !row_has_dest)
    if (length(bad)) stop("Some origins have move_rates>0 but no allowed destinations.", call. = FALSE)
  } else {
    if (!is.numeric(move_rate) || length(move_rate) != 1 || !is.finite(move_rate) || move_rate < 0) {
      stop("move_rate must be a single finite >=0 number.", call. = FALSE)
    }
    lam <- rep(move_rate, n)
    lam[!row_has_dest] <- 0
  }

  active <- which(lam > 0 & row_has_dest)
  if (!length(active)) {
    return(list(mu_min = NA_real_, mu_beta0 = NA_real_, active_origins = integer(0),
                origin_event_weights = rep(0, n)))
  }

  # Origin weights
  if (is.null(origin_event_weights)) {
    w <- rep(0, n); w[active] <- 1
  } else {
    w <- as.numeric(origin_event_weights)
    if (length(w) != n) stop("origin_event_weights must have length n.", call. = FALSE)
    if (any(!is.finite(w)) || any(w < 0)) stop("origin_event_weights must be finite and >=0.", call. = FALSE)
    badw <- which(w > 0 & !(lam > 0 & row_has_dest))
    if (length(badw)) stop("origin_event_weights puts mass on inactive origins.", call. = FALSE)
  }
  w <- w / sum(w[active])

  D_exp <- D
  D_exp[!is.finite(D_exp)] <- Inf
  D_exp[!base_allowed] <- Inf
  diag(D_exp) <- Inf

  D_num <- D
  D_num[!is.finite(D_num)] <- 0
  D_num[!base_allowed] <- 0
  diag(D_num) <- 0

  dmin <- apply(D_exp, 1, min)
  mu_min <- sum(w[active] * dmin[active])

  weights0 <- sweep(base_allowed * 1.0, 2, A, `*`)
  denom0 <- rowSums(weights0)
  denom0_safe <- denom0; denom0_safe[denom0_safe == 0] <- 1
  P0 <- weights0 / denom0_safe
  P0[denom0 == 0, ] <- 0

  mu_beta0 <- sum(w[active] * rowSums(P0[active, , drop=FALSE] * D_num[active, , drop=FALSE]))

  list(mu_min = mu_min, mu_beta0 = mu_beta0,
       active_origins = active,
       origin_event_weights = w)
}

pick_mu_from_fraction <- function(mu_min, mu_beta0, fraction = 0.5) {
  if (!is.finite(mu_min) || !is.finite(mu_beta0)) stop("mu_min and mu_beta0 must be finite.", call. = FALSE)
  fraction <- max(0, min(1, fraction))
  mu_min + fraction * (mu_beta0 - mu_min)
}
