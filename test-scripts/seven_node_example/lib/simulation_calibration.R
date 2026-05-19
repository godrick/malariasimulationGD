# ------------------------------------------------------------------------------
# simulation_calibration.R
# ------------------------------------------------------------------------------
# Simulation-based calibration of init_EIR to a target PfPR2-10 under
# *seasonal* dynamics. Wraps the production helpers
# `msimGD_build_baseline_checkpoint()` (warmup runner) and
# `msimGD_run_truth()` (used here to extract realised post-warmup PfPR over
# the final year).
#
# Why this exists: `calibrate_eir_from_pfpr()` in calibrate_eir.R uses the
# *non-seasonal* analytical equilibrium from `malariaEquilibrium::human_equilibrium()`.
# With seasonality switched on, realised mean PfPR can differ noticeably from
# the analytical value at the same init_EIR. This helper calibrates against the
# simulator configuration actually used by the example with a small
# coarse-grid-then-refine search.
#
# Cost: each iteration runs a short warmup (default 2 years) at the
# candidate EIR -- ~25-30 s on a typical desktop. The search reuses the
# analytical-equilibrium EIR as the anchor for the initial grid, so 4
# coarse points plus 1-2 refinements typically suffice.
#
# Returns a list with `init_EIR`, `realised_pfpr`, the search trajectory,
# and the analytical-equilibrium starting EIR.
# ------------------------------------------------------------------------------

#' Run one short warmup at a candidate init_EIR and return the realised
#' mean PfPR over a measurement window.
#'
#' @param init_EIR candidate annual EIR
#' @param setup,cube,NH,mu,p_move,theta,parameter_modifier passed straight to
#'   msimGD_build_baseline_checkpoint(); same semantics as in 01_.
#' @param warmup_days short-warmup length (default 730, two years)
#' @param measurement_window_days how many trailing days from warmup_days to
#'   average PfPR over (default 365 -- the final year of warmup)
#' @param age_band_days integer length-2 vector (min, max) in days; passes to
#'   msimGD's prevalence rendering as prevalence_rendering_min_age and
#'   prevalence_rendering_max_age. Default = c(730, 3650) = PfPR2-10.
#' @return list(init_EIR, realised_pfpr, ...)
.scal_one_warmup_pfpr <- function(init_EIR,
                                  setup, cube, NH, mu, p_move, theta,
                                  parameter_modifier,
                                  warmup_days = 730L,
                                  measurement_window_days = 365L,
                                  age_band_days = c(730L, 3650L),
                                  seed = 20260514L) {
  if (!exists("msimGD_run_truth", mode = "function")) {
    stop("Source lib/msimGD_truth_generation.R before calling this helper.",
         call. = FALSE)
  }

  res <- msimGD_run_truth(
    setup        = setup,
    cube         = cube,
    NF           = NULL,
    NH           = NH,
    tmax         = warmup_days,
    mu           = mu,
    p_move       = p_move,
    release      = NULL,
    theta        = theta,
    init_EIR     = init_EIR,
    prevalence_rendering_min_age = age_band_days[[1]],
    prevalence_rendering_max_age = age_band_days[[2]],
    warmup_days  = 0L,
    parameter_modifier = parameter_modifier,
    seed = as.integer(seed)
  )
  # msimGD_run_truth returns either a list-of-per-node-data.frames directly
  # or a wrapped list. Match the same defensive logic used in 02_.
  data_list <- if (is.list(res) && all(vapply(res, is.data.frame, logical(1)))) res
               else if (!is.null(res$data)) res$data
               else if (!is.null(res$sim$data)) res$sim$data
               else stop("Cannot locate per-node data.frames.", call. = FALSE)

  n_col <- sprintf("n_detect_lm_%d_%d", age_band_days[[1]], age_band_days[[2]])
  d_col <- sprintf("n_age_%d_%d",      age_band_days[[1]], age_band_days[[2]])
  for (col in c(n_col, d_col, "timestep")) {
    for (d in data_list) {
      if (!(col %in% names(d))) {
        stop(sprintf("Expected column '%s' missing from per-node output. ",
                     col),
             "Did you pass prevalence_rendering_min_age / _max_age correctly?",
             call. = FALSE)
      }
    }
  }
  # Aggregate numerator and denominator across nodes per day, then average
  # over the trailing measurement window.
  combined <- do.call(rbind, data_list)
  agg <- aggregate(combined[, c(n_col, d_col)],
                   by = list(timestep = combined$timestep), FUN = sum)
  agg$prev <- agg[[n_col]] / pmax(agg[[d_col]], 1)
  meas_start <- max(1L, warmup_days - measurement_window_days + 1L)
  meas <- agg[agg$timestep >= meas_start & agg$timestep <= warmup_days, ]
  list(
    init_EIR     = init_EIR,
    realised_pfpr = mean(meas$prev),
    meas_min     = min(meas$prev),
    meas_max     = max(meas$prev),
    n_days       = nrow(meas)
  )
}


#' Simulation-based init_EIR calibration to a target *seasonal-mean* PfPR2-10.
#'
#' Uses a small log-spaced coarse grid + log-linear interpolation. Each grid
#' point runs one short warmup via msimGD_run_truth().
#'
#' @param target_pfpr scalar in (0, 1)
#' @param ... see .scal_one_warmup_pfpr() for the per-warmup args
#' @param starting_init_EIR optional anchor (e.g. analytical-equilibrium result);
#'   the coarse grid spans approximately [0.25 * anchor, 1.20 * anchor]. If
#'   NULL the grid spans [0.5, 30] absolute EIR.
#' @param tolerance_rel desired relative tolerance |achieved-target|/target
#' @param max_refine_iter cap on refinement iterations after the coarse grid
#' @return list with init_EIR, realised_pfpr, full search log, analytical_starting_EIR
simulation_calibrate_init_eir <- function(
    target_pfpr,
    setup, cube, NH, mu, p_move, theta, parameter_modifier,
    warmup_days = 730L,
    measurement_window_days = 365L,
    age_band_days = c(730L, 3650L),
    starting_init_EIR = NULL,
    tolerance_rel = 0.05,
    max_refine_iter = 3L,
    seed = 20260514L,
    verbose = TRUE) {

  stopifnot(is.numeric(target_pfpr), length(target_pfpr) == 1L,
            target_pfpr > 0, target_pfpr < 1)

  # Coarse log-spaced grid spanning ~[0.25, 1.2] of the analytical EIR (or
  # [0.5, 30] absolute if no anchor was provided).
  if (is.null(starting_init_EIR)) {
    grid <- exp(seq(log(0.5), log(30), length.out = 4L))
  } else {
    # Wider downward range than upward because the analytical anchor can
    # overstate the EIR needed under the seasonal simulator used here.
    grid <- exp(seq(log(starting_init_EIR * 0.10),
                    log(starting_init_EIR * 1.20),
                    length.out = 4L))
  }

  log_rows <- list()
  one <- function(eir, kind) {
    r <- .scal_one_warmup_pfpr(
      init_EIR = eir,
      setup = setup, cube = cube, NH = NH, mu = mu, p_move = p_move,
      theta = theta, parameter_modifier = parameter_modifier,
      warmup_days = warmup_days,
      measurement_window_days = measurement_window_days,
      age_band_days = age_band_days,
      seed = seed
    )
    if (verbose) {
      cat(sprintf("  [%s] init_EIR=%7.4f  realised_PfPR=%.4f  (target %.4f)\n",
                  kind, eir, r$realised_pfpr, target_pfpr))
    }
    log_rows[[length(log_rows) + 1L]] <<-
      data.frame(kind = kind, init_EIR = eir,
                 realised_pfpr = r$realised_pfpr,
                 stringsAsFactors = FALSE)
    r
  }

  if (verbose) cat("Coarse grid:\n")
  coarse <- lapply(grid, function(e) one(e, "coarse"))

  # Sort by realised PfPR so we can interpolate monotonically
  pfprs <- vapply(coarse, function(r) r$realised_pfpr, numeric(1))
  eirs  <- vapply(coarse, function(r) r$init_EIR,    numeric(1))
  ord <- order(eirs)
  eirs <- eirs[ord]; pfprs <- pfprs[ord]

  if (target_pfpr < min(pfprs)) {
    stop(sprintf(
      "Target PfPR %.4f is below the lowest grid PfPR %.4f (at EIR=%.3f). Lower the grid.",
      target_pfpr, min(pfprs), eirs[which.min(pfprs)]), call. = FALSE)
  }
  if (target_pfpr > max(pfprs)) {
    stop(sprintf(
      "Target PfPR %.4f is above the highest grid PfPR %.4f (at EIR=%.3f). Raise the grid.",
      target_pfpr, max(pfprs), eirs[which.max(pfprs)]), call. = FALSE)
  }

  # Log-linear interpolation on log(EIR) vs PfPR
  log_eirs <- log(eirs)
  predict_eir_at <- function(p) {
    exp(stats::approx(x = pfprs, y = log_eirs, xout = p, rule = 2)$y)
  }

  init_EIR_hat <- predict_eir_at(target_pfpr)

  # Refinement: add the interpolated point, re-fit, repeat until converged.
  for (it in seq_len(max_refine_iter)) {
    refine <- one(init_EIR_hat, sprintf("refine_%d", it))
    log_eirs <- c(log_eirs, log(init_EIR_hat))
    pfprs    <- c(pfprs,    refine$realised_pfpr)
    if (abs(refine$realised_pfpr - target_pfpr) / target_pfpr < tolerance_rel) {
      if (verbose) cat(sprintf("Converged within %.1f%% relative on iter %d.\n",
                               100 * tolerance_rel, it))
      break
    }
    ord <- order(log_eirs)
    log_eirs <- log_eirs[ord]; pfprs <- pfprs[ord]
    # De-duplicate identical PfPRs for approx()
    keep <- !duplicated(pfprs)
    init_EIR_hat <- predict_eir_at(target_pfpr)
  }

  list(
    init_EIR        = init_EIR_hat,
    realised_pfpr   = refine$realised_pfpr,
    target_pfpr     = target_pfpr,
    starting_init_EIR = starting_init_EIR,
    search_log      = do.call(rbind, log_rows)
  )
}
