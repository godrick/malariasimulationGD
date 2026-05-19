# ------------------------------------------------------------------------------
# baseline_library.R
# ------------------------------------------------------------------------------
# Compact example-local implementation of the production audit-cell
# baseline-library promotion pattern. It keeps the same architecture, defaults,
# and diagnostic signals as the Ochomo helpers in customMGDrive2's
# scripts/msimGD/test-scripts/malaria-sim-epi/ochomo_baseline_calibration_utils.R
# (functions like `ochomo_run_library_grid_search`,
# `stationary_checkpoint_library_select`, the stationarity battery, and
# `ochomo_compute_operational_objective`), expressed as example-local code
# instead of vendoring the full Ochomo helper stack. Tolerance choices are
# documented below and are deliberately example-local where needed for a
# single-replicate share example.
#
# Functions exposed:
#   build_baseline_checkpoint_library(): run burnin + take N snapshots
#   score_stationarity_per_snapshot():   per-snapshot Y/Y stationarity metrics
#   promote_best_snapshot():              pick best-scoring snapshot
#
# Production defaults baked in (overridable):
#   burnin_timesteps          2190     (6 years; matches audit-cell baseline)
#   n_snapshots                  6     (matches checkpoint_n_snapshots = 6)
#   snapshot_spacing           365     (matches checkpoint_snapshot_spacing = 365)
#   stationarity scope        "aggregate_only"
#   tol_rel_change_pfpr       0.02     (looser than production's 0.01 because
#                                       the share example is single-rep; the
#                                       audit cells run M~240 and average, so
#                                       can afford 0.01. 0.02 here keeps the
#                                       pass rate sensible at M = 1.)
#   tol_rel_change_M          0.05
#   tol_cycle_rmse_rel_pfpr   0.10     (matches production)
#   tol_cycle_rmse_rel_M      0.10
#
# Signals checked Y/Y on the validation run:
#   - PfPR in the calibration age band (microscopy positivity)
#   - regional total_M_by_node
#
# Each snapshot is one (state, metadata) pair. The library is a list of
# snapshots plus the promotion record.
# ------------------------------------------------------------------------------

# Sanity: needs the vendored production helpers in scope (msimGD_run_truth,
# msimGD_build_baseline_checkpoint).
.bl_needs <- function() {
  for (h in c("msimGD_build_baseline_checkpoint", "msimGD_run_truth")) {
    if (!exists(h, mode = "function")) {
      stop(sprintf(
        "%s() not in scope. Source lib/msimGD_truth_generation.R first.", h
      ), call. = FALSE)
    }
  }
}

# Helper: pull the per-node data.frame list from an msimGD_run_truth result
# (handles either bare list-of-data.frames or wrapped $data).
.bl_get_data_list <- function(res) {
  if (is.list(res) && all(vapply(res, is.data.frame, logical(1)))) res
  else if (!is.null(res$data)) res$data
  else if (!is.null(res$sim$data)) res$sim$data
  else stop("Cannot locate per-node data.frames in msimGD_run_truth output.",
            call. = FALSE)
}

# Helper: regional PfPR (numerator/denominator summed across nodes per day).
.bl_regional_pfpr <- function(data_list, age_band_days) {
  n_col <- sprintf("n_detect_lm_%d_%d", age_band_days[[1]], age_band_days[[2]])
  d_col <- sprintf("n_age_%d_%d",       age_band_days[[1]], age_band_days[[2]])
  combined <- do.call(rbind, data_list)
  agg <- aggregate(combined[, c(n_col, d_col)],
                   by = list(timestep = combined$timestep), FUN = sum)
  data.frame(
    timestep = agg$timestep,
    prev     = agg[[n_col]] / pmax(agg[[d_col]], 1)
  )
}

.bl_regional_total_M <- function(data_list) {
  combined <- do.call(rbind, data_list)
  agg <- aggregate(combined$total_M_gamb,
                   by = list(timestep = combined$timestep), FUN = sum)
  names(agg)[2] <- "total_M"
  agg
}


# ------------------------------------------------------------------------------
# build_baseline_checkpoint_library
# ------------------------------------------------------------------------------
#' Run a burnin then take N snapshots spaced `snapshot_spacing` days apart.
#'
#' Snapshots are captured by chaining `msimGD_run_truth(baseline_checkpoint=)`:
#' the first snapshot is the state at `burnin_timesteps`, snapshot k>1 is the
#' state at `burnin_timesteps + (k-1) * snapshot_spacing`.
#'
#' @param burnin_timesteps integer (days); pre-snapshot warmup
#' @param n_snapshots integer >= 1
#' @param snapshot_spacing integer (days); days between successive snapshots
#' @param setup,cube,NH,mu,p_move,theta,init_EIR,parameter_modifier passed to
#'   msimGD_build_baseline_checkpoint() / msimGD_run_truth()
#' @param seed master seed (snapshot k uses seed + k - 1)
#' @return list(snapshots = list of {state, metadata, snapshot_index,
#'   timesteps_from_t0}, metadata)
build_baseline_checkpoint_library <- function(
    burnin_timesteps = 2190L,
    n_snapshots = 6L,
    snapshot_spacing = 365L,
    setup, cube, NH, mu, p_move, theta, init_EIR,
    parameter_modifier, seed = 20260514L,
    verbose = TRUE) {
  .bl_needs()
  if (!exists("msimGD_build_stationary_checkpoint_library", mode = "function")) {
    stop("msimGD_build_stationary_checkpoint_library() not in scope.",
         call. = FALSE)
  }
  burnin_timesteps <- as.integer(burnin_timesteps)
  n_snapshots      <- as.integer(n_snapshots)
  snapshot_spacing <- as.integer(snapshot_spacing)
  stopifnot(burnin_timesteps >= 0L, n_snapshots >= 1L, snapshot_spacing >= 1L)

  if (verbose) {
    cat(sprintf("[%s] baseline library: burnin %d days, then %d snapshots %d days apart\n",
                format(Sys.time(), "%F %T"),
                burnin_timesteps, n_snapshots, snapshot_spacing))
  }
  t0 <- Sys.time()

  # Delegate to the production helper, which simulates ONE continuous
  # `run_metapop_simulation` chain across all snapshot timesteps with
  # restore_random_state = TRUE between chunks. Much cheaper than calling
  # msimGD_build_baseline_checkpoint multiple times, and bit-identical to
  # the audit-cell pipeline.
  prod_lib <- msimGD_build_stationary_checkpoint_library(
    setup        = setup,
    cube         = cube,
    NF           = NULL,
    NH           = NH,
    mu           = mu,
    p_move       = p_move,
    theta        = theta,
    init_EIR     = init_EIR,
    parameter_modifier = parameter_modifier,
    burnin_timesteps   = burnin_timesteps,
    n_snapshots        = n_snapshots,
    snapshot_spacing   = snapshot_spacing,
    seed         = as.integer(seed),
    library_path = NULL
  )

  # Adapt to the {snapshots, metadata} shape this lib uses internally.
  # Carry through the production checkpoint metadata (especially the
  # `baseline_time_dependent_signature` and `baseline_contact_signature`
  # fields) so downstream `msimGD_run_truth(baseline_checkpoint = ...)`
  # calls accept the snapshot under seasonality.
  snapshots <- lapply(seq_along(prod_lib$checkpoints), function(i) {
    ck <- prod_lib$checkpoints[[i]]
    md <- ck$metadata
    md$snapshot_index    <- as.integer(i)
    md$timesteps_from_t0 <- as.integer(prod_lib$metadata$snapshot_timesteps[[i]])
    md$init_EIR          <- init_EIR
    md$seed              <- as.integer(seed)
    list(state = ck$state, metadata = md)
  })

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (verbose) cat(sprintf("[%s] library built in %.1f s (%.1f min).\n",
                           format(Sys.time(), "%F %T"),
                           elapsed, elapsed / 60))

  list(
    snapshots = snapshots,
    metadata = list(
      burnin_timesteps = burnin_timesteps,
      n_snapshots      = n_snapshots,
      snapshot_spacing = snapshot_spacing,
      init_EIR         = init_EIR,
      master_seed      = as.integer(seed),
      elapsed_seconds  = elapsed,
      prod_metadata    = prod_lib$metadata
    )
  )
}


# ------------------------------------------------------------------------------
# score_stationarity_per_snapshot
# ------------------------------------------------------------------------------
#' For each snapshot, run a 2-year validation, compute Y/Y rel_change and
#' cycle RMSE on PfPR and total_M, score each snapshot.
#'
#' A snapshot's score is a per-signal weighted sum of distances from
#' tolerance — lower is better. Each signal's contribution is max(0,
#' metric / tolerance - 1)^2; sum over signals.
#'
#' @param library result of build_baseline_checkpoint_library()
#' @return list(snapshots_scored = list of {snapshot, metrics, score},
#'              ...)
score_stationarity_per_snapshot <- function(
    library,
    setup, cube, NH, mu, p_move, theta, init_EIR, parameter_modifier,
    age_band_days = c(730L, 3650L),
    validation_years = 2L,
    tol_rel_change_pfpr = 0.02,
    tol_rel_change_M    = 0.05,
    tol_cycle_rmse_rel_pfpr = 0.10,
    tol_cycle_rmse_rel_M    = 0.10,
    seed_offset_per_snapshot = 1000L,
    verbose = TRUE) {
  .bl_needs()
  validation_years <- as.integer(validation_years)
  stopifnot(validation_years >= 2L)
  validation_days <- validation_years * 365L

  metric_one <- function(values, tol_rc, tol_rmse) {
    yN  <- values[(validation_days - 365L + 1L):validation_days]
    yN1 <- values[1L:365L]
    m_yN  <- mean(yN);  m_yN1 <- mean(yN1)
    rel_change      <- abs(m_yN - m_yN1) / max(abs(m_yN1), .Machine$double.eps)
    cycle_rmse_rel  <- sqrt(mean((yN - yN1)^2)) /
      max(abs(m_yN1), .Machine$double.eps)
    list(
      mean_yearN = m_yN, mean_yearN1 = m_yN1,
      rel_change = rel_change,
      cycle_rmse_rel = cycle_rmse_rel,
      pass_rel = rel_change <= tol_rc,
      pass_rmse = cycle_rmse_rel <= tol_rmse
    )
  }
  contribution <- function(m, tol_rc, tol_rmse) {
    max(0, m$rel_change / tol_rc - 1)^2 +
      max(0, m$cycle_rmse_rel / tol_rmse - 1)^2
  }

  scored <- vector("list", length(library$snapshots))
  for (k in seq_along(library$snapshots)) {
    snap <- library$snapshots[[k]]
    if (verbose) cat(sprintf(
      "[%s] scoring snapshot %d/%d (t = %d days)...\n",
      format(Sys.time(), "%F %T"), k, length(library$snapshots),
      snap$metadata$timesteps_from_t0
    ))
    t0 <- Sys.time()
    res <- msimGD_run_truth(
      setup = setup, cube = cube, NF = NULL, NH = NH,
      tmax = validation_days, mu = mu, p_move = p_move, release = NULL,
      theta = theta, init_EIR = snap$metadata$init_EIR,
      prevalence_rendering_min_age = age_band_days[[1]],
      prevalence_rendering_max_age = age_band_days[[2]],
      warmup_days = 0L,
      parameter_modifier = parameter_modifier,
      baseline_checkpoint = list(state = snap$state, metadata = snap$metadata),
      seed = snap$metadata$seed + seed_offset_per_snapshot
    )
    data_list <- .bl_get_data_list(res)
    pfpr_df <- .bl_regional_pfpr(data_list, age_band_days)
    M_df    <- .bl_regional_total_M(data_list)
    m_pfpr <- metric_one(pfpr_df$prev,
                         tol_rel_change_pfpr, tol_cycle_rmse_rel_pfpr)
    m_M    <- metric_one(M_df$total_M,
                         tol_rel_change_M,    tol_cycle_rmse_rel_M)
    score  <- contribution(m_pfpr, tol_rel_change_pfpr, tol_cycle_rmse_rel_pfpr) +
              contribution(m_M,    tol_rel_change_M,    tol_cycle_rmse_rel_M)
    pass   <- m_pfpr$pass_rel && m_pfpr$pass_rmse &&
              m_M$pass_rel    && m_M$pass_rmse
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (verbose) cat(sprintf(
      "    pfpr rel_change=%.4f (%s)  rmse_rel=%.4f (%s); M rel=%.4f (%s) rmse_rel=%.4f (%s); score=%.4f; overall %s; %.1f s\n",
      m_pfpr$rel_change, ifelse(m_pfpr$pass_rel, "PASS", "FAIL"),
      m_pfpr$cycle_rmse_rel, ifelse(m_pfpr$pass_rmse, "PASS", "FAIL"),
      m_M$rel_change, ifelse(m_M$pass_rel, "PASS", "FAIL"),
      m_M$cycle_rmse_rel, ifelse(m_M$pass_rmse, "PASS", "FAIL"),
      score, ifelse(pass, "PASS", "FAIL"), elapsed
    ))
    scored[[k]] <- list(
      snapshot_index = snap$metadata$snapshot_index,
      timesteps_from_t0 = snap$metadata$timesteps_from_t0,
      pfpr = m_pfpr, total_M = m_M,
      score = score, pass = pass,
      elapsed_seconds = elapsed
    )
  }
  list(
    snapshots_scored = scored,
    tolerances = list(
      tol_rel_change_pfpr = tol_rel_change_pfpr,
      tol_rel_change_M    = tol_rel_change_M,
      tol_cycle_rmse_rel_pfpr = tol_cycle_rmse_rel_pfpr,
      tol_cycle_rmse_rel_M    = tol_cycle_rmse_rel_M
    ),
    age_band_days = age_band_days,
    validation_years = validation_years
  )
}


# ------------------------------------------------------------------------------
# promote_best_snapshot
# ------------------------------------------------------------------------------
#' Pick the best-scoring snapshot from a scored library.
#'
#' Selection rule: among snapshots that PASS the stationarity battery, the
#' one with the lowest score. If scores tie, prefer the later snapshot so
#' the promoted checkpoint has more burn-in. If none passes, use the same
#' lowest-score/latest-snapshot rule across all snapshots.
#'
#' @return list with the promoted snapshot, the promotion record, and a
#'   data.frame summarising all candidates.
promote_best_snapshot <- function(library, scored, min_promotion_snapshots = 3L) {
  scores_vec <- vapply(scored$snapshots_scored, function(s) s$score, numeric(1))
  pass_vec   <- vapply(scored$snapshots_scored, function(s) s$pass, logical(1))
  time_vec    <- vapply(scored$snapshots_scored, function(s) s$timesteps_from_t0, integer(1))
  passers    <- which(pass_vec)

  choose_best <- function(candidates) {
    cand_scores <- scores_vec[candidates]
    min_score <- min(cand_scores)
    tied <- candidates[abs(cand_scores - min_score) <= sqrt(.Machine$double.eps)]
    tied[which.max(time_vec[tied])]
  }

  if (length(passers) >= min_promotion_snapshots) {
    chosen_idx <- choose_best(passers)
    quorum_met <- TRUE
  } else {
    chosen_idx <- choose_best(seq_along(scored$snapshots_scored))
    quorum_met <- FALSE
  }

  summary_df <- data.frame(
    snapshot_index    = vapply(scored$snapshots_scored, function(s) s$snapshot_index, integer(1)),
    timesteps_from_t0 = vapply(scored$snapshots_scored, function(s) s$timesteps_from_t0, integer(1)),
    pfpr_rel_change   = vapply(scored$snapshots_scored, function(s) s$pfpr$rel_change, numeric(1)),
    pfpr_cycle_rmse_rel = vapply(scored$snapshots_scored, function(s) s$pfpr$cycle_rmse_rel, numeric(1)),
    M_rel_change      = vapply(scored$snapshots_scored, function(s) s$total_M$rel_change, numeric(1)),
    M_cycle_rmse_rel  = vapply(scored$snapshots_scored, function(s) s$total_M$cycle_rmse_rel, numeric(1)),
    score             = scores_vec,
    pass              = pass_vec,
    chosen            = seq_along(scored$snapshots_scored) == chosen_idx,
    stringsAsFactors = FALSE
  )

  list(
    promoted_snapshot = library$snapshots[[chosen_idx]],
    promoted_index    = chosen_idx,
    quorum_met        = quorum_met,
    n_passers         = length(passers),
    min_promotion_snapshots = min_promotion_snapshots,
    summary           = summary_df,
    tolerances        = scored$tolerances
  )
}
