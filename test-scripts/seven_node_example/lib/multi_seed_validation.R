# ------------------------------------------------------------------------------
# multi_seed_validation.R
# ------------------------------------------------------------------------------
# Slim faithful re-implementation of the production multi-seed stationarity
# validation step. Re-runs the promoted baseline checkpoint under M distinct
# RNG seeds for `validation_years` years each and verifies that the
# stationarity metrics (Y/Y rel_change, cycle RMSE) computed earlier are
# stable across seeds.
#
# "Stable across seeds" gate: across the M reruns, the mean rel_change per
# signal must satisfy the same tolerance as the per-snapshot battery, AND
# the SD of rel_change across seeds must be below tol_seed_sd_rel (default
# 0.02). This is the seed-robustness check.
# ------------------------------------------------------------------------------

#' Re-run the promoted snapshot under M seeds and confirm stationarity holds.
#'
#' @param promoted_snapshot one of the elements of `library$snapshots`
#'   (a list with `state` and `metadata`)
#' @param n_seeds integer M; number of seeds to retry under (default 3)
#' @param validation_years integer (default 2)
#' @param age_band_days integer length-2 (default c(730, 3650))
#' @param tol_rel_change_pfpr, tol_rel_change_M, tol_cycle_rmse_rel_pfpr,
#'   tol_cycle_rmse_rel_M same as the per-snapshot battery
#' @param tol_seed_sd_rel tolerance on across-seed SD of rel_change (default 0.02)
#' @param master_seed_offset offset added to promoted_snapshot$metadata$seed
#'   for each seed; default 100000 (large to avoid colliding with the
#'   library-build seeds)
#' @return list with per-seed metrics, across-seed summary, and pass/fail
multi_seed_validate_promoted_snapshot <- function(
    promoted_snapshot,
    setup, cube, NH, mu, p_move, theta, parameter_modifier,
    n_seeds = 3L, validation_years = 2L,
    age_band_days = c(730L, 3650L),
    tol_rel_change_pfpr = 0.02, tol_rel_change_M = 0.05,
    tol_cycle_rmse_rel_pfpr = 0.10, tol_cycle_rmse_rel_M = 0.10,
    tol_seed_sd_rel = 0.02,
    master_seed_offset = 100000L,
    verbose = TRUE) {
  for (h in c("msimGD_run_truth")) {
    if (!exists(h, mode = "function")) {
      stop(sprintf("%s() not in scope. Source lib/msimGD_truth_generation.R first.", h),
           call. = FALSE)
    }
  }
  validation_days <- as.integer(validation_years) * 365L

  # Per-seed runs
  per_seed <- vector("list", n_seeds)
  base_seed <- promoted_snapshot$metadata$seed
  for (s in seq_len(n_seeds)) {
    this_seed <- as.integer(base_seed + master_seed_offset + s)
    if (verbose) cat(sprintf(
      "[%s] seed %d/%d (= %d)...\n", format(Sys.time(), "%F %T"),
      s, n_seeds, this_seed
    ))
    t0 <- Sys.time()
    res <- msimGD_run_truth(
      setup = setup, cube = cube, NF = NULL, NH = NH,
      tmax = validation_days, mu = mu, p_move = p_move, release = NULL,
      theta = theta, init_EIR = promoted_snapshot$metadata$init_EIR,
      prevalence_rendering_min_age = age_band_days[[1]],
      prevalence_rendering_max_age = age_band_days[[2]],
      warmup_days = 0L,
      parameter_modifier = parameter_modifier,
      baseline_checkpoint = list(state = promoted_snapshot$state,
                                  metadata = promoted_snapshot$metadata),
      seed = this_seed
    )
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    data_list <- if (is.list(res) && all(vapply(res, is.data.frame, logical(1)))) res
                 else if (!is.null(res$data)) res$data
                 else stop("Could not locate data list.", call. = FALSE)

    n_col <- sprintf("n_detect_lm_%d_%d", age_band_days[[1]], age_band_days[[2]])
    d_col <- sprintf("n_age_%d_%d",       age_band_days[[1]], age_band_days[[2]])
    combined <- do.call(rbind, data_list)
    agg_p <- aggregate(combined[, c(n_col, d_col)],
                       by = list(t = combined$timestep), FUN = sum)
    pfpr_series <- agg_p[[n_col]] / pmax(agg_p[[d_col]], 1)
    agg_M <- aggregate(combined$total_M_gamb,
                       by = list(t = combined$timestep), FUN = sum)
    M_series <- agg_M$x

    yN_p   <- pfpr_series[(validation_days - 365L + 1L):validation_days]
    yN1_p  <- pfpr_series[1L:365L]
    yN_M   <- M_series[(validation_days - 365L + 1L):validation_days]
    yN1_M  <- M_series[1L:365L]

    per_seed[[s]] <- list(
      seed = this_seed,
      pfpr_rel_change = abs(mean(yN_p) - mean(yN1_p)) / max(abs(mean(yN1_p)), .Machine$double.eps),
      pfpr_cycle_rmse_rel = sqrt(mean((yN_p - yN1_p)^2)) /
        max(abs(mean(yN1_p)), .Machine$double.eps),
      M_rel_change = abs(mean(yN_M) - mean(yN1_M)) / max(abs(mean(yN1_M)), .Machine$double.eps),
      M_cycle_rmse_rel = sqrt(mean((yN_M - yN1_M)^2)) /
        max(abs(mean(yN1_M)), .Machine$double.eps),
      elapsed_seconds = elapsed
    )
  }

  # Across-seed aggregates
  pf_rc   <- vapply(per_seed, `[[`, numeric(1), "pfpr_rel_change")
  pf_rmse <- vapply(per_seed, `[[`, numeric(1), "pfpr_cycle_rmse_rel")
  M_rc    <- vapply(per_seed, `[[`, numeric(1), "M_rel_change")
  M_rmse  <- vapply(per_seed, `[[`, numeric(1), "M_cycle_rmse_rel")

  pass <- (mean(pf_rc) <= tol_rel_change_pfpr) &&
          (mean(M_rc)  <= tol_rel_change_M) &&
          (mean(pf_rmse) <= tol_cycle_rmse_rel_pfpr) &&
          (mean(M_rmse)  <= tol_cycle_rmse_rel_M) &&
          (sd(pf_rc) <= tol_seed_sd_rel) &&
          (sd(M_rc)  <= tol_seed_sd_rel)

  summary_df <- data.frame(
    metric = c("pfpr_rel_change", "pfpr_cycle_rmse_rel",
               "M_rel_change",    "M_cycle_rmse_rel"),
    across_seed_mean = c(mean(pf_rc), mean(pf_rmse), mean(M_rc), mean(M_rmse)),
    across_seed_sd   = c(sd(pf_rc),   sd(pf_rmse),   sd(M_rc),   sd(M_rmse)),
    tol_mean = c(tol_rel_change_pfpr, tol_cycle_rmse_rel_pfpr,
                 tol_rel_change_M,    tol_cycle_rmse_rel_M),
    tol_sd   = c(tol_seed_sd_rel, NA_real_, tol_seed_sd_rel, NA_real_),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat(sprintf("[%s] multi-seed validation: %s\n",
                format(Sys.time(), "%F %T"),
                ifelse(pass, "PASS", "FAIL")))
    print(summary_df)
  }

  list(
    pass = pass,
    n_seeds = as.integer(n_seeds),
    per_seed = per_seed,
    summary = summary_df,
    tolerances = list(
      tol_rel_change_pfpr = tol_rel_change_pfpr,
      tol_rel_change_M    = tol_rel_change_M,
      tol_cycle_rmse_rel_pfpr = tol_cycle_rmse_rel_pfpr,
      tol_cycle_rmse_rel_M    = tol_cycle_rmse_rel_M,
      tol_seed_sd_rel     = tol_seed_sd_rel
    )
  )
}
