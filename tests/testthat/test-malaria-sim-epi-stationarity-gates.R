scripts_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

required_full_repo_paths <- c(
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_seasonality_utils.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_baseline_calibration_utils.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_calibration_spec.R")
)
if (!all(file.exists(required_full_repo_paths))) {
  testthat::skip(
    paste(
      "Full customMGDrive2 malaria-sim-epi integration fixtures are not",
      "bundled with this isolated package share."
    )
  )
}

source(file.path(
  scripts_root,
  "msimGD",
  "test-scripts",
  "malaria-sim-epi",
  "ochomo_seasonality_utils.R"
), local = TRUE)
source(file.path(
  scripts_root,
  "msimGD",
  "test-scripts",
  "malaria-sim-epi",
  "ochomo_baseline_calibration_utils.R"
), local = TRUE)
source(file.path(
  scripts_root,
  "msimGD",
  "test-scripts",
  "malaria-sim-epi",
  "ochomo_calibration_spec.R"
), local = TRUE)


build_stationarity_gate_spec <- function(mode = c("flat", "annual_cycle")) {
  mode <- match.arg(mode)
  list(
    stationarity = list(
      mode = mode,
      cycle_days = 365L,
      gate_policy = "practical_effect_size",
      gate_scope = "aggregate_only",
      strict_aggregate_metrics = c("total_M_per_person", "target_prevalence"),
      strict_node_metrics = character(0),
      strict_mean_slope_metrics = character(0),
      strict_sign_test_alpha = 0.05,
      strict_max_abs_mean_slope_rel_per_year = 0.01,
      strict_max_nodes_all_same_direction = 0L,
      practical_effect_metrics = c("total_M_per_person", "target_prevalence"),
      practical_max_abs_mean_rel_change = c(
        total_M_per_person = 0.05,
        target_prevalence = 0.02
      ),
      practical_max_cycle_rmse_rel = c(
        total_M_per_person = 0.10,
        target_prevalence = 0.05
      )
    )
  )
}


test_that("annual-cycle stationarity summary aligns repeated yearly phases", {
  spec <- build_stationarity_gate_spec("annual_cycle")
  time <- seq(14, 728, by = 14)
  value <- 10 + sin(2 * pi * (time - 14) / 365)

  annual <- ochomo_stationarity_metric_summary(time, value, spec)
  flat <- ochomo_stationarity_trend_summary(time, value)

  testthat::expect_lt(abs(annual$rel_change), 0.01)
  testthat::expect_lt(annual$cycle_rmse_rel, 0.01)
  testthat::expect_gt(abs(flat$rel_change), 0.02)
  testthat::expect_gt(annual$cycle_overlap_n, 20L)
})


test_that("annual-cycle gate uses cycle-aligned thresholds", {
  spec <- build_stationarity_gate_spec("annual_cycle")
  battery <- list(
    aggregate_summary = data.frame(
      metric = c("total_M_per_person", "target_prevalence"),
      n_seeds = c(24L, 24L),
      slope_pos = c(0L, 0L),
      slope_neg = c(0L, 0L),
      sign_test_p = c(NA_real_, NA_real_),
      mean_rel_change = c(0.03, 0.01),
      median_rel_change = c(0.03, 0.01),
      mean_slope_rel_per_year = c(NA_real_, NA_real_),
      median_slope_rel_per_year = c(NA_real_, NA_real_),
      ci_low_slope_rel_per_year = c(NA_real_, NA_real_),
      ci_high_slope_rel_per_year = c(NA_real_, NA_real_),
      mean_spearman_rho = c(NA_real_, NA_real_),
      mean_cycle_rmse_rel = c(0.08, 0.03),
      median_cycle_rmse_rel = c(0.08, 0.03),
      mean_cycle_phase_cor = c(0.99, 0.99),
      stringsAsFactors = FALSE
    ),
    node_consistency = data.frame(metric = character(0)),
    node_summary = data.frame(metric = character(0))
  )

  gates <- ochomo_evaluate_stationarity_gates(battery, spec)
  testthat::expect_true(gates$overall_pass)

  battery$aggregate_summary$mean_cycle_rmse_rel[[1L]] <- 0.15
  gates_fail <- ochomo_evaluate_stationarity_gates(battery, spec)
  testthat::expect_false(gates_fail$overall_pass)
  testthat::expect_false(gates_fail$aggregate$pass[[1L]])
})


test_that("Ochomo seasonal spec switches checkpoint acceptance to annual-cycle mode", {
  old_csv <- Sys.getenv("MSIMGD_OCHOMO_SEASONALITY_CSV", unset = NA_character_)
  old_harmonics <- Sys.getenv("MSIMGD_OCHOMO_SEASONALITY_HARMONICS", unset = NA_character_)
  on.exit({
    if (is.na(old_csv)) {
      Sys.unsetenv("MSIMGD_OCHOMO_SEASONALITY_CSV")
    } else {
      Sys.setenv(MSIMGD_OCHOMO_SEASONALITY_CSV = old_csv)
    }
    if (is.na(old_harmonics)) {
      Sys.unsetenv("MSIMGD_OCHOMO_SEASONALITY_HARMONICS")
    } else {
      Sys.setenv(MSIMGD_OCHOMO_SEASONALITY_HARMONICS = old_harmonics)
    }
  }, add = TRUE)

  csv_path <- tempfile(fileext = ".csv")
  rainfall_df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), as.Date("2021-12-31"), by = "day"),
    rainfall_mm = rep(c(0.5, 1.0, 2.0, 3.0), length.out = 731L)
  )
  utils::write.csv(rainfall_df, csv_path, row.names = FALSE)

  Sys.setenv(
    MSIMGD_OCHOMO_SEASONALITY_CSV = csv_path,
    MSIMGD_OCHOMO_SEASONALITY_HARMONICS = "3"
  )
  spec <- ochomo_calibration_spec()

  testthat::expect_identical(spec$stationarity$mode, "annual_cycle")
  testthat::expect_equal(
    spec$stationarity$strict_aggregate_metrics,
    c("total_M_per_person", "target_prevalence")
  )
  testthat::expect_true(all(
    c("total_M_per_person", "target_prevalence") %in%
      names(spec$stationarity$practical_max_abs_mean_rel_change)
  ))
  testthat::expect_true(all(
    c("total_M_per_person", "target_prevalence") %in%
      names(spec$stationarity$practical_max_cycle_rmse_rel)
  ))
})
