testthat::skip_if_not_installed("pkgload")

scripts_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)

required_full_repo_paths <- c(
  file.path(scripts_root, "msimGD", "DESCRIPTION"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_seasonality_utils.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_calibration_spec.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_baseline_calibration_utils.R")
)
if (!all(file.exists(required_full_repo_paths))) {
  testthat::skip(
    paste(
      "Full customMGDrive2 malaria-sim-epi integration fixtures are not",
      "bundled with this isolated package share."
    )
  )
}

pkgload::load_all(file.path(scripts_root, "msimGD"), quiet = TRUE, export_all = TRUE)
source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_seasonality_utils.R"), local = TRUE)
source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_calibration_spec.R"), local = TRUE)
source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_baseline_calibration_utils.R"), local = TRUE)


test_that("Ochomo Fourier seasonality fitter reproduces the synthetic rainfall cycle", {
  truth <- ochomo_fourier_rainfall_seasonality(
    g0 = 1.2,
    g = c(0.25, -0.06, 0.03),
    h = c(0.10, 0.04, -0.02),
    rainfall_floor = 0.001
  )
  days <- 1:365
  target <- vapply(
    days,
    function(day) {
      rainfall(
        t = as.integer(day),
        g0 = truth$g0,
        g = truth$g,
        h = truth$h,
        floor = truth$rainfall_floor
      )
    },
    numeric(1)
  )

  fit <- ochomo_fit_fourier_rainfall_seasonality(
    daily_rainfall = target,
    day = days,
    harmonics = 3L,
    rainfall_floor = truth$rainfall_floor
  )

  expect_lt(fit$rmse, 1e-6)
  expect_equal(
    fit$fitted_daily_cycle$fitted,
    fit$fitted_daily_cycle$target,
    tolerance = 1e-6
  )
})


test_that("Ochomo baseline modifier applies configured seasonality and keeps the default spec flat", {
  spec_default <- ochomo_calibration_spec()
  modifier_default <- ochomo_build_baseline_parameter_modifier(
    spec = spec_default,
    treatment_coverage = 0,
    detectability_d1 = spec_default$microscopy_detection$d1
  )
  params_default <- modifier_default(
    get_parameters(list(progress_bar = FALSE)),
    node_index = 1L,
    warmup_days = 365L
  )

  expect_identical(spec_default$release_ready_baseline$seasonality$type, "none")
  expect_false(isTRUE(params_default$model_seasonality))

  spec_seasonal <- ochomo_calibration_spec()
  seasonal_cfg <- ochomo_fourier_rainfall_seasonality(
    g0 = 1.5,
    g = c(0.22, 0.05, -0.03),
    h = c(0.08, -0.02, 0.01),
    rainfall_floor = 0.005
  )
  spec_seasonal$release_ready_baseline$explicit_bednets <- FALSE
  spec_seasonal$release_ready_baseline$seasonality <- seasonal_cfg

  modifier_seasonal <- ochomo_build_baseline_parameter_modifier(
    spec = spec_seasonal,
    treatment_coverage = 0,
    detectability_d1 = spec_seasonal$microscopy_detection$d1
  )
  params_seasonal <- modifier_seasonal(
    get_parameters(list(progress_bar = FALSE)),
    node_index = 1L,
    warmup_days = 365L
  )

  expect_true(isTRUE(params_seasonal$model_seasonality))
  expect_equal(params_seasonal$g0, seasonal_cfg$g0, tolerance = 1e-12)
  expect_equal(params_seasonal$g, seasonal_cfg$g, tolerance = 1e-12)
  expect_equal(params_seasonal$h, seasonal_cfg$h, tolerance = 1e-12)
  expect_equal(params_seasonal$rainfall_floor, seasonal_cfg$rainfall_floor, tolerance = 1e-12)
})


test_that("Ochomo seasonality resolver builds a shared Fourier forcing from CSV rainfall data", {
  days <- 1:365
  truth <- ochomo_fourier_rainfall_seasonality(
    g0 = 1.1,
    g = c(0.18, -0.04, 0.02),
    h = c(0.07, 0.03, -0.01),
    rainfall_floor = 0.001
  )
  daily_target <- vapply(
    days,
    function(day) {
      rainfall(
        t = as.integer(day),
        g0 = truth$g0,
        g = truth$g,
        h = truth$h,
        floor = truth$rainfall_floor
      )
    },
    numeric(1)
  )
  synthetic <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 365 * 2L),
    rain_mm = rep(daily_target, 2L)
  )
  csv_path <- tempfile(fileext = ".csv")
  utils::write.csv(synthetic, csv_path, row.names = FALSE)
  on.exit(unlink(csv_path), add = TRUE)

  resolved <- ochomo_resolve_release_ready_seasonality(
    csv_path = csv_path,
    date_col = "date",
    rainfall_col = "rain_mm",
    harmonics = 3L,
    rainfall_floor = truth$rainfall_floor
  )

  expect_identical(resolved$type, "fourier_rainfall")
  expect_lt(resolved$rmse, 2e-4)
  expect_match(resolved$source, "csv:", fixed = TRUE)
})
