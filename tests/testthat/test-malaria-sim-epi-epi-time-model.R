testthat::skip_if_not_installed("deSolve")
testthat::skip_if_not_installed("nloptr")
testthat::skip_if_not_installed("MGDrivE")
testthat::skip_if_not_installed("pkgload")
if (!requireNamespace("customMGDrive2", quietly = TRUE)) {
  testthat::skip("customMGDrive2 package is required for malaria-sim-epi epi-time-model tests.")
}

scripts_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
repo_root <- dirname(scripts_root)

required_full_repo_paths <- c(
  file.path(repo_root, "DESCRIPTION"),
  file.path(scripts_root, "landscape-movement", "movement_mu.R"),
  file.path(scripts_root, "archive", "env-codes", "cube.R"),
  file.path(scripts_root, "busia-design-study", "monitoring-design", "gp-baseline", "fit_epi_gp.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_forward_model.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_joint_fit.R")
)
if (!all(file.exists(required_full_repo_paths))) {
  testthat::skip(
    paste(
      "Full customMGDrive2 malaria-sim-epi integration fixtures are not",
      "bundled with this isolated package share."
    )
  )
}

pkgload::load_all(repo_root, quiet = TRUE, export_all = TRUE)
pkgload::load_all(file.path(scripts_root, "msimGD"), quiet = TRUE, export_all = TRUE)
source(file.path(scripts_root, "landscape-movement", "movement_mu.R"), local = TRUE)
source(file.path(scripts_root, "archive", "env-codes", "cube.R"), local = TRUE)
source(file.path(scripts_root, "busia-design-study", "monitoring-design", "gp-baseline", "fit_epi_gp.R"), local = TRUE)
source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_forward_model.R"), local = TRUE)
source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_joint_fit.R"), local = TRUE)

build_tp13_test_cube <- function() {
  mating_comp <- c(
    BB = 1, BH = 0.9, BR = 1, BW = 1, HH = 0.9,
    HR = 0.9, HW = 0.9, RR = 1, RW = 1, WW = 1
  )
  eta_mc <- split(cbind(names(mating_comp), mating_comp), seq_along(mating_comp))
  cube <- cubeTP13(
    gtype = c("BB", "BH", "BR", "BW", "HH", "HR", "HW", "RR", "RW", "WW"),
    p1_germline_M = 0, p2_germline_M = 0.9787776157, p3_germline_M = 0,
    p1_germline_F = 0, p2_germline_F = 0.9850854772, p3_germline_F = 0,
    p1_MDHH = 0.0437169674, p2_MDHH = 0, p3_MDHH = 0.0127526596,
    p1_MDH = 0.0017235774, p2_MDH = 0, p3_MDH = 0.0006749124,
    fc = mating_comp,
    eta = eta_mc, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL
  )
  cube$c <- setNames(rep(1, cube$genotypesN), cube$genotypesID)
  cube
}

theta_test <- list(
  qE = 1 / 3, nE = 2L, qL = 1 / 7, nL = 3L, qP = 1, nP = 2L,
  muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.132, muM = 0.132,
  beta = 16, nu = 1 / (4 / 24)
)

test_that("seasonal Fourier epi basis respects checkpoint phase and harmonic count", {
  times <- c(0, 10, 20, 30)
  basis <- .gp_build_time_basis(
    time = times,
    mode = "seasonal_fourier",
    checkpoint_phase_day = 5L,
    seasonal_harmonics = 1L,
    seasonal_cycle_days = 365L
  )

  expect_identical(basis$info$mode, "seasonal_fourier")
  expect_identical(basis$info$n_columns, 2L)
  expect_identical(colnames(basis$X), c("season_cos1", "season_sin1"))

  phase <- (5 + times) %% 365
  expect_equal(
    as.numeric(basis$X[, "season_cos1"]),
    cos(2 * pi * phase / 365),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(basis$X[, "season_sin1"]),
    sin(2 * pi * phase / 365),
    tolerance = 1e-12
  )
})

test_that("seasonal Fourier epi basis is lower-dimensional than time dummies", {
  times <- rep(seq(0, 84, by = 14), each = 2L)
  basis_dummy <- .gp_build_time_basis(time = times, mode = "time_dummies")
  basis_season <- .gp_build_time_basis(
    time = times,
    mode = "seasonal_fourier",
    checkpoint_phase_day = 0L,
    seasonal_harmonics = 1L,
    seasonal_cycle_days = 365L
  )

  expect_gt(basis_dummy$info$n_columns, basis_season$info$n_columns)
  expect_identical(basis_season$info$n_columns, 2L)
})

test_that("epi node covariates are transformed once at node level", {
  cov_df <- data.frame(
    node_index = 1:3,
    open_eaves_fraction = c(0.55, 0.70, 0.82),
    windows_mean = c(1.8, 2.3, 3.0)
  )

  prepared <- msimGD_prepare_epi_node_covariates(cov_df)

  expect_identical(prepared$node_col, "node_index")
  expect_identical(prepared$extra_effect_names, c("z_eaves", "z_windows"))
  expect_equal(mean(prepared$data$z_eaves), 0, tolerance = 1e-12)
  expect_equal(stats::sd(prepared$data$z_eaves), 1, tolerance = 1e-12)
  expect_equal(mean(prepared$data$z_windows), 0, tolerance = 1e-12)
  expect_equal(stats::sd(prepared$data$z_windows), 1, tolerance = 1e-12)
})

test_that("msimGD_joint_fit supports the seasonal Fourier epi time model", {
  set.seed(101)

  cube <- build_tp13_test_cube()
  NF <- c(35, 55, 45)
  NH <- c(100L, 100L, 100L)
  setup <- list(D = matrix(c(
    0, 1, 5,
    1, 0, 4,
    5, 4, 0
  ), nrow = 3, byrow = TRUE))
  release <- list(nodes = 1L, time = 3L, size = 200L, stage = "M", genotype = "HH")

  prebuilt <- msimGD_precompute_base(
    cube = cube,
    NF = NF,
    NH = NH,
    theta = theta_test,
    init_EIR = 5
  )

  design <- expand.grid(node = 1:3, time = seq(0, 28, by = 7))
  pred <- msimGD_det_predict_carrier(
    setup = setup,
    design = design,
    mu = 2.36,
    p_move = 0.04,
    tmax = 28L,
    release = release,
    prebuilt_base = prebuilt
  )$pred

  ento_obs <- within(pred, {
    n <- 2000L
    x <- stats::rbinom(length(carrier_freq), size = n, prob = pmax(1e-6, pmin(1 - 1e-6, carrier_freq)))
  })
  ento_obs <- ento_obs[, c("node", "time", "n", "x")]

  epi_df <- pred
  epi_df$w <- 14L
  epi_df$offset_pop <- NH[epi_df$node]
  phase <- (17 + epi_df$time) %% 365
  seasonal_effect <- 0.4 * cos(2 * pi * phase / 365) - 0.2 * sin(2 * pi * phase / 365)
  lambda <- exp(log(epi_df$w * epi_df$offset_pop) - 4.2 + 2.0 * epi_df$carrier_freq + seasonal_effect)
  epi_df$Y <- stats::rpois(nrow(epi_df), lambda)
  epi_obs <- epi_df[, c("node", "time", "w", "Y", "offset_pop")]

  fit <- msimGD_joint_fit(
    ento_obs = ento_obs,
    epi_obs = epi_obs,
    setup = setup,
    release = release,
    prebuilt_base = prebuilt,
    NH = NH,
    tmax = 28L,
    bounds_mu = c(2.359, 2.361),
    bounds_p_move = c(0.039, 0.041),
    start_mu = 2.36,
    start_p_move = 0.04,
    start_kappa = 50,
    start_log_tau = log(0.5),
    theta_gp_fixed = 8.0,
    fit_mode = "joint",
    epi_time_model = "seasonal_fourier",
    checkpoint_phase_day = 17L,
    seasonal_harmonics = 1L,
    maxeval = 200,
    verbose = FALSE,
    D_monitor = as.matrix(setup$D)
  )

  expect_identical(fit$epi_time_model, "seasonal_fourier")
  expect_identical(fit$checkpoint_phase_day, 17L)
  expect_identical(fit$seasonal_harmonics, 1L)
  expect_identical(fit$opt$convergence, 0L)
  expect_match(fit$model_name, "seasonal_fourier")
})

test_that("msimGD_joint_fit includes node covariates in the epi GP fixed effects", {
  set.seed(102)

  cube <- build_tp13_test_cube()
  NF <- c(35, 55, 45)
  NH <- c(100L, 100L, 100L)
  setup <- list(D = matrix(c(
    0, 1, 5,
    1, 0, 4,
    5, 4, 0
  ), nrow = 3, byrow = TRUE))
  release <- list(nodes = 1L, time = 3L, size = 200L, stage = "M", genotype = "HH")

  prebuilt <- msimGD_precompute_base(
    cube = cube,
    NF = NF,
    NH = NH,
    theta = theta_test,
    init_EIR = 5
  )

  design <- expand.grid(node = 1:3, time = seq(0, 28, by = 7))
  pred <- msimGD_det_predict_carrier(
    setup = setup,
    design = design,
    mu = 2.36,
    p_move = 0.04,
    tmax = 28L,
    release = release,
    prebuilt_base = prebuilt
  )$pred

  ento_obs <- within(pred, {
    n <- 2000L
    x <- stats::rbinom(length(carrier_freq), size = n, prob = pmax(1e-6, pmin(1 - 1e-6, carrier_freq)))
  })
  ento_obs <- ento_obs[, c("node", "time", "n", "x")]

  cov_df <- data.frame(
    node = 1:3,
    open_eaves_fraction = c(0.55, 0.70, 0.82),
    windows_mean = c(1.8, 2.3, 3.0)
  )
  cov_prepared <- msimGD_prepare_epi_node_covariates(cov_df)$data

  epi_df <- pred
  epi_df$w <- 14L
  epi_df$offset_pop <- NH[epi_df$node]
  phase <- (17 + epi_df$time) %% 365
  seasonal_effect <- 0.4 * cos(2 * pi * phase / 365) - 0.2 * sin(2 * pi * phase / 365)
  cov_match <- match(epi_df$node, cov_prepared$node)
  lambda <- exp(
    log(epi_df$w * epi_df$offset_pop) - 4.2 +
      2.0 * epi_df$carrier_freq +
      0.35 * cov_prepared$z_eaves[cov_match] -
      0.25 * cov_prepared$z_windows[cov_match] +
      seasonal_effect
  )
  epi_df$Y <- stats::rpois(nrow(epi_df), lambda)
  epi_obs <- epi_df[, c("node", "time", "w", "Y", "offset_pop")]

  fit <- msimGD_joint_fit(
    ento_obs = ento_obs,
    epi_obs = epi_obs,
    epi_node_covariates = cov_df,
    setup = setup,
    release = release,
    prebuilt_base = prebuilt,
    NH = NH,
    tmax = 28L,
    bounds_mu = c(2.359, 2.361),
    bounds_p_move = c(0.039, 0.041),
    start_mu = 2.36,
    start_p_move = 0.04,
    start_kappa = 50,
    start_log_tau = log(0.5),
    theta_gp_fixed = 8.0,
    fit_mode = "joint",
    epi_time_model = "seasonal_fourier",
    checkpoint_phase_day = 17L,
    seasonal_harmonics = 1L,
    maxeval = 200,
    verbose = FALSE,
    D_monitor = as.matrix(setup$D)
  )

  expect_identical(fit$opt$convergence, 0L)
  expect_match(fit$model_name, "epi_covariates")
  expect_false(is.null(fit$epi_node_covariates))
  expect_identical(
    fit$epi_node_covariates$extra_effect_names,
    c("z_eaves", "z_windows")
  )
})
