testthat::skip_if_not_installed("deSolve")
testthat::skip_if_not_installed("nloptr")
testthat::skip_if_not_installed("MGDrivE")
testthat::skip_if_not_installed("pkgload")
if (!requireNamespace("customMGDrive2", quietly = TRUE)) {
  testthat::skip("customMGDrive2 package is required for malaria-sim-epi forward-model tests.")
}

scripts_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
repo_root <- dirname(scripts_root)

required_full_repo_paths <- c(
  file.path(repo_root, "DESCRIPTION"),
  file.path(scripts_root, "landscape-movement", "movement_mu.R"),
  file.path(scripts_root, "archive", "env-codes", "cube.R"),
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

test_that("relative abundance surfaces are mean-one relative multipliers", {
  a_i <- msimGD_normalise_relative_abundance_surface(c(2, 4, 6), n_nodes = 3L)
  expect_equal(mean(a_i), 1, tolerance = 1e-12)
  expect_equal(sum(a_i), 3, tolerance = 1e-12)
  expect_equal(as.numeric(a_i), c(0.5, 1, 1.5), tolerance = 1e-12)
})

test_that("scaled-surface baseline reconstructs node abundances from s_N * a_i", {
  cube <- build_tp13_test_cube()
  NF <- c(30, 60)
  NH <- c(100L, 100L)

  prebuilt <- msimGD_precompute_base(
    cube = cube,
    NF = NF,
    NH = NH,
    theta = theta_test,
    init_EIR = 5
  )
  a_i <- msimGD_relative_abundance_from_total_M(NF, n_nodes = length(NF))
  rebuilt <- msimGD_prebuilt_base_with_scaled_surface(
    prebuilt_base = prebuilt,
    s_N = mean(NF),
    a_i = a_i,
    baseline_surface_source = "test_surface"
  )

  expect_equal(unname(rebuilt$NF), NF, tolerance = 1e-12)
  expect_equal(
    vapply(rebuilt$node_params, function(bp) as.numeric(bp$total_M), numeric(1)),
    NF,
    tolerance = 1e-12
  )
  expect_equal(mean(rebuilt$baseline_surface), 1, tolerance = 1e-12)
  expect_identical(rebuilt$baseline_mode, "scaled_surface")
})

test_that("forward precompute records but does not mechanistically propagate contact_multiplier", {
  cube <- build_tp13_test_cube()
  NF <- c(30, 60)
  NH <- c(100L, 100L)

  prebuilt <- msimGD_precompute_base(
    cube = cube,
    NF = NF,
    NH = NH,
    theta = theta_test,
    init_EIR = 5,
    parameter_modifier = function(parameters, node_index, warmup_days) {
      parameters$contact_multiplier <- c(0.8, 1.2)[[node_index]]
      parameters
    }
  )

  expect_equal(prebuilt$contact_multiplier_by_node, c(0.8, 1.2), tolerance = 1e-12)
  expect_false(prebuilt$contact_multiplier_supported)
})

test_that("joint s_N estimation recovers the oracle mean under exact validation surface", {
  cube <- build_tp13_test_cube()
  NF <- c(30, 60, 45)
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
  a_i <- msimGD_relative_abundance_from_total_M(NF, n_nodes = length(NF))
  true_base <- msimGD_prebuilt_base_with_scaled_surface(
    prebuilt_base = prebuilt,
    s_N = mean(NF),
    a_i = a_i,
    baseline_surface_source = "simulation_validation_oracle_nf_eq"
  )

  design <- expand.grid(node = 1:3, time = 0:10)
  pred <- msimGD_det_predict_carrier(
    setup = setup,
    design = design,
    mu = 2.36,
    p_move = 0.04,
    tmax = 10L,
    release = release,
    prebuilt_base = true_base
  )$pred

  ento_obs <- within(pred, {
    n <- 100000L
    x <- pmax(0L, pmin(n, as.integer(round(n * carrier_freq))))
  })
  ento_obs <- ento_obs[, c("node", "time", "n", "x")]

  fit <- msimGD_joint_fit(
    ento_obs = ento_obs,
    epi_obs = NULL,
    setup = setup,
    release = release,
    prebuilt_base = prebuilt,
    NH = NH,
    tmax = 10L,
    bounds_mu = c(2.3599, 2.3601),
    bounds_p_move = c(0.0399, 0.0401),
    start_mu = 2.36,
    start_p_move = 0.04,
    start_kappa = 50,
    start_s_N = mean(NF),
    baseline_mode = "scaled_surface",
    baseline_surface = a_i,
    baseline_surface_source = "simulation_validation_oracle_nf_eq",
    bounds_log_s_N = log(mean(NF)) + log(c(0.5, 2)),
    maxeval = 100,
    fit_mode = "ento_only",
    verbose = FALSE
  )

  expect_true(is.finite(fit$est$s_N) && fit$est$s_N > 0)
  expect_equal(
    as.numeric(fit$baseline$total_M),
    as.numeric(fit$est$s_N * a_i),
    tolerance = 1e-8
  )
  expect_identical(fit$baseline$mode, "scaled_surface")
  expect_identical(
    fit$baseline$source,
    "simulation_validation_oracle_nf_eq"
  )
})

test_that("forward carrier predictions respect time-varying carrying capacity after day 0", {
  testthat::skip_if_not_installed("callr")

  maxdiff <- callr::r(
    function(repo_root, scripts_root) {
      library(pkgload)

      pkgload::load_all(repo_root, quiet = TRUE, export_all = TRUE)
      pkgload::load_all(file.path(scripts_root, "msimGD"), quiet = TRUE, export_all = TRUE)
      source(file.path(scripts_root, "landscape-movement", "movement_mu.R"), local = TRUE)
      source(file.path(scripts_root, "archive", "env-codes", "cube.R"), local = TRUE)
      source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_forward_model.R"), local = TRUE)

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

      cube <- build_tp13_test_cube()
      theta_test <- list(
        qE = 1 / 3, nE = 2L, qL = 1 / 7, nL = 3L, qP = 1, nP = 2L,
        muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.132, muM = 0.132,
        beta = 16, nu = 1 / (4 / 24)
      )
      NF <- c(40, 55)
      NH <- c(100L, 100L)
      setup <- list(D = matrix(c(
        0, 1,
        1, 0
      ), nrow = 2, byrow = TRUE))
      release <- list(nodes = 1L, time = 3L, size = 20000L, stage = "M", genotype = "HH")
      design <- expand.grid(node = 1:2, time = 0:120)

      prebuilt_plain <- msimGD_precompute_base(
        cube = cube,
        NF = NF,
        NH = NH,
        theta = theta_test,
        init_EIR = 5
      )
      prebuilt_shifted <- msimGD_precompute_base(
        cube = cube,
        NF = NF,
        NH = NH,
        theta = theta_test,
        init_EIR = 5,
        parameter_modifier = function(parameters, node_index, warmup_days) {
          set_carrying_capacity(
            parameters = parameters,
            timesteps = 4L,
            carrying_capacity_scalers = matrix(0.01, nrow = 1L, ncol = length(parameters$species))
          )
        }
      )

      pred_plain <- msimGD_det_predict_carrier(
        setup = setup,
        design = design,
        mu = 1.0,
        p_move = 0.04,
        tmax = 120L,
        release = release,
        prebuilt_base = prebuilt_plain
      )$pred
      pred_shifted <- msimGD_det_predict_carrier(
        setup = setup,
        design = design,
        mu = 1.0,
        p_move = 0.04,
        tmax = 120L,
        release = release,
        prebuilt_base = prebuilt_shifted
      )$pred

      post_change <- pred_plain$time >= 4
      stopifnot(any(post_change))
      max(abs(pred_shifted$carrier_freq[post_change] - pred_plain$carrier_freq[post_change]))
    },
    args = list(repo_root = repo_root, scripts_root = scripts_root)
  )

  expect_gt(maxdiff, 1e-3)
})


test_that("forward carrier predictions are periodic in start_phase_day", {
  testthat::skip_if_not_installed("callr")

  phase_check <- callr::r(
    function(repo_root, scripts_root) {
      library(pkgload)

      pkgload::load_all(repo_root, quiet = TRUE, export_all = TRUE)
      pkgload::load_all(file.path(scripts_root, "msimGD"), quiet = TRUE, export_all = TRUE)
      source(file.path(scripts_root, "landscape-movement", "movement_mu.R"), local = TRUE)
      source(file.path(scripts_root, "archive", "env-codes", "cube.R"), local = TRUE)
      source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_forward_model.R"), local = TRUE)

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

      cube <- build_tp13_test_cube()
      theta_test <- list(
        qE = 1 / 3, nE = 2L, qL = 1 / 7, nL = 3L, qP = 1, nP = 2L,
        muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.132, muM = 0.132,
        beta = 16, nu = 1 / (4 / 24)
      )
      NF <- c(40, 55)
      NH <- c(100L, 100L)
      setup <- list(D = matrix(c(
        0, 1,
        1, 0
      ), nrow = 2, byrow = TRUE))
      release <- list(nodes = 1L, time = 3L, size = 4000L, stage = "M", genotype = "HH")
      design <- expand.grid(node = 1:2, time = 0:60)

      prebuilt <- msimGD_precompute_base(
        cube = cube,
        NF = NF,
        NH = NH,
        theta = theta_test,
        init_EIR = 5,
        parameter_modifier = function(parameters, node_index, warmup_days) {
          parameters$model_seasonality <- TRUE
          parameters$g0 <- 2.1
          parameters$g <- c(0.35, -0.18, 0.11)
          parameters$h <- c(0.27, 0.09, -0.16)
          parameters
        }
      )

      pred_phase_0 <- msimGD_det_predict_carrier(
        setup = setup,
        design = design,
        mu = 1.0,
        p_move = 0.04,
        tmax = 60L,
        release = release,
        start_phase_day = 0L,
        prebuilt_base = prebuilt
      )$pred
      pred_phase_17 <- msimGD_det_predict_carrier(
        setup = setup,
        design = design,
        mu = 1.0,
        p_move = 0.04,
        tmax = 60L,
        release = release,
        start_phase_day = 17L,
        prebuilt_base = prebuilt
      )$pred
      pred_phase_382 <- msimGD_det_predict_carrier(
        setup = setup,
        design = design,
        mu = 1.0,
        p_move = 0.04,
        tmax = 60L,
        release = release,
        start_phase_day = 382L,
        prebuilt_base = prebuilt
      )$pred

      list(
        same_cycle_maxdiff = max(abs(pred_phase_17$carrier_freq - pred_phase_382$carrier_freq)),
        shifted_phase_maxdiff = max(abs(pred_phase_0$carrier_freq - pred_phase_17$carrier_freq))
      )
    },
    args = list(repo_root = repo_root, scripts_root = scripts_root)
  )

  expect_lt(phase_check$same_cycle_maxdiff, 1e-9)
  expect_gt(phase_check$shifted_phase_maxdiff, 1e-4)
})
