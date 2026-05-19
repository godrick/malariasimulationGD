testthat::skip_if_not_installed("pkgload")
testthat::skip_if_not_installed("callr")
if (!requireNamespace("customMGDrive2", quietly = TRUE)) {
  testthat::skip("customMGDrive2 package is required for malaria-sim-epi checkpoint tests.")
}

scripts_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
repo_root <- dirname(scripts_root)

required_full_repo_paths <- c(
  file.path(repo_root, "DESCRIPTION"),
  file.path(scripts_root, "landscape-movement", "movement_mu.R"),
  file.path(scripts_root, "archive", "env-codes", "cube.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_truth_generation.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_seasonality_utils.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_baseline_calibration_utils.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_study_start_utils.R")
)
if (!all(file.exists(required_full_repo_paths))) {
  testthat::skip(
    paste(
      "Full customMGDrive2 malaria-sim-epi integration fixtures are not",
      "bundled with this isolated package share."
    )
  )
}

run_checkpoint_subprocess <- function(case_name) {
  callr::r(
    function(repo_root, scripts_root, case_name) {
      library(pkgload)
      setwd(repo_root)

      pkgload::load_all(repo_root, quiet = TRUE, export_all = TRUE)
      pkgload::load_all(file.path(scripts_root, "msimGD"), quiet = TRUE, export_all = TRUE)
      source(file.path(scripts_root, "landscape-movement", "movement_mu.R"), local = TRUE)
      source(file.path(scripts_root, "archive", "env-codes", "cube.R"), local = TRUE)
      source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_truth_generation.R"), local = TRUE)
      source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_seasonality_utils.R"), local = TRUE)
      source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_baseline_calibration_utils.R"), local = TRUE)
      source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_study_start_utils.R"), local = TRUE)

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
        cube$b <- setNames(rep(1, cube$genotypesN), cube$genotypesID)
        cube$c <- setNames(rep(1, cube$genotypesN), cube$genotypesID)
        cube
      }

      theta_test <- list(
        qE = 1 / 3, nE = 2L, qL = 1 / 7, nL = 3L, qP = 1, nP = 2L,
        muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.132, muM = 0.132,
        beta = 16, nu = 1 / (4 / 24)
      )
      setup <- list(D = matrix(c(
        0, 1,
        1, 0
      ), nrow = 2, byrow = TRUE))
      cube <- build_tp13_test_cube()
      NF <- c(20, 25)
      NH <- c(50L, 60L)

      seasonal_a <- function(parameters, node_index, warmup_days) {
        parameters$model_seasonality <- TRUE
        parameters$g0 <- 2.1
        parameters$g <- c(0.25, 0.35, 0.15)
        parameters$h <- c(0.1, 0.5, 0.8)
        parameters
      }
      seasonal_b <- function(parameters, node_index, warmup_days) {
        parameters$model_seasonality <- TRUE
        parameters$g0 <- 2.4
        parameters$g <- c(0.25, 0.35, 0.15)
        parameters$h <- c(0.1, 0.5, 0.8)
        parameters
      }

      if (identical(case_name, "metadata")) {
        checkpoint <- suppressWarnings(msimGD_build_baseline_checkpoint(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 14L,
          parameter_modifier = seasonal_a,
          seed = 1L
        ))

        return(list(
          phase = checkpoint$metadata$seasonal_phase_day,
          cycle_days = checkpoint$metadata$seasonal_cycle_days,
          has_signature = !is.null(checkpoint$metadata$baseline_time_dependent_signature),
          node1_is_seasonal = checkpoint$metadata$baseline_time_dependent_signature$nodes[[1L]]$species[[1L]]$model_seasonality
        ))
      }

      if (identical(case_name, "mismatch")) {
        checkpoint <- suppressWarnings(msimGD_build_baseline_checkpoint(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 14L,
          parameter_modifier = seasonal_a,
          seed = 1L
        ))

        return(tryCatch(
          {
            suppressWarnings(msimGD_run_truth(
              setup = setup,
              cube = cube,
              NF = NF,
              NH = NH,
              tmax = 3L,
              mu = 1.0,
              p_move = 0.04,
              release = NULL,
              theta = theta_test,
              init_EIR = 5,
              warmup_days = 0L,
              parameter_modifier = seasonal_b,
              baseline_checkpoint = checkpoint,
              seed = 1L
            ))
            "ok"
          },
          error = function(e) conditionMessage(e)
        ))
      }

      if (identical(case_name, "legacy_static")) {
        checkpoint <- suppressWarnings(msimGD_build_baseline_checkpoint(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 14L,
          parameter_modifier = NULL,
          seed = 1L
        ))
        checkpoint$metadata$baseline_time_dependent_signature <- NULL
        checkpoint$metadata$seasonal_phase_day <- NULL
        checkpoint$metadata$seasonal_cycle_days <- NULL

        result <- suppressWarnings(msimGD_run_truth(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          tmax = 3L,
          mu = 1.0,
          p_move = 0.04,
          release = NULL,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 0L,
          parameter_modifier = NULL,
          baseline_checkpoint = checkpoint,
          seed = 1L
        ))

        return(list(
          n_nodes = length(result),
          rows = vapply(result, nrow, integer(1))
        ))
      }

      if (identical(case_name, "resolve_ochomo_seasonality")) {
        checkpoint <- suppressWarnings(msimGD_build_baseline_checkpoint(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 14L,
          parameter_modifier = seasonal_a,
          seed = 1L
        ))

        seasonality <- ochomo_resolve_baseline_checkpoint_seasonality(checkpoint)
        return(list(
          type = seasonality$type,
          g_len = length(seasonality$g),
          g0 = seasonality$g0
        ))
      }

      if (identical(case_name, "forward_phase_guard")) {
        checkpoint <- suppressWarnings(msimGD_build_baseline_checkpoint(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 14L,
          parameter_modifier = seasonal_a,
          seed = 1L
        ))

        return(ochomo_validate_forward_checkpoint_phase(checkpoint))
      }

      if (identical(case_name, "forward_phase_expected_mismatch")) {
        checkpoint <- suppressWarnings(msimGD_build_baseline_checkpoint(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 14L,
          parameter_modifier = seasonal_a,
          seed = 1L
        ))

        return(tryCatch(
          {
            ochomo_validate_forward_checkpoint_phase(
              checkpoint,
              expected_phase_day = 0L
            )
            "ok"
          },
          error = function(e) conditionMessage(e)
        ))
      }

      if (identical(case_name, "derive_study_start_checkpoint")) {
        checkpoint <- suppressWarnings(msimGD_build_baseline_checkpoint(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 14L,
          parameter_modifier = seasonal_a,
          seed = 1L
        ))
        checkpoint_path <- file.path(tempdir(), "accepted_checkpoint.rds")
        saveRDS(checkpoint, checkpoint_path)

        derived <- suppressWarnings(ochomo_resolve_study_start_checkpoint(
          accepted_checkpoint_path = checkpoint_path,
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          parameter_modifier = seasonal_a,
          study_start_offset_days = 9L,
          output_dir = tempdir()
        ))
        cached <- suppressWarnings(ochomo_resolve_study_start_checkpoint(
          accepted_checkpoint_path = checkpoint_path,
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          parameter_modifier = seasonal_a,
          study_start_offset_days = 9L,
          output_dir = tempdir()
        ))

        return(list(
          source = derived$source,
          cached_source = cached$source,
          parent_phase = checkpoint$metadata$seasonal_phase_day,
          phase = derived$checkpoint$metadata$seasonal_phase_day,
          offset = derived$checkpoint$metadata$study_start_offset_days,
          warmup_days = derived$checkpoint$metadata$warmup_days,
          accepted_checkpoint_path = basename(derived$checkpoint$metadata$accepted_checkpoint_path),
          same_signature = isTRUE(all.equal(
            derived$checkpoint$metadata$baseline_time_dependent_signature,
            checkpoint$metadata$baseline_time_dependent_signature,
            tolerance = 1e-8,
            check.attributes = FALSE
          )),
          same_contact_signature = isTRUE(all.equal(
            derived$checkpoint$metadata$baseline_contact_signature,
            checkpoint$metadata$baseline_contact_signature,
            tolerance = 1e-8,
            check.attributes = FALSE
          ))
        ))
      }

      if (identical(case_name, "derive_study_start_checkpoint_rng_key")) {
        checkpoint <- suppressWarnings(msimGD_build_baseline_checkpoint(
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          warmup_days = 14L,
          parameter_modifier = seasonal_a,
          seed = 1L
        ))
        checkpoint_path <- file.path(tempdir(), "accepted_checkpoint_rng.rds")
        output_dir <- file.path(tempdir(), "study_start_rng_key")
        unlink(output_dir, recursive = TRUE, force = TRUE)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(checkpoint, checkpoint_path)

        derived_preserve <- suppressWarnings(ochomo_resolve_study_start_checkpoint(
          accepted_checkpoint_path = checkpoint_path,
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          parameter_modifier = seasonal_a,
          study_start_offset_days = 9L,
          output_dir = output_dir,
          preserve_rng_stream = TRUE
        ))
        derived_fresh <- suppressWarnings(ochomo_resolve_study_start_checkpoint(
          accepted_checkpoint_path = checkpoint_path,
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          parameter_modifier = seasonal_a,
          study_start_offset_days = 9L,
          output_dir = output_dir,
          preserve_rng_stream = FALSE
        ))
        cached_fresh <- suppressWarnings(ochomo_resolve_study_start_checkpoint(
          accepted_checkpoint_path = checkpoint_path,
          setup = setup,
          cube = cube,
          NF = NF,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          parameter_modifier = seasonal_a,
          study_start_offset_days = 9L,
          output_dir = output_dir,
          preserve_rng_stream = FALSE
        ))

        return(list(
          preserve_path = basename(derived_preserve$checkpoint_path),
          fresh_path = basename(derived_fresh$checkpoint_path),
          preserve_source = derived_preserve$source,
          fresh_source = derived_fresh$source,
          cached_fresh_source = cached_fresh$source,
          preserve_flag_preserve = isTRUE(derived_preserve$checkpoint$metadata$preserve_rng_stream),
          preserve_flag_fresh = isTRUE(derived_fresh$checkpoint$metadata$preserve_rng_stream)
        ))
      }

      stop(sprintf("Unknown checkpoint test case: %s", case_name), call. = FALSE)
    },
    args = list(
      repo_root = repo_root,
      scripts_root = scripts_root,
      case_name = case_name
    )
  )
}


test_that("baseline checkpoints record seasonal phase and forcing signature", {
  metadata <- run_checkpoint_subprocess("metadata")

  expect_identical(metadata$phase, 14L)
  expect_identical(metadata$cycle_days, 365L)
  expect_true(metadata$has_signature)
  expect_true(metadata$node1_is_seasonal)
})


test_that("seasonal checkpoints reject mismatched forcing on resume", {
  message_text <- run_checkpoint_subprocess("mismatch")

  expect_match(
    message_text,
    "baseline_checkpoint baseline_time_dependent_signature does not match",
    fixed = TRUE
  )
})


test_that("legacy static checkpoints without forcing metadata remain usable", {
  result <- run_checkpoint_subprocess("legacy_static")

  expect_identical(result$n_nodes, 2L)
  expect_true(all(result$rows > 0L))
})


test_that("Ochomo release helpers reconstruct shared checkpoint seasonality", {
  resolved <- run_checkpoint_subprocess("resolve_ochomo_seasonality")

  expect_identical(resolved$type, "fourier_rainfall")
  expect_identical(resolved$g_len, 3L)
  expect_equal(resolved$g0, 2.1, tolerance = 1e-12)
})


test_that("forward checkpoint phase validator accepts non-zero seasonal phases", {
  phase_day <- run_checkpoint_subprocess("forward_phase_guard")

  expect_identical(phase_day, 14L)
})


test_that("forward checkpoint phase validator still enforces explicit mismatches", {
  message_text <- run_checkpoint_subprocess("forward_phase_expected_mismatch")

  expect_match(
    message_text,
    "Checkpoint seasonal_phase_day (14) does not match the expected seasonal_phase_day (0).",
    fixed = TRUE
  )
})


test_that("study-start checkpoint derivation preserves baseline signatures", {
  derived <- run_checkpoint_subprocess("derive_study_start_checkpoint")

  expect_identical(derived$source, "derived_checkpoint")
  expect_identical(derived$cached_source, "cached_derived_checkpoint")
  expect_identical(derived$parent_phase, 14L)
  expect_identical(derived$phase, 23L)
  expect_identical(derived$offset, 9L)
  expect_identical(derived$warmup_days, 23L)
  expect_identical(derived$accepted_checkpoint_path, "accepted_checkpoint.rds")
  expect_true(derived$same_signature)
  expect_true(derived$same_contact_signature)
})


test_that("study-start checkpoint cache keys include RNG continuation mode", {
  derived <- run_checkpoint_subprocess("derive_study_start_checkpoint_rng_key")

  expect_identical(derived$preserve_source, "derived_checkpoint")
  expect_identical(derived$fresh_source, "derived_checkpoint")
  expect_identical(derived$cached_fresh_source, "cached_derived_checkpoint")
  expect_false(identical(derived$preserve_path, derived$fresh_path))
  expect_match(derived$preserve_path, "preserve_rng", fixed = TRUE)
  expect_match(derived$fresh_path, "fresh_rng", fixed = TRUE)
  expect_true(derived$preserve_flag_preserve)
  expect_false(derived$preserve_flag_fresh)
})
