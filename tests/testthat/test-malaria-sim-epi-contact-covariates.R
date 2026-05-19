testthat::skip_if_not_installed("pkgload")
testthat::skip_if_not_installed("callr")
if (!requireNamespace("customMGDrive2", quietly = TRUE)) {
  testthat::skip("customMGDrive2 package is required for malaria-sim-epi contact-covariate tests.")
}

scripts_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
repo_root <- dirname(scripts_root)

required_full_repo_paths <- c(
  file.path(repo_root, "DESCRIPTION"),
  file.path(scripts_root, "landscape-movement", "movement_mu.R"),
  file.path(scripts_root, "archive", "env-codes", "cube.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_truth_generation.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_baseline_calibration_utils.R"),
  file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_covariate_utils.R")
)
if (!all(file.exists(required_full_repo_paths))) {
  testthat::skip(
    paste(
      "Full customMGDrive2 malaria-sim-epi integration fixtures are not",
      "bundled with this isolated package share."
    )
  )
}

load_test_busia_geometry <- function() {
  raw_lines <- readLines(
    file.path(scripts_root, "busia-landscape", "busia-study.csv"),
    warn = FALSE
  )
  csv_lines <- grep("^[0-9]+,", raw_lines, value = TRUE)
  csv_lines <- gsub("\\\\$", "", csv_lines)
  csv_lines <- gsub("\\}$", "", csv_lines)
  csv_lines <- trimws(csv_lines)
  raw <- read.csv(
    textConnection(paste(c("Name,Latitude,Longitude", csv_lines), collapse = "\n")),
    stringsAsFactors = FALSE
  )

  lat0 <- mean(raw$Latitude)
  lon0 <- mean(raw$Longitude)
  nodes <- data.frame(
    node = seq_len(nrow(raw)),
    x = (raw$Longitude - lon0) * 111.32 * cos(lat0 * pi / 180),
    y = (raw$Latitude - lat0) * 110.57,
    original_node = raw$Name,
    stringsAsFactors = FALSE
  )

  list(
    nodes = nodes,
    D = as.matrix(stats::dist(nodes[, c("x", "y"), drop = FALSE]))
  )
}


run_contact_subprocess <- function(case_name) {
  callr::r(
    function(repo_root, scripts_root, case_name) {
      library(pkgload)

      setwd(repo_root)
      pkgload::load_all(repo_root, quiet = TRUE, export_all = TRUE)
      pkgload::load_all(file.path(scripts_root, "msimGD"), quiet = TRUE, export_all = TRUE)
      source(file.path(scripts_root, "landscape-movement", "movement_mu.R"), local = TRUE)
      source(file.path(scripts_root, "archive", "env-codes", "cube.R"), local = TRUE)
      source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "msimGD_truth_generation.R"), local = TRUE)
      source(file.path(scripts_root, "msimGD", "test-scripts", "malaria-sim-epi", "ochomo_baseline_calibration_utils.R"), local = TRUE)

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

      surface_a <- resolve_node_contact_effect(
        node_data = data.frame(
          node_index = 1:2,
          open_eaves = c("no", "yes"),
          roof_type = c("metal", "tile"),
          stringsAsFactors = FALSE
        ),
        effect = node_contact_effect_bundle(
          effects = list(
            node_contact_effect_binary(
              covariate = "open_eaves",
              reference_level = "no",
              exposed_level = "yes",
              exposed_multiplier = 1.5,
              normalization = "none"
            ),
            node_contact_effect_categorical(
              covariate = "roof_type",
              level_multipliers = c(metal = 0.95, tile = 1.05),
              normalization = "none"
            )
          ),
          normalization = "mean_one",
          label = "surface_a",
          source = "test"
        )
      )
      surface_b <- resolve_node_contact_effect(
        node_data = data.frame(
          node_index = 1:2,
          open_eaves = c("yes", "no"),
          roof_type = c("metal", "tile"),
          stringsAsFactors = FALSE
        ),
        effect = node_contact_effect_bundle(
          effects = list(
            node_contact_effect_binary(
              covariate = "open_eaves",
              reference_level = "no",
              exposed_level = "yes",
              exposed_multiplier = 1.5,
              normalization = "none"
            ),
            node_contact_effect_categorical(
              covariate = "roof_type",
              level_multipliers = c(metal = 0.95, tile = 1.05),
              normalization = "none"
            )
          ),
          normalization = "mean_one",
          label = "surface_b",
          source = "test"
        )
      )
      modifier_a <- function(parameters, node_index, warmup_days) {
        apply_node_contact_effect(parameters, surface_a, node_index)
      }
      modifier_b <- function(parameters, node_index, warmup_days) {
        apply_node_contact_effect(parameters, surface_b, node_index)
      }

      if (identical(case_name, "checkpoint_metadata")) {
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
          parameter_modifier = modifier_a,
          seed = 1L
        ))

        return(list(
          has_contact_signature = !is.null(checkpoint$metadata$baseline_contact_signature),
          contact_hook = checkpoint$metadata$baseline_contact_signature$contact_hook,
          contact_multiplier_by_node = checkpoint$metadata$contact_multiplier_by_node
        ))
      }

      if (identical(case_name, "checkpoint_mismatch")) {
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
          parameter_modifier = modifier_a,
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
              parameter_modifier = modifier_b,
              baseline_checkpoint = checkpoint,
              seed = 1L
            ))
            "ok"
          },
          error = function(e) conditionMessage(e)
        ))
      }

      if (identical(case_name, "legacy_identity")) {
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
        checkpoint$metadata$baseline_contact_signature <- NULL
        checkpoint$metadata$contact_multiplier_by_node <- NULL

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

      if (identical(case_name, "human_library_mismatch")) {
        library <- suppressWarnings(msimGD_build_node_conditioned_human_library(
          setup = setup,
          cube = cube,
          NH = NH,
          mu = 1.0,
          p_move = 0.04,
          theta = theta_test,
          init_EIR = 5,
          parameter_modifier = modifier_a,
          burnin_timesteps = 30L,
          n_snapshots = 1L,
          snapshot_spacing = 0L,
          seed = 1L
        ))

        return(tryCatch(
          {
            suppressWarnings(msimGD_run_truth(
              setup = setup,
              cube = cube,
              NF = NULL,
              NH = NH,
              tmax = 3L,
              mu = 1.0,
              p_move = 0.04,
              release = NULL,
              theta = theta_test,
              init_EIR = 5,
              warmup_days = 0L,
              parameter_modifier = modifier_b,
              human_initialization_library = library,
              baseline_checkpoint = NULL,
              seed = 1L
            ))
            "ok"
          },
          error = function(e) conditionMessage(e)
        ))
      }

      if (identical(case_name, "checkpoint_contact_surface_resolve")) {
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
          parameter_modifier = modifier_a,
          seed = 1L
        ))

        resolved <- ochomo_resolve_baseline_checkpoint_contact_surface(checkpoint)
        return(list(
          type = resolved$type,
          hook = resolved$hook,
          node_index = resolved$node_index,
          contact_multiplier = unname(resolved$contact_multiplier),
          contact_multiplier_expected = unname(checkpoint$metadata$contact_multiplier_by_node)
        ))
      }

      stop(sprintf("Unknown contact test case: %s", case_name), call. = FALSE)
    },
    args = list(
      repo_root = repo_root,
      scripts_root = scripts_root,
      case_name = case_name
    )
  )
}


test_that("baseline checkpoints record contact-surface metadata", {
  metadata <- run_contact_subprocess("checkpoint_metadata")

  expect_true(metadata$has_contact_signature)
  expect_identical(metadata$contact_hook, "human_blood_meal_rate")
  expect_equal(mean(unname(metadata$contact_multiplier_by_node)), 1, tolerance = 1e-10)
  expect_true(all(unname(metadata$contact_multiplier_by_node) > 0))
})


test_that("contact-conditioned checkpoints reject mismatched surfaces on resume", {
  message_text <- run_contact_subprocess("checkpoint_mismatch")

  expect_match(
    message_text,
    "baseline_checkpoint baseline_contact_signature does not match",
    fixed = TRUE
  )
})


test_that("legacy checkpoints without contact metadata remain usable under identity contact", {
  result <- run_contact_subprocess("legacy_identity")

  expect_identical(result$n_nodes, 2L)
  expect_true(all(result$rows > 0L))
})


test_that("stationary human libraries reject mismatched contact surfaces", {
  message_text <- run_contact_subprocess("human_library_mismatch")

  expect_match(
    message_text,
    "human_initialization_library baseline_contact_signature does not match",
    fixed = TRUE
  )
})


test_that("Ochomo resolves a contact surface from checkpoint metadata", {
  result <- run_contact_subprocess("checkpoint_contact_surface_resolve")

  expect_identical(result$type, "contact_surface")
  expect_identical(result$hook, "human_blood_meal_rate")
  expect_equal(result$node_index, c(1L, 2L))
  expect_equal(
    result$contact_multiplier,
    result$contact_multiplier_expected,
    tolerance = 1e-12
  )
})


test_that("Ochomo prototype housing covariate generator is reproducible and practically centered", {
  source(
    file.path(
      scripts_root,
      "msimGD",
      "test-scripts",
      "malaria-sim-epi",
      "ochomo_covariate_utils.R"
    ),
    local = TRUE
  )

  busia_geometry <- load_test_busia_geometry()
  cov_a <- ochomo_generate_housing_contact_covariates(
    nodes = busia_geometry$nodes,
    D = busia_geometry$D,
    seed = 321L
  )
  cov_b <- ochomo_generate_housing_contact_covariates(
    nodes = busia_geometry$nodes,
    D = busia_geometry$D,
    seed = 321L
  )

  expect_equal(cov_a$node_covariates, cov_b$node_covariates, tolerance = 1e-12)
  expect_identical(nrow(cov_a$node_covariates), 58L)
  expect_true(all(cov_a$node_covariates$open_eaves_fraction > 0))
  expect_true(all(cov_a$node_covariates$open_eaves_fraction < 1))
  expect_true(all(cov_a$node_covariates$windows_mean > 0))
  expect_true(all(c(
    "housing_shared_latent",
    "housing_eaves_latent",
    "housing_windows_latent"
  ) %in% names(cov_a$node_covariates)))

  expect_lt(abs(mean(cov_a$node_covariates$open_eaves_fraction) - 0.6915), 0.03)
  expect_lt(abs(stats::sd(cov_a$node_covariates$open_eaves_fraction) - 0.10), 0.03)
  expect_lt(abs(mean(cov_a$node_covariates$windows_mean) - 2.32), 0.12)
  expect_lt(abs(stats::sd(cov_a$node_covariates$windows_mean) - 0.35), 0.10)
  expect_gt(stats::cor(
    cov_a$node_covariates$open_eaves_fraction,
    cov_a$node_covariates$windows_mean
  ), 0.2)
})


test_that("Ochomo housing contact translator builds one joint openness score", {
  source(
    file.path(
      scripts_root,
      "msimGD",
      "test-scripts",
      "malaria-sim-epi",
      "ochomo_covariate_utils.R"
    ),
    local = TRUE
  )

  covariates <- data.frame(
    node_index = 1:4,
    open_eaves_fraction = c(0.50, 0.60, 0.72, 0.82),
    windows_mean = c(1.4, 1.9, 2.6, 3.2),
    stringsAsFactors = FALSE
  )

  surface <- ochomo_housing_contact_surface_from_covariates(
    covariates = covariates,
    contact_multiplier_per_sd = 1.08
  )

  expect_equal(mean(unname(surface$contact_multiplier)), 1, tolerance = 1e-12)
  expect_true(all(unname(surface$contact_multiplier) > 0))
  expect_true(all(diff(surface$node_covariates$housing_openness_score) > 0))
  expect_true(all(diff(unname(surface$contact_multiplier)) > 0))
  expect_identical(
    surface$effect_spec$combination,
    "single_joint_housing_openness_score"
  )
  expect_match(
    surface$effect_spec$note,
    "not treated as mechanistically separable multiplicative effects",
    fixed = TRUE
  )
})


test_that("Ochomo housing_openness env mode returns resolved raw covariates and multiplier", {
  source(
    file.path(
      scripts_root,
      "msimGD",
      "test-scripts",
      "malaria-sim-epi",
      "ochomo_covariate_utils.R"
    ),
    local = TRUE
  )

  busia_small <- list(
    nodes = data.frame(
      node = 1:5,
      original_node = c(21L, 22L, 23L, 24L, 25L),
      x = c(0, 1, 2, 1.5, 0.5),
      y = c(0, 0.5, 1.0, 1.8, 2.2)
    ),
    D = as.matrix(stats::dist(rbind(
      c(0, 0),
      c(1, 0.5),
      c(2, 1.0),
      c(1.5, 1.8),
      c(0.5, 2.2)
    )))
  )

  old_env <- Sys.getenv(c(
    "MSIMGD_OCHOMO_CONTACT_SURFACE_MODE",
    "MSIMGD_OCHOMO_CONTACT_SURFACE_SEED",
    "MSIMGD_OCHOMO_CONTACT_SURFACE_RANGE",
    "MSIMGD_OCHOMO_CONTACT_SURFACE_SHARED_RHO",
    "MSIMGD_OCHOMO_EAVES_MEAN",
    "MSIMGD_OCHOMO_EAVES_SD_BETWEEN",
    "MSIMGD_OCHOMO_WINDOWS_MEAN",
    "MSIMGD_OCHOMO_WINDOWS_SD_BETWEEN",
    "MSIMGD_OCHOMO_CONTACT_MULTIPLIER_PER_SD",
    "MSIMGD_OCHOMO_EAVES_WEIGHT",
    "MSIMGD_OCHOMO_WINDOWS_WEIGHT"
  ), unset = NA_character_)
  on.exit({
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        Sys.setenv(structure(old_env[[nm]], names = nm))
      }
    }
  }, add = TRUE)

  Sys.setenv(
    MSIMGD_OCHOMO_CONTACT_SURFACE_MODE = "housing_openness",
    MSIMGD_OCHOMO_CONTACT_SURFACE_SEED = "777",
    MSIMGD_OCHOMO_CONTACT_SURFACE_RANGE = "4",
    MSIMGD_OCHOMO_CONTACT_SURFACE_SHARED_RHO = "0.6",
    MSIMGD_OCHOMO_EAVES_MEAN = "0.6915",
    MSIMGD_OCHOMO_EAVES_SD_BETWEEN = "0.10",
    MSIMGD_OCHOMO_WINDOWS_MEAN = "2.32",
    MSIMGD_OCHOMO_WINDOWS_SD_BETWEEN = "0.35",
    MSIMGD_OCHOMO_CONTACT_MULTIPLIER_PER_SD = "1.08",
    MSIMGD_OCHOMO_EAVES_WEIGHT = "0.5",
    MSIMGD_OCHOMO_WINDOWS_WEIGHT = "0.5"
  )

  surface <- ochomo_contact_surface_from_env(
    busia = busia_small,
    seed_base = 42L
  )

  expect_equal(mean(unname(surface$contact_multiplier)), 1, tolerance = 1e-12)
  expect_true(all(c(
    "open_eaves_fraction",
    "windows_mean",
    "housing_openness_score",
    "contact_multiplier"
  ) %in% names(surface$node_covariates)))
  expect_identical(
    surface$effect_spec$type,
    "ochomo_housing_openness_contact_effect"
  )
})
