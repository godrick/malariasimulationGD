testthat::skip_if_not_installed("pkgload")

scripts_root <- normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
repo_root <- dirname(scripts_root)

required_full_repo_paths <- c(
  file.path(repo_root, "DESCRIPTION"),
  file.path(scripts_root, "msimGD", "DESCRIPTION"),
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

pkgload::load_all(repo_root, quiet = TRUE, export_all = TRUE)
pkgload::load_all(file.path(scripts_root, "msimGD"), quiet = TRUE, export_all = TRUE)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)
source(
  file.path(
    scripts_root,
    "msimGD",
    "test-scripts",
    "malaria-sim-epi",
    "ochomo_baseline_calibration_utils.R"
  ),
  local = TRUE
)

test_that("ochomo_resolve_baseline_total_M prefers checkpoint parameter surfaces", {
  checkpoint_path <- tempfile(fileext = ".rds")
  on.exit(unlink(checkpoint_path), add = TRUE)

  saveRDS(
    list(metadata = list(
      parameter_total_M = c(10, 20),
      baseline_total_M_by_node = c(100, 200),
      effective_total_M = c(11, 21)
    )),
    checkpoint_path
  )

  baseline_fit <- list(
    fitted_baseline = list(
      checkpoint_parameter_total_M = c(12, 22),
      NF_eq = c(300, 400),
      baseline_total_M_by_node = c(500, 600)
    )
  )

  resolved <- ochomo_resolve_baseline_total_M(
    checkpoint_path = checkpoint_path,
    baseline_fit = baseline_fit,
    NH = c(100, 100),
    theta = list(),
    init_EIR = 1
  )

  expect_equal(resolved$total_M, c(10, 20))
  expect_identical(resolved$source, "checkpoint_parameter_surface")
})


test_that("ochomo_resolve_baseline_total_M falls back to fit checkpoint surfaces", {
  checkpoint_path <- tempfile(fileext = ".rds")
  on.exit(unlink(checkpoint_path), add = TRUE)

  saveRDS(
    list(metadata = list(
      baseline_total_M_by_node = c(100, 200)
    )),
    checkpoint_path
  )

  baseline_fit <- list(
    fitted_baseline = list(
      checkpoint_parameter_total_M = c(12, 22),
      NF_eq = c(300, 400),
      baseline_total_M_by_node = c(500, 600)
    )
  )

  resolved <- ochomo_resolve_baseline_total_M(
    checkpoint_path = checkpoint_path,
    baseline_fit = baseline_fit,
    NH = c(100, 100),
    theta = list(),
    init_EIR = 1
  )

  expect_equal(resolved$total_M, c(12, 22))
  expect_identical(resolved$source, "baseline_fit_checkpoint_parameter_surface")
})
