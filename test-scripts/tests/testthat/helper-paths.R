# Resolve paths from this test directory, independently of the caller's cwd.
developer_test_repo <- function() {
  normalizePath(file.path(testthat::test_path(), "../../.."), mustWork = TRUE)
}

developer_integration_scripts_root <- function() {
  normalizePath(Sys.getenv(
    "MSIMGD_INTEGRATION_SCRIPTS_ROOT",
    unset = dirname(developer_test_repo())
  ), mustWork = TRUE)
}

developer_integration_repo_root <- function() {
  normalizePath(Sys.getenv(
    "MSIMGD_INTEGRATION_REPO_ROOT",
    unset = dirname(developer_integration_scripts_root())
  ), mustWork = TRUE)
}
