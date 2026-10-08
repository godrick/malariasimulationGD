#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Deterministic-mosquito screen for the stationary-ITN Busia calibration.
#
# Human state remains individual-based. Only the native mosquito backend uses
# count-based deterministic dynamics, so results locate plausible (delta,
# init_EIR) regions and must be confirmed by the stochastic calibration.
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/Busia_landscape/careful_normalized_contact/00_calibrate_init_eir_incidence_normalized_stationary_itn_deterministic_screen.R",
    mustWork = TRUE
  )
}
careful_dir <- dirname(script_path)
example_dir <- normalizePath(file.path(careful_dir, ".."), mustWork = TRUE)

Sys.setenv(
  BUSIA_INDIVIDUAL_MOSQUITOES = "false",
  BUSIA_STABLE_NORMALIZED_CAREFUL_OUT_DIR = file.path(
    example_dir,
    "output",
    "careful_normalized_contact_stationary_itn_deterministic_screen"
  )
)

source(
  file.path(careful_dir,
            "00_calibrate_init_eir_incidence_normalized_stationary_itn_coarse.R"),
  local = globalenv()
)
