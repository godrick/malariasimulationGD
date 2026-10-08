#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Coarse Busia calibration with a Godrick/P20-style stationary ITN background.
#
# This is a screening run: it searches broadly for the approximate joint
# (delta, raw counterfactual init_EIR) region under fixed ITN use before a
# narrower, longer final calibration is run. Apart from ITN timing/decay and
# these search controls, it uses the careful_normalized_contact settings.
# ------------------------------------------------------------------------------

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "test-scripts/Busia_landscape/careful_normalized_contact/00_calibrate_init_eir_incidence_normalized_stationary_itn_coarse.R",
    mustWork = TRUE
  )
}
careful_dir <- dirname(script_path)
example_dir <- normalizePath(file.path(careful_dir, ".."), mustWork = TRUE)

set_default_env <- function(name, value) {
  current <- Sys.getenv(name, unset = "")
  if (!nzchar(current)) {
    do.call(Sys.setenv, stats::setNames(list(as.character(value)), name))
  }
}

# This file's defining scientific setting: initialise 68.4% active ITN use in
# place, skip the nominal distribution, and use near-infinite retention and
# efficacy half-life. Do not allow an ambient campaign setting to override it.
Sys.setenv(
  BUSIA_ITN_REGIME = "stationary_background",
  BUSIA_ITN_DN0 = "0.437",
  BUSIA_ITN_RN = "0.555",
  BUSIA_ITN_RNM = "0.24"
)

# Broad screening defaults. Each can be overridden at the command line through
# its environment variable when a different search range is warranted.
set_default_env("BUSIA_DELTA_GRID", "0.22,0.25,0.28,0.31")
set_default_env("BUSIA_STABLE_USE_EIR_ANCHORS", "false")
set_default_env("BUSIA_STABLE_REUSE_NEAREST_EIR", "false")
set_default_env("BUSIA_STABLE_EIR_START", "85")
set_default_env("BUSIA_STABLE_EIR_MIN", "70")
set_default_env("BUSIA_STABLE_EIR_MAX", "110")
set_default_env("BUSIA_EIR_MAX_BRACKET_ITER", "1")
set_default_env("BUSIA_CALIBRATION_MAX_REFINE", "1")
set_default_env("BUSIA_DELTA_REFINE_STEP", "0")

set_default_env("BUSIA_STABLE_PARALLEL_INITIAL_DELTAS", "true")
set_default_env("BUSIA_STABLE_PARALLEL_WORKERS", "4")

# Six simulated years permits two final annual windows while avoiding the
# 22-year diagnostic horizon used by the campaign-cycle stability workflow.
set_default_env("BUSIA_STABILITY_DIAGNOSTIC_WARMUP_DAYS", 10L * 365L)
set_default_env("BUSIA_STABILITY_YEARS", "2")
set_default_env(
  "BUSIA_STABLE_NORMALIZED_CAREFUL_OUT_DIR",
  file.path(example_dir, "output", "careful_normalized_contact_stationary_itn_coarse")
)

source(
  file.path(careful_dir, "00_calibrate_init_eir_incidence_normalized_stable.R"),
  local = globalenv()
)
