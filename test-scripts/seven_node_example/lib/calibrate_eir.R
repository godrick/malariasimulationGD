# ------------------------------------------------------------------------------
# calibrate_eir.R
# ------------------------------------------------------------------------------
# Verbatim vendor of scripts/msimGD/test-scripts/malaria-sim-epi/calibrate_eir.R
# from the customMGDrive2 working repository, with two minor dependency
# adjustments at the top of the file (see below).
#
# Utilities for calibrating baseline msimGD transmission from observable
# equilibrium prevalence targets.

# ----- Vendored-file dependency stubs ---------------------------------------
# `EQUILIBRIUM_AGES` is an internal constant in the malariasimulationGD
# package (R/compatibility.R: `EQUILIBRIUM_AGES <- 0:999 / 10`); we redefine
# it here at file scope so the vendored functions below can reference it
# without `:::`.
if (!exists("EQUILIBRIUM_AGES", inherits = TRUE)) {
  EQUILIBRIUM_AGES <- 0:999 / 10
}
# ----------------------------------------------------------------------------
#
# Standard use:
#   1. Build the baseline parameter object you actually want to simulate.
#   2. Choose an equilibrium prevalence target (usually PfPR2-10, microscopy).
#   3. Solve for the init_EIR that yields that target under malariaEquilibrium.
#   4. Call set_equilibrium() with the calibrated init_EIR.
#
# This is the standard malariasimulationGD pattern described in the package
# vignette "Matching EIR to PfPR2-10". The key idea is:
#
#   target prevalence -> calibrated init_EIR -> set_equilibrium(init_EIR)
#
# Rather than choosing init_EIR arbitrarily.
#
# Important scope:
#   - These helpers are for baseline / equilibrium calibration.
#   - They are appropriate when the target is a steady-state microscopy
#     prevalence in a specified age band.
#   - They are not a substitute for full simulation-based calibration when
#     dynamic interventions, seasonality changes, or non-equilibrium targets are
#     the quantities of interest.
# ------------------------------------------------------------------------------


#' Build a baseline msimGD parameter object for equilibrium calibration.
#'
#' Prefer passing the actual parameter object you intend to simulate. If
#' `parameters` is `NULL`, a fresh object is created with `get_parameters()`.
#'
#' The optional `theta` argument mirrors the mosquito lifecycle overrides used by
#' the msimGD inference prototype so that the abundance implied by
#' `set_equilibrium()` is consistent with the mosquito biology used elsewhere in
#' the pipeline.
#'
#' @param parameters Existing msimGD parameter list. If supplied, it is used as
#'   the baseline object to calibrate.
#' @param NH Human population size used only when `parameters` is `NULL`.
#' @param theta Optional customMGDrive2 lifecycle parameter list. Supported keys:
#'   `qE`, `nE`, `qL`, `nL`, `qP`, `nP`, `muE`, `muL`, `muP`, `muF`, `beta`,
#'   and optionally `nu`.
#' @param overrides Optional named list of additional parameter overrides applied
#'   after the baseline object and `theta` overrides are set.
#' @param individual_mosquitoes Passed to `get_parameters()` only when
#'   `parameters` is `NULL`.
#' @param native_mosquito_backend Passed to `get_parameters()` only when
#'   `parameters` is `NULL`.
#' @return A parameter list suitable for `pfpr_from_eir()` and
#'   `calibrate_eir_from_pfpr()`.
msimGD_build_calibration_parameters <- function(
    parameters = NULL,
    NH = 770L,
    theta = NULL,
    overrides = NULL,
    individual_mosquitoes = TRUE,
    native_mosquito_backend = TRUE
) {
  if (is.null(parameters)) {
    parameters <- get_parameters(list(
      human_population = NH,
      individual_mosquitoes = individual_mosquitoes,
      native_mosquito_backend = native_mosquito_backend,
      progress_bar = FALSE
    ))
  }

  if (!is.null(theta)) {
    parameters$del <- theta$nE / theta$qE
    parameters$dl  <- theta$nL / theta$qL
    parameters$dpl <- theta$nP / theta$qP
    parameters$me  <- theta$muE
    parameters$ml  <- theta$muL
    parameters$mup <- theta$muP
    parameters$mum <- theta$muF
    parameters$beta <- theta$beta
    if (!is.null(theta$nE)) parameters$native_mosquito_nE <- theta$nE
    if (!is.null(theta$nL)) parameters$native_mosquito_nL <- theta$nL
    if (!is.null(theta$nP)) parameters$native_mosquito_nP <- theta$nP
    if (!is.null(theta$nu)) parameters$native_mosquito_nu <- theta$nu
  }

  if (!is.null(overrides)) {
    stopifnot(is.list(overrides), !is.null(names(overrides)))
    for (nm in names(overrides)) {
      parameters[[nm]] <- overrides[[nm]]
    }
  }

  parameters
}


.msimGD_equilibrium_treatment_coverage <- function(parameters, ft = NULL, timestep = 1L) {
  if (!is.null(ft)) {
    if (!is.numeric(ft) || length(ft) != 1L || !is.finite(ft) || ft < 0 || ft > 1) {
      stop("ft must be a single finite number in [0, 1].", call. = FALSE)
    }
    return(as.numeric(ft))
  }

  sum(malariasimulationGD:::get_treatment_coverages(parameters, timestep))
}


.msimGD_validate_age_band <- function(age_min, age_max) {
  if (!is.numeric(age_min) || length(age_min) != 1L || !is.finite(age_min) || age_min < 0) {
    stop("age_min must be a single finite number >= 0, in years.", call. = FALSE)
  }
  if (!is.numeric(age_max) || length(age_max) != 1L || !is.finite(age_max) || age_max <= age_min) {
    stop("age_max must be a single finite number > age_min, in years.", call. = FALSE)
  }
  invisible(TRUE)
}


.msimGD_validate_eir <- function(eir, name = "eir") {
  if (!is.numeric(eir) || length(eir) != 1L || !is.finite(eir) || eir < 0) {
    stop(sprintf("%s must be a single finite number >= 0.", name), call. = FALSE)
  }
  invisible(TRUE)
}


.msimGD_validate_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("%s must be a single finite number > 0.", name), call. = FALSE)
  }
  invisible(TRUE)
}


.msimGD_human_population_scalar <- function(parameters) {
  hp <- parameters$human_population
  if (!is.numeric(hp) || length(hp) != 1L || !is.finite(hp) || hp <= 0) {
    stop(
      "parameters$human_population must be a single finite number > 0 for baseline calibration.",
      call. = FALSE
    )
  }
  as.numeric(hp)
}


.msimGD_scale_contact_parameter <- function(
    parameters,
    contact_scale,
    contact_parameter = c("blood_meal_rates", "Q0")
) {
  contact_parameter <- match.arg(contact_parameter)
  .msimGD_validate_positive_scalar(contact_scale, "contact_scale")

  vals <- parameters[[contact_parameter]]
  if (is.null(vals)) {
    stop(sprintf("parameters does not contain `%s`.", contact_parameter), call. = FALSE)
  }

  scale_one <- function(x) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop(
        sprintf("`%s` must contain single finite numeric values.", contact_parameter),
        call. = FALSE
      )
    }
    as.numeric(x) * contact_scale
  }

  scaled <- if (is.list(vals)) {
    out <- lapply(vals, scale_one)
    names(out) <- names(vals)
    out
  } else if (is.numeric(vals)) {
    out <- vals * contact_scale
    names(out) <- names(vals)
    out
  } else {
    stop(
      sprintf("`%s` must be numeric or a list of numeric scalars.", contact_parameter),
      call. = FALSE
    )
  }

  scaled_vec <- if (is.list(scaled)) unlist(scaled, use.names = FALSE) else as.numeric(scaled)
  if (contact_parameter == "Q0") {
    if (any(scaled_vec < 0 | scaled_vec > 1)) {
      stop(
        sprintf(
          "Scaling `%s` by %.6g produced values outside [0, 1].",
          contact_parameter, contact_scale
        ),
        call. = FALSE
      )
    }
  } else {
    if (any(scaled_vec <= 0)) {
      stop(
        sprintf(
          "Scaling `%s` by %.6g produced non-positive values.",
          contact_parameter, contact_scale
        ),
        call. = FALSE
      )
    }
  }

  parameters[[contact_parameter]] <- scaled
  parameters
}


#' Compute equilibrium microscopy prevalence from EIR.
#'
#' This uses `malariaEquilibrium::human_equilibrium()` rather than running a full
#' msimGD simulation. Ages are in years because `EQUILIBRIUM_AGES` in msimGD are
#' defined on a yearly scale.
#'
#' @param eir Equilibrium EIR in infectious bites per person per year.
#' @param parameters msimGD parameter list representing the intended baseline
#'   setting.
#' @param ft Optional effective treatment coverage in `[0, 1]`. If `NULL`,
#'   the coverage is derived from `parameters` using `malariasimulationGD:::get_treatment_coverages()`.
#' @param age_min Lower bound of the target age band in years. Default `2`.
#' @param age_max Upper bound of the target age band in years. Default `10`.
#' @param prevalence_col Column of the equilibrium state table to aggregate.
#'   Default `"pos_M"` which corresponds to microscopy-detectable prevalence and
#'   therefore matches PfPR.
#' @param ages Age support for the analytical equilibrium. Defaults to
#'   `EQUILIBRIUM_AGES`.
#' @param quadrature Optional heterogeneity quadrature object. If `NULL`, uses
#'   `malariaEquilibrium::gq_normal(parameters$n_heterogeneity_groups)`.
#' @return A scalar equilibrium prevalence in the requested age band.
pfpr_from_eir <- function(
    eir,
    parameters,
    ft = NULL,
    age_min = 2,
    age_max = 10,
    prevalence_col = "pos_M",
    ages = EQUILIBRIUM_AGES,
    quadrature = NULL
) {
  .msimGD_validate_eir(eir)
  .msimGD_validate_age_band(age_min, age_max)

  if (is.null(parameters) || !is.list(parameters)) {
    stop("parameters must be an msimGD parameter list.", call. = FALSE)
  }

  # translate_parameters is an internal helper in malariasimulationGD; accessed
  # via ::: because it is not exported.
  eq_params <- malariasimulationGD:::translate_parameters(parameters)
  if (is.null(quadrature)) {
    quadrature <- malariaEquilibrium::gq_normal(parameters$n_heterogeneity_groups)
  }

  ft_eff <- .msimGD_equilibrium_treatment_coverage(parameters, ft = ft, timestep = 1L)

  eq <- malariaEquilibrium::human_equilibrium(
    EIR = eir,
    ft = ft_eff,
    p = eq_params,
    age = ages,
    h = quadrature
  )

  states <- eq$states
  state_cols <- colnames(states)
  if (is.null(state_cols)) {
    stop("human_equilibrium() output did not contain named state columns.", call. = FALSE)
  }
  if (!("age" %in% state_cols)) {
    stop("human_equilibrium() output did not contain an `age` column.", call. = FALSE)
  }
  if (!(prevalence_col %in% state_cols)) {
    stop(sprintf("Column `%s` was not found in equilibrium state output.", prevalence_col), call. = FALSE)
  }
  if (!("prop" %in% state_cols)) {
    stop("human_equilibrium() output did not contain a `prop` column.", call. = FALSE)
  }

  age_idx <- which(states[, "age"] >= age_min & states[, "age"] <= age_max)
  if (length(age_idx) == 0L) {
    stop("No equilibrium age bins fell inside the requested age band.", call. = FALSE)
  }

  denom <- sum(states[age_idx, "prop"])
  if (!is.finite(denom) || denom <= 0) {
    stop("Equilibrium age-band denominator was non-positive.", call. = FALSE)
  }

  sum(states[age_idx, prevalence_col]) / denom
}


#' Evaluate the equilibrium EIR-to-prevalence relationship on a grid.
#'
#' This is useful for diagnosing the calibration problem before solving it:
#' checking monotonicity, feasible prevalence ranges, and sensible search bounds.
#'
#' @param eir_values Numeric vector of candidate EIR values.
#' @param parameters msimGD parameter list representing the intended baseline
#'   setting.
#' @param ft Optional effective treatment coverage. If `NULL`, derived from
#'   `parameters`.
#' @param age_min Lower age bound in years.
#' @param age_max Upper age bound in years.
#' @param prevalence_col Equilibrium state column to aggregate. Default `"pos_M"`.
#' @return A data.frame with columns `init_EIR` and `prevalence`.
pfpr_curve_from_eir <- function(
    eir_values,
    parameters,
    ft = NULL,
    age_min = 2,
    age_max = 10,
    prevalence_col = "pos_M"
) {
  if (!is.numeric(eir_values) || length(eir_values) < 1L || any(!is.finite(eir_values)) || any(eir_values < 0)) {
    stop("eir_values must be a numeric vector of finite values >= 0.", call. = FALSE)
  }

  data.frame(
    init_EIR = as.numeric(eir_values),
    prevalence = vapply(
      eir_values,
      function(eir) {
        pfpr_from_eir(
          eir = eir,
          parameters = parameters,
          ft = ft,
          age_min = age_min,
          age_max = age_max,
          prevalence_col = prevalence_col
        )
      },
      numeric(1)
    )
  )
}


#' Calibrate baseline init_EIR from a target equilibrium prevalence.
#'
#' This function inverts the equilibrium prevalence map using `uniroot()` and
#' then calls `set_equilibrium()` to recover the matching `init_foim` and
#' mosquito abundance.
#'
#' Recommended workflow:
#'   1. Build the baseline parameter list you actually intend to simulate
#'      using `msimGD_build_calibration_parameters()` or by passing your own
#'      parameter object directly.
#'   2. Choose a target baseline prevalence, usually PfPR2-10.
#'   3. Calibrate `init_EIR` with this function.
#'   4. Use the returned `parameters_eq` or the returned `init_EIR` in your
#'      truth/forward pipeline.
#'
#' Notes:
#'   - `target_pfpr` is interpreted as equilibrium microscopy prevalence in the
#'     requested age band when `prevalence_col = "pos_M"`.
#'   - Age bounds are in years.
#'   - `theta` overrides mainly matter for the recovered `total_M`; they usually
#'     do not strongly change the equilibrium PfPR-to-EIR map unless they alter
#'     parameters used by `translate_parameters()`.
#'
#' @param target_pfpr Target equilibrium prevalence in `[0, 1]`.
#' @param parameters Optional existing msimGD parameter list. If `NULL`, one is
#'   built using `NH`, `theta`, and `overrides`.
#' @param NH Human population size used only when `parameters` is `NULL`.
#' @param theta Optional lifecycle parameter list used to override mosquito
#'   biology in the baseline parameter object.
#' @param overrides Optional named list of extra parameter overrides applied when
#'   building the baseline parameter object.
#' @param ft Optional effective treatment coverage in `[0, 1]`. If `NULL`,
#'   treatment coverage is derived from `parameters`.
#' @param age_min Lower age bound in years. Default `2`.
#' @param age_max Upper age bound in years. Default `10`.
#' @param prevalence_col Equilibrium state column to match. Default `"pos_M"`.
#' @param eir_range Search interval for `uniroot()`, in infectious bites per
#'   person per year. Default `c(0.01, 500)`.
#' @param tol Root-finding tolerance.
#' @param return_parameters If `TRUE`, include the fully equilibrated parameter
#'   list as `parameters_eq` in the return value.
#' @return A list containing the calibrated `init_EIR`, the achieved prevalence,
#'   derived `total_M`, derived `init_foim`, treatment coverage used, age band,
#'   and optionally the equilibrated parameter list.
calibrate_eir_from_pfpr <- function(
    target_pfpr,
    parameters = NULL,
    NH = 770L,
    theta = NULL,
    overrides = NULL,
    ft = NULL,
    age_min = 2,
    age_max = 10,
    prevalence_col = "pos_M",
    eir_range = c(0.01, 500),
    tol = 1e-6,
    return_parameters = FALSE
) {
  if (!is.numeric(target_pfpr) || length(target_pfpr) != 1L ||
      !is.finite(target_pfpr) || target_pfpr < 0 || target_pfpr > 1) {
    stop("target_pfpr must be a single finite number in [0, 1].", call. = FALSE)
  }
  .msimGD_validate_age_band(age_min, age_max)

  if (!is.numeric(eir_range) || length(eir_range) != 2L || any(!is.finite(eir_range)) ||
      eir_range[[1]] < 0 || eir_range[[2]] <= eir_range[[1]]) {
    stop("eir_range must be a numeric vector c(lower, upper) with 0 <= lower < upper.", call. = FALSE)
  }

  baseline_params <- msimGD_build_calibration_parameters(
    parameters = parameters,
    NH = NH,
    theta = theta,
    overrides = overrides
  )

  ft_eff <- .msimGD_equilibrium_treatment_coverage(baseline_params, ft = ft, timestep = 1L)

  pfpr_target_fn <- function(eir) {
    pfpr_from_eir(
      eir = eir,
      parameters = baseline_params,
      ft = ft_eff,
      age_min = age_min,
      age_max = age_max,
      prevalence_col = prevalence_col
    ) - target_pfpr
  }

  f_lo <- pfpr_target_fn(eir_range[[1]])
  f_hi <- pfpr_target_fn(eir_range[[2]])

  if (isTRUE(all.equal(f_lo, 0, tolerance = tol))) {
    root <- eir_range[[1]]
  } else if (isTRUE(all.equal(f_hi, 0, tolerance = tol))) {
    root <- eir_range[[2]]
  } else if (sign(f_lo) == sign(f_hi)) {
    curve_diag <- pfpr_curve_from_eir(
      eir_values = seq(eir_range[[1]], eir_range[[2]], length.out = 50L),
      parameters = baseline_params,
      ft = ft_eff,
      age_min = age_min,
      age_max = age_max,
      prevalence_col = prevalence_col
    )
    stop(
      sprintf(
        paste(
          "target_pfpr %.4f was not bracketed by eir_range [%.4f, %.4f].",
          "Over this range, equilibrium prevalence spans approximately [%.4f, %.4f].",
          "Expand eir_range or revise the baseline parameter set."
        ),
        target_pfpr,
        eir_range[[1]],
        eir_range[[2]],
        min(curve_diag$prevalence),
        max(curve_diag$prevalence)
      ),
      call. = FALSE
    )
  } else {
    root <- uniroot(pfpr_target_fn, interval = eir_range, tol = tol)$root
  }

  parameters_eq <- set_equilibrium(baseline_params, init_EIR = root)
  achieved_prev <- pfpr_from_eir(
    eir = root,
    parameters = baseline_params,
    ft = ft_eff,
    age_min = age_min,
    age_max = age_max,
    prevalence_col = prevalence_col
  )

  out <- list(
    init_EIR = unname(as.numeric(root)),
    prevalence = unname(as.numeric(achieved_prev)),
    target_prevalence = unname(as.numeric(target_pfpr)),
    prevalence_col = prevalence_col,
    age_min_years = age_min,
    age_max_years = age_max,
    ft = ft_eff,
    total_M = unname(as.numeric(parameters_eq$total_M)),
    init_foim = unname(as.numeric(parameters_eq$init_foim))
  )

  if (isTRUE(return_parameters)) {
    out$parameters_eq <- parameters_eq
  }

  out
}


#' Jointly calibrate baseline prevalence and adult females per human.
#'
#' This is a two-target equilibrium calibration. It is useful when you want the
#' baseline epidemiology to match an observed equilibrium prevalence, while also
#' constraining the baseline mosquito abundance to a concrete entomological
#' target expressed as adult females per human.
#'
#' The helper solves a nested root-finding problem:
#'   1. Outer root: find a uniform `contact_scale` applied to a chosen mosquito
#'      contact parameter (`blood_meal_rates` by default, or `Q0`).
#'   2. Inner root: for that scaled contact process, solve for the `init_EIR`
#'      that matches the requested equilibrium prevalence.
#'   3. Evaluate the implied equilibrium `total_M`; the outer root stops when it
#'      matches `target_females_per_human * human_population`.
#'
#' This is more principled than manually overriding `total_M` after prevalence
#' calibration because the final baseline remains internally consistent under a
#' modified contact process.
#'
#' Important interpretation:
#'   - This helper does not give you "more mosquitoes with everything else
#'     unchanged". It changes both mosquito abundance and the effective contact
#'     process so that the requested prevalence and the requested females-per-
#'     human target can coexist at equilibrium.
#'   - The returned `contact_scale` is therefore a model change, not a cosmetic
#'     multiplier.
#'   - `contact_parameter = "blood_meal_rates"` is the safest default because it
#'     is positive-valued. `contact_parameter = "Q0"` is supported, but the
#'     scaled values must remain in `[0, 1]`.
#'
#' Recommended workflow:
#'   1. Choose a target equilibrium prevalence, usually PfPR2-10.
#'   2. Choose a target baseline adult-female density in females per human.
#'   3. Call this function to recover the joint baseline calibration.
#'   4. If desired, use the returned `parameters_eq` directly in later work.
#'
#' @param target_pfpr Target equilibrium prevalence in `[0, 1]`.
#' @param target_females_per_human Target adult female mosquitoes per human at
#'   equilibrium. This is converted to `target_total_M` using
#'   `parameters$human_population`.
#' @param parameters Optional existing msimGD parameter list. If `NULL`, one is
#'   built using `NH`, `theta`, and `overrides`.
#' @param NH Human population size used only when `parameters` is `NULL`.
#' @param theta Optional lifecycle parameter list used to override mosquito
#'   biology in the baseline parameter object.
#' @param overrides Optional named list of extra parameter overrides applied when
#'   building the baseline parameter object.
#' @param ft Optional effective treatment coverage in `[0, 1]`. If `NULL`,
#'   treatment coverage is derived from `parameters`.
#' @param age_min Lower age bound in years. Default `2`.
#' @param age_max Upper age bound in years. Default `10`.
#' @param prevalence_col Equilibrium state column to match. Default `"pos_M"`.
#' @param contact_parameter Which contact term to scale uniformly across species.
#'   Default `"blood_meal_rates"`. Alternative `"Q0"`.
#' @param eir_range Search interval for the inner `init_EIR` solve, in
#'   infectious bites per person per year. Default `c(0.01, 500)`.
#' @param contact_scale_range Search interval for the outer contact-scale root.
#'   Default `c(0.05, 20)`.
#' @param eir_tol Inner root-finding tolerance for `init_EIR`.
#' @param contact_tol Outer root-finding tolerance for `contact_scale`.
#' @param return_parameters If `TRUE`, include the fully equilibrated parameter
#'   list as `parameters_eq` in the return value.
#' @return A list containing the calibrated `init_EIR`, achieved prevalence,
#'   target females per human, implied `total_M`, calibrated `contact_scale`,
#'   chosen contact parameter, and optionally the equilibrated parameter list.
calibrate_baseline_from_pfpr_and_fph <- function(
    target_pfpr,
    target_females_per_human,
    parameters = NULL,
    NH = 770L,
    theta = NULL,
    overrides = NULL,
    ft = NULL,
    age_min = 2,
    age_max = 10,
    prevalence_col = "pos_M",
    contact_parameter = c("blood_meal_rates", "Q0"),
    eir_range = c(0.01, 500),
    contact_scale_range = c(0.05, 20),
    eir_tol = 1e-6,
    contact_tol = 1e-6,
    return_parameters = FALSE
) {
  contact_parameter <- match.arg(contact_parameter)

  if (!is.numeric(target_pfpr) || length(target_pfpr) != 1L ||
      !is.finite(target_pfpr) || target_pfpr < 0 || target_pfpr > 1) {
    stop("target_pfpr must be a single finite number in [0, 1].", call. = FALSE)
  }
  .msimGD_validate_positive_scalar(target_females_per_human, "target_females_per_human")
  .msimGD_validate_age_band(age_min, age_max)

  if (!is.numeric(contact_scale_range) || length(contact_scale_range) != 2L ||
      any(!is.finite(contact_scale_range)) || contact_scale_range[[1]] <= 0 ||
      contact_scale_range[[2]] <= contact_scale_range[[1]]) {
    stop(
      paste(
        "contact_scale_range must be a numeric vector c(lower, upper)",
        "with 0 < lower < upper."
      ),
      call. = FALSE
    )
  }

  baseline_params <- msimGD_build_calibration_parameters(
    parameters = parameters,
    NH = NH,
    theta = theta,
    overrides = overrides
  )

  human_population <- .msimGD_human_population_scalar(baseline_params)
  target_total_M <- target_females_per_human * human_population

  eval_cache <- new.env(parent = emptyenv())

  evaluate_scale <- function(contact_scale, return_parameters_local = FALSE) {
    key <- paste(
      formatC(contact_scale, digits = 17L, format = "fg", flag = "#"),
      return_parameters_local,
      sep = "|"
    )
    if (exists(key, envir = eval_cache, inherits = FALSE)) {
      return(get(key, envir = eval_cache, inherits = FALSE))
    }

    scaled_params <- .msimGD_scale_contact_parameter(
      parameters = baseline_params,
      contact_scale = contact_scale,
      contact_parameter = contact_parameter
    )

    cal <- calibrate_eir_from_pfpr(
      target_pfpr = target_pfpr,
      parameters = scaled_params,
      ft = ft,
      age_min = age_min,
      age_max = age_max,
      prevalence_col = prevalence_col,
      eir_range = eir_range,
      tol = eir_tol,
      return_parameters = return_parameters_local
    )

    cal$contact_parameter <- contact_parameter
    cal$contact_scale <- unname(as.numeric(contact_scale))
    cal$target_females_per_human <- unname(as.numeric(target_females_per_human))
    cal$target_total_M <- unname(as.numeric(target_total_M))
    cal$human_population <- human_population
    cal$females_per_human <- unname(as.numeric(cal$total_M / human_population))

    assign(key, cal, envir = eval_cache)
    cal
  }

  scale_target_fn <- function(contact_scale) {
    evaluate_scale(contact_scale, return_parameters_local = FALSE)$total_M - target_total_M
  }

  f_lo <- scale_target_fn(contact_scale_range[[1]])
  f_hi <- scale_target_fn(contact_scale_range[[2]])

  if (isTRUE(all.equal(f_lo, 0, tolerance = contact_tol))) {
    root <- contact_scale_range[[1]]
  } else if (isTRUE(all.equal(f_hi, 0, tolerance = contact_tol))) {
    root <- contact_scale_range[[2]]
  } else if (sign(f_lo) == sign(f_hi)) {
    diag_scales <- exp(seq(
      log(contact_scale_range[[1]]),
      log(contact_scale_range[[2]]),
      length.out = 25L
    ))
    diag_vals <- vapply(
      diag_scales,
      function(s) evaluate_scale(s, return_parameters_local = FALSE)$total_M,
      numeric(1)
    )
    stop(
      sprintf(
        paste(
          "target_total_M %.4f was not bracketed by contact_scale_range [%.4f, %.4f].",
          "Over this range, the implied equilibrium total_M spans approximately [%.4f, %.4f].",
          "Expand contact_scale_range or revise the target prevalence / females-per-human pair."
        ),
        target_total_M,
        contact_scale_range[[1]],
        contact_scale_range[[2]],
        min(diag_vals),
        max(diag_vals)
      ),
      call. = FALSE
    )
  } else {
    root <- uniroot(scale_target_fn, interval = contact_scale_range, tol = contact_tol)$root
  }

  out <- evaluate_scale(root, return_parameters_local = return_parameters)
  baseline_ref <- evaluate_scale(1, return_parameters_local = FALSE)

  out$baseline_contact_scale <- 1
  out$baseline_init_EIR <- baseline_ref$init_EIR
  out$baseline_total_M <- baseline_ref$total_M
  out$baseline_females_per_human <- baseline_ref$females_per_human

  out
}
