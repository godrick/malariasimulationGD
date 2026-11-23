# Deterministic, mosquito-only comparison of installed MGDrivE2 and the native
# malariasimulationGD entomology engine. No replacement lifecycle equations.

mc_native <- function(name) get(name, asNamespace("malariasimulationGD"))

mc_require <- function() {
  for (package in c("malariasimulationGD", "MGDrivE2", "MGDrivE", "deSolve")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Comparison requires installed package '", package, "'.", call. = FALSE)
    }
  }
}

mc_theta <- function() {
  list(qE = 1 / 4, nE = 2L, qL = 1 / 3, nL = 3L,
       qP = 1 / 6, nP = 2L, muE = 0.05, muL = 0.15,
       muP = 0.05, muF = 0.09, muM = 0.09, beta = 16, nu = 6)
}

mc_parameters <- function(theta, cube, rtol, atol) {
  if (theta$muM != theta$muF) {
    stop("The standalone native interface requires equal baseline sex mortality.", call. = FALSE)
  }
  malariasimulationGD::get_mosquito_parameters(list(
    total_M = 500, native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    native_mosquito_nE = theta$nE, native_mosquito_nL = theta$nL,
    native_mosquito_nP = theta$nP, native_mosquito_nu = theta$nu,
    # q is inverse mean duration, NOT the progression rate per Erlang substage.
    del = 1 / theta$qE, dl = 1 / theta$qL, dpl = 1 / theta$qP,
    me = theta$muE, ml = theta$muL, mup = theta$muP,
    mum = theta$muF, beta = theta$beta, model_seasonality = FALSE,
    r_tol = rtol, a_tol = atol, ode_max_steps = 100000L
  ), cube = cube)
}

mc_release <- function(days, stage, genotype = "aa", count = 100,
                       substage = NA_integer_, mate = NA_character_) {
  n <- length(days)
  data.frame(timestep = as.integer(days), node = rep("node1", n), species = rep("gamb", n),
             stage = rep_len(stage, n), genotype = rep_len(genotype, n), count = rep_len(count, n),
             substage = rep_len(substage, n), mate_genotype = rep_len(mate, n),
             stringsAsFactors = FALSE)
}

mc_scenarios <- function() {
  base <- MGDrivE::cubeMendelian()
  sterile <- base
  # A female inseminated by an aa male produces no viable offspring.
  sterile$tau[, "aa", ] <- 0
  costly <- base
  costly$omega["aa"] <- 0.8
  costly$s["aa"] <- 0.6
  costly$xiF["aa"] <- 0.7
  costly$xiM["aa"] <- 0.9
  costly$eta["AA", "aa"] <- 0.6
  costly$eta["aa", "AA"] <- 0.8
  homing <- base
  # Illustrative symmetric homing: Aa transmits a with probability 0.95.
  gamete <- c(AA = 0, Aa = 0.95, aa = 1)
  for (gf in names(gamete)) for (gm in names(gamete)) {
    pf <- gamete[[gf]]
    pm <- gamete[[gm]]
    homing$ih[gf, gm, ] <- c((1 - pf) * (1 - pm),
                            pf * (1 - pm) + (1 - pf) * pm, pf * pm)
  }
  scenario <- function(id, label, group, releases = NULL, intervention = "none",
                       cube = base, initial = "equilibrium", expected = "overlap",
                       native_beta_multiplier = 1) {
    list(id = id, label = label, group = group, releases = releases,
         intervention = intervention, cube = cube, initial = initial,
         expected = expected, native_beta_multiplier = native_beta_multiplier)
  }
  repeated <- mc_release(seq.int(20L, 60L, 10L), "adult_male", count = 250)
  list(
    scenario("baseline", "Unperturbed equilibrium", "Baseline"),
    scenario("adult_depletion", "70% initial adult depletion", "Baseline",
             initial = "adult_depletion"),
    scenario("single_male", "Single male release", "Releases",
             mc_release(20L, "adult_male", count = 250)),
    scenario("repeated_male", "Repeated male releases", "Releases", repeated),
    scenario("mated_female", "Mated female release", "Releases",
             mc_release(20L, "adult_female", count = 100, mate = "AA")),
    scenario("mixed_sex", "Mixed-sex release", "Releases",
             rbind(mc_release(20L, "adult_male", count = 100),
                   mc_release(20L, "adult_female", count = 100))),
    scenario("egg_release", "Egg release (substage 2)", "Releases",
             mc_release(20L, "egg", count = 2000, substage = 2L)),
    scenario("larval_release", "Larval release (substage 2)", "Releases",
             mc_release(20L, "larva", count = 1000, substage = 2L)),
    scenario("pupal_release", "Pupal release (substage 1)", "Releases",
             mc_release(20L, "pupa", count = 250, substage = 1L)),
    scenario("unmated_recovery", "Unmated females; males arrive day 20", "Releases",
             mc_release(20L, "adult_male", genotype = "AA", count = 200),
             initial = "unmated_only"),
    scenario("sterile_male", "Repeated sterile-male releases", "Genetics",
             transform(repeated, count = 500), cube = sterile),
    scenario("fitness_cost", "Release with genotype fitness costs", "Genetics",
             rbind(mc_release(20L, "adult_male", count = 100),
                   mc_release(20L, "adult_female", count = 100, mate = "AA")), cube = costly),
    scenario("homing_drive", "Illustrative homing-drive release", "Genetics",
             mc_release(20L, "adult_male", genotype = "Aa", count = 100), cube = homing),
    scenario("adult_mortality", "Adult mortality x3, days 30-59", "Interventions",
             intervention = "adult_mortality"),
    scenario("larvicide", "Larval mortality x3, days 30-59", "Interventions",
             intervention = "larvicide"),
    scenario("habitat", "Habitat K reduced 60%, days 30-59", "Interventions",
             intervention = "habitat"),
    scenario("fecundity", "Egg production reduced 60%, days 30-59", "Interventions",
             intervention = "fecundity"),
    scenario("combined", "Male releases + mortality + habitat", "Interventions", repeated,
             intervention = "combined"),
    scenario("seasonality", "Matched daily seasonal habitat", "Interventions",
             intervention = "seasonality"),
    scenario("control_beta_mismatch", "Control: native egg rate deliberately +5%", "Control",
             expected = "different", native_beta_multiplier = 1.05)
  )
}

mc_rates <- function(theta, K, day, intervention) {
  rates <- c(theta, list(K = K))
  active <- day >= 30L && day < 60L
  if (active && intervention %in% c("adult_mortality", "combined")) {
    rates$muM <- rates$muM * 3
    rates$muF <- rates$muF * 3
  }
  if (active && intervention == "larvicide") rates$muL <- rates$muL * 3
  if (active && intervention %in% c("habitat", "combined")) rates$K <- K * 0.4
  if (active && intervention == "fecundity") rates$beta <- rates$beta * 0.4
  if (intervention == "seasonality") {
    rates$K <- K * (1 + 0.6 * sin(2 * pi * day / 60))
  }
  rates
}

mc_mapping <- function(backend, places) {
  ix4 <- backend$index
  ix2 <- places$ix[[1L]]
  G <- backend$cube_info$G
  map <- rep(NA_integer_, length(places$u))
  for (g in seq_len(G)) {
    map[ix2$egg[, g]] <- ix4$egg_ix[, g, 1L]
    map[ix2$larvae[, g]] <- ix4$larv_ix[, g, 1L]
    map[ix2$pupae[, g]] <- ix4$pup_ix[, g, 1L]
    map[ix2$males[g]] <- ix4$male_ix[g, 1L]
    map[ix2$females_unmated[g]] <- ix4$unm_ix[g, 1L]
    for (gm in seq_len(G)) {
      map[ix2$females[g, gm, 1L]] <- ix4$fem_ix[g + (gm - 1L) * G, 1L]
    }
  }
  stopifnot(!anyNA(map), !anyDuplicated(map),
            identical(sort(map), seq_len(length(map))))
  map
}

mc_new_native_engine <- function(backend, rates) {
  ix <- backend$index
  cfg <- backend$cfg
  mc_native("ento_engine_create")(
    nNodes = 1L, nG = backend$cube_info$G,
    nE = cfg$nE, nL = cfg$nL, nP = cfg$nP,
    egg_ix = ix$egg_ix, larv_ix = ix$larv_ix, pup_ix = ix$pup_ix,
    male_ix = ix$male_ix, unm_ix = ix$unm_ix, fem_ix = ix$fem_ix,
    rE = rates$nE * rates$qE, rL = rates$nL * rates$qL,
    rP = rates$nP * rates$qP, muE = rates$muE, muL = rates$muL,
    muP = rates$muP, muM = rates$muM, muF = rates$muF,
    log_dd = TRUE, K = rates$K, gamma_dd = 0,
    beta = rates$beta, nu = rates$nu, omega = backend$omega,
    phi_cube = backend$phi, xiF = backend$xiF, xiM = backend$xiM,
    eta = backend$eta, B_mat = backend$B_mat,
    has_move = FALSE, move_probs_dense = matrix(0, 1L, 1L), move_rates = 0,
    muM_node_base = 1, muF_node_base = 1, tol = 1e-12,
    nState = ix$total_state_len
  )
}

mc_mg2_rhs <- function(places, transitions, S, cube, rates) {
  hazards <- MGDrivE2::spn_hazards(
    spn_P = places, spn_T = transitions, cube = cube, params = rates,
    type = "life", log_dd = TRUE, exact = FALSE, tol = 1e-12, verbose = FALSE
  )$hazards
  function(t, y, parms = NULL) {
    h <- vapply(hazards, function(hazard) hazard(t = t, M = y), numeric(1))
    list(as.numeric(S %*% h))
  }
}

mc_apply_mg2_release <- function(state, release, places, cube) {
  ix <- places$ix[[1L]]
  g <- match(release$genotype, cube$genotypesID)
  stage <- release$stage
  substage <- if (is.na(release$substage)) 1L else release$substage
  if (stage %in% c("egg", "larva", "pupa")) {
    block <- switch(stage, egg = ix$egg, larva = ix$larvae, pupa = ix$pupae)
    target <- block[substage, g]
  } else if (stage == "adult_male") {
    target <- ix$males[g]
  } else if (!is.na(release$mate_genotype)) {
    gm <- match(release$mate_genotype, cube$genotypesID)
    target <- ix$females[g, gm, 1L]
  } else {
    # Immediate mating at release, matching the native standalone API.
    weights <- state[ix$males] * cube$eta[g, ]
    if (sum(weights) > 0) {
      target <- ix$females[g, , 1L]
      state[target] <- state[target] + release$count * weights / sum(weights)
      return(state)
    }
    target <- ix$females_unmated[g]
  }
  state[target] <- state[target] + release$count
  state
}

mc_state_metadata <- function(places, cube) {
  ix <- places$ix[[1L]]
  out <- data.frame(compartment = places$u, stage = NA_character_,
                    genotype = NA_character_, stringsAsFactors = FALSE)
  for (g in seq_len(cube$genotypesN)) {
    for (stage in c("egg", "larva", "pupa", "adult_male", "adult_female")) {
      ids <- switch(stage, egg = ix$egg[, g], larva = ix$larvae[, g],
                    pupa = ix$pupae[, g], adult_male = ix$males[g],
                    adult_female = c(ix$females_unmated[g], ix$females[g, , 1L]))
      out$stage[ids] <- stage
      out$genotype[ids] <- cube$genotypesID[g]
    }
  }
  stopifnot(!anyNA(out))
  out
}

mc_run_scenario <- function(scenario, days = 120L, rtol = 1e-8, atol = 1e-10,
                            comparison_rtol = 1e-5, comparison_atol = 1e-4) {
  mc_require()
  theta <- mc_theta()
  cube <- scenario$cube
  p4 <- mc_parameters(theta, cube, rtol, atol)
  backend <- mc_native("mosquito_build_species_backend")(
    list(node1 = p4), 1L, "ode", days
  )
  places <- MGDrivE2::spn_P_lifecycle_node(theta, cube)
  transitions <- MGDrivE2::spn_T_lifecycle_node(places, theta, cube)
  S <- MGDrivE2::spn_S(places, transitions)
  eq2 <- MGDrivE2::equilibrium_lifeycle(theta, NF = 500, phi = 0.5,
                                      log_dd = TRUE, spn_P = places, cube = cube)
  map <- mc_mapping(backend, places)
  x4 <- backend$state
  x2 <- stats::setNames(x4[map], places$u)
  initial_difference <- max(abs(eq2$M0 - x2))
  K <- mc_native("mosquito_carrying_at")(backend$K_lookup, 0L)[[1L]]
  K_difference <- abs(eq2$params$K[[1L]] - K)
  stopifnot(initial_difference < 1e-7, K_difference < 1e-7)

  if (scenario$initial == "adult_depletion") {
    ids <- c(places$ix[[1L]]$males, places$ix[[1L]]$females,
             places$ix[[1L]]$females_unmated)
    x2[ids] <- x2[ids] * 0.3
  } else if (scenario$initial == "unmated_only") {
    x2[] <- 0
    x2[places$ix[[1L]]$females_unmated["AA"]] <- 200
  }
  x4[map] <- x2
  releases <- scenario$releases
  if (is.null(releases)) releases <- mc_release(integer(0), character(0))
  out2 <- out4 <- matrix(NA_real_, days + 1L, length(x2),
                         dimnames = list(NULL, places$u))
  out2[1L, ] <- x2
  out4[1L, ] <- x4[map]
  rhs_audit <- list()
  release_audit <- list()
  previous2 <- previous4 <- NULL
  engine <- solver <- rhs2 <- NULL
  runtime <- system.time({
    for (day in seq_len(days) - 1L) {
      rates2 <- mc_rates(theta, K, day, scenario$intervention)
      rates4 <- rates2
      rates4$beta <- rates4$beta * scenario$native_beta_multiplier
      if (!identical(rates2, previous2)) {
        rhs2 <- mc_mg2_rhs(places, transitions, S, cube, rates2)
        previous2 <- rates2
      }
      if (!identical(rates4, previous4)) {
        # Recreate the actual engine at parameter-change boundaries, carrying
        # forward its complete state. K is not recalibrated after interventions.
        engine <- mc_new_native_engine(backend, rates4)
        solver <- mc_native("create_ento_solver")(engine, x4, rtol, atol, 100000L)
        mc_native("solver_set_states")(solver, day, x4)
        previous4 <- rates4
      }
      # Compare both equations at the SAME state before every daily step.
      probe4 <- x4
      probe4[map] <- pmax(x2, 0)
      dx4 <- mc_native("ento_engine_dxdt")(engine, day, probe4)[map]
      dx2 <- rhs2(day, pmax(x2, 0))[[1L]]
      rhs_audit[[day + 1L]] <- data.frame(
        day = day, max_abs_difference = max(abs(dx4 - dx2)),
        max_scaled_difference = max(abs(dx4 - dx2) / pmax(1, abs(dx2)))
      )
      solved <- deSolve::ode(y = x2, times = c(day, day + 1), func = rhs2,
                            parms = NULL, method = "lsoda", rtol = rtol,
                            atol = atol, maxsteps = 100000L)
      x2 <- stats::setNames(as.numeric(solved[2L, -1L]), places$u)
      mc_native("solver_step")(solver)
      x4 <- as.numeric(mc_native("solver_get_states")(solver))
      due <- releases[releases$timestep == day + 1L, , drop = FALSE]
      if (nrow(due)) for (i in seq_len(nrow(due))) {
        before2 <- x2
        before4 <- x4
        x2 <- mc_apply_mg2_release(x2, due[i, ], places, cube)
        x4 <- mc_native("mosquito_apply_release")(x4, backend, due[i, ], 1L, FALSE)
        delta2 <- x2 - before2
        delta4 <- (x4 - before4)[map]
        release_audit[[length(release_audit) + 1L]] <- data.frame(
          day = day + 1L, stage = due$stage[i], genotype = due$genotype[i],
          requested = due$count[i], added_mg2 = sum(delta2), added_mg4 = sum(delta4),
          max_allocation_difference = max(abs(delta2 - delta4))
        )
        stopifnot(abs(sum(delta2) - due$count[i]) < 1e-7,
                  abs(sum(delta4) - due$count[i]) < 1e-7)
      }
      mc_native("solver_set_states")(solver, day + 1L, x4)
      out2[day + 2L, ] <- x2
      out4[day + 2L, ] <- x4[map]
    }
  })[["elapsed"]]
  stopifnot(all(is.finite(out2)), all(is.finite(out4)),
            min(out2) > -1e-6, min(out4) > -1e-6)
  difference <- abs(out4 - out2)
  allowance <- comparison_atol + comparison_rtol * pmax(abs(out2), abs(out4))
  passes <- all(difference <= allowance)
  metadata <- mc_state_metadata(places, cube)
  states <- data.frame(
    day = rep(0:days, times = ncol(out2)),
    compartment = rep(metadata$compartment, each = days + 1L),
    stage = rep(metadata$stage, each = days + 1L),
    genotype = rep(metadata$genotype, each = days + 1L),
    mg2 = as.vector(out2), mg4 = as.vector(out4),
    abs_difference = as.vector(difference),
    tolerance_ratio = as.vector(difference / allowance), stringsAsFactors = FALSE
  )
  totals <- stats::aggregate(cbind(mg2, mg4) ~ day + stage, states, sum)
  totals$abs_difference <- abs(totals$mg4 - totals$mg2)
  totals$relative_difference <- totals$abs_difference / pmax(1, abs(totals$mg2))
  genotypes <- stats::aggregate(cbind(mg2, mg4) ~ day + stage + genotype, states, sum)
  rhs_audit <- do.call(rbind, rhs_audit)
  summary <- data.frame(
    scenario = scenario$id, label = scenario$label, group = scenario$group,
    expected = scenario$expected, rtol = rtol, atol = atol, days = days,
    compartments = ncol(out2), K = K,
    initial_max_abs_difference = initial_difference, K_abs_difference = K_difference,
    rhs_max_abs_difference = max(rhs_audit$max_abs_difference),
    rhs_max_scaled_difference = max(rhs_audit$max_scaled_difference),
    max_compartment_abs_difference = max(difference),
    max_tolerance_ratio = max(difference / allowance),
    max_stage_abs_difference = max(totals$abs_difference),
    max_stage_relative_difference = max(totals$relative_difference),
    overlap_pass = passes,
    expectation_pass = if (scenario$expected == "overlap") passes else !passes,
    elapsed_seconds = runtime, stringsAsFactors = FALSE
  )
  list(summary = summary, states = states, totals = totals, genotypes = genotypes,
       rhs = rhs_audit, releases = if (length(release_audit)) do.call(rbind, release_audit),
       parameters = p4, theta = c(theta, list(K = K)), cube = cube,
       initial_state = stats::setNames(out4[1L, ], places$u))
}

mc_plot_results <- function(results, output) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 unavailable: CSV results saved without plots.", call. = FALSE)
    return(invisible(NULL))
  }
  totals <- do.call(rbind, lapply(results, function(x) {
    transform(x$totals, scenario = x$summary$scenario, label = x$summary$label,
              group = x$summary$group)
  }))
  long <- rbind(transform(totals, engine = "MGDrivE2", count = mg2),
                transform(totals, engine = "MGDrivE4 native", count = mg4))
  colors <- c(MGDrivE2 = "#16608d", "MGDrivE4 native" = "#e47b26")
  styles <- c(MGDrivE2 = "solid", "MGDrivE4 native" = "longdash")
  theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", size = 10),
                   plot.title = ggplot2::element_text(face = "bold"))
  plot_counts <- function(d, title, filename, ncol = 3L, height = 9) {
    p <- ggplot2::ggplot(d, ggplot2::aes(day, count, color = engine, linetype = engine)) +
      ggplot2::geom_line(linewidth = 0.85) +
      ggplot2::expand_limits(y = 0) +
      ggplot2::facet_wrap(~ label, scales = "free_y", ncol = ncol) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_linetype_manual(values = styles) +
      ggplot2::labs(title = title, subtitle = "Matched deterministic equations, full initial states and boundary events",
                    x = "Day", y = "Adult females", color = NULL, linetype = NULL) + theme
    ggplot2::ggsave(file.path(output, filename), p, width = 13, height = height, dpi = 170)
  }
  female <- long[long$stage == "adult_female", ]
  overview <- c("adult_depletion", "mated_female", "sterile_male", "homing_drive",
                "adult_mortality", "larvicide", "habitat", "fecundity", "combined")
  plot_counts(female[female$scenario %in% overview, ],
              "Mosquito population trajectories: MGDrivE2 versus native MGDrivE4",
              "trajectory_overview.png", height = 10)
  plot_counts(female[female$group %in% c("Releases", "Genetics"), ],
              "Release scenarios: adult female population", "release_trajectories.png", height = 12)
  plot_counts(female[female$group == "Interventions", ],
              "Mosquito-rate interventions: adult female population",
              "intervention_trajectories.png", height = 8)
  plot_counts(female[female$group %in% c("Baseline", "Control"), ],
              "Baseline checks and deliberate mismatch control", "baseline_and_control.png", height = 4.5)
  detail <- long[long$scenario %in% c("mated_female", "sterile_male", "larvicide", "combined"), ]
  p <- ggplot2::ggplot(detail, ggplot2::aes(day, count, color = engine, linetype = engine)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::facet_grid(stage ~ label, scales = "free_y") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_linetype_manual(values = styles) +
    ggplot2::labs(title = "Aquatic and adult trajectories", x = "Day", y = "Mosquito count",
                  color = NULL, linetype = NULL) + theme
  ggplot2::ggsave(file.path(output, "all_stage_trajectories.png"), p, width = 15, height = 12, dpi = 160)
  errors <- do.call(rbind, lapply(results, function(x) {
    d <- stats::aggregate(tolerance_ratio ~ day, x$states, max)
    transform(d, label = x$summary$label, group = x$summary$group)
  }))
  errors <- errors[errors$group != "Control", ]
  p <- ggplot2::ggplot(errors, ggplot2::aes(day, pmax(tolerance_ratio, 1e-10))) +
    ggplot2::geom_line(color = "#16608d", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 1, color = "#c23b38", linetype = "dashed") +
    ggplot2::scale_y_log10() + ggplot2::facet_wrap(~ label, ncol = 4) +
    ggplot2::labs(title = "Maximum compartment error relative to the acceptance limit",
                  subtitle = "Below the red line passes: |difference| <= 0.0001 + 0.00001 x population",
                  x = "Day", y = "Error / acceptance limit (log scale)") + theme
  ggplot2::ggsave(file.path(output, "compartment_errors.png"), p, width = 15, height = 13, dpi = 150)
  genetics <- do.call(rbind, lapply(results, function(x) {
    transform(x$genotypes, scenario = x$summary$scenario, label = x$summary$label)
  }))
  genetics <- genetics[genetics$scenario %in% c("repeated_male", "sterile_male", "fitness_cost", "homing_drive") &
                         genetics$stage == "adult_female", ]
  genetics <- rbind(transform(genetics, engine = "MGDrivE2", count = mg2),
                    transform(genetics, engine = "MGDrivE4 native", count = mg4))
  p <- ggplot2::ggplot(genetics, ggplot2::aes(day, count, color = genotype, linetype = engine)) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::facet_wrap(~ label, ncol = 2, scales = "free_y") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_linetype_manual(values = styles) +
    ggplot2::labs(title = "Adult females by genotype", x = "Day", y = "Adult females", linetype = NULL) + theme
  ggplot2::ggsave(file.path(output, "genotype_trajectories.png"), p, width = 13, height = 8, dpi = 170)
}
