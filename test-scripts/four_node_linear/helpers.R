four_node_movement <- function(cfg) {
  validate_four_node_config(cfg)
  distance <- abs(outer(cfg$x_km, cfg$x_km, "-"))
  dimnames(distance) <- list(cfg$sites, cfg$sites)
  weights <- exp(-distance / cfg$movement_decay_km)
  diag(weights) <- 0
  # Rows are origins, columns are destinations, conditional on a move.
  probs <- weights / rowSums(weights)
  rate <- cfg$p_move * 0.132 / (1 - cfg$p_move)
  stopifnot(all(abs(rowSums(probs) - 1) < 1e-12), all(diag(probs) == 0),
            all(probs[row(probs) != col(probs)] > 0))
  list(distance = distance, probs = probs, rates = rep(rate, 4L),
       mean_distance_km = mean(rowSums(probs * distance)))
}

four_node_cube <- function() {
  cube <- MGDrivE::cubeHomingDrive(
    cM = 1, cF = 1, chM = 0.95, chF = 0.95, crM = 0.05, crF = 0.05
  )
  gid <- as.character(cube$genotypesID)
  allele_s <- c(W = 1, H = 0.8, B = 0.5, R = 1)
  cube$s <- stats::setNames(vapply(strsplit(gid, ""), function(a) {
    if (any(!a %in% names(allele_s))) stop("Unexpected drive allele.", call. = FALSE)
    prod(allele_s[a])
  }, numeric(1)), gid)
  cube$b <- stats::setNames(ifelse(grepl("H", gid, fixed = TRUE), 0, 0.55), gid)
  cube$c <- stats::setNames(rep(1, length(gid)), gid)
  cube
}

four_node_parameters <- function(cfg, movement) {
  p <- malariasimulationGD::get_parameters(list(
    human_population = cfg$population,
    native_mosquito_backend = TRUE,
    individual_mosquitoes = TRUE, # native stochastic count-based tau-leap
    model_seasonality = FALSE,
    carrying_capacity = FALSE, # no time-dependent capacity rescaling
    progress_bar = FALSE,
    contact_multiplier = 1,
    human_slot_contact_multiplier = 1,
    human_mobility_enabled = FALSE,
    dem = 11,
    del = 6, dl = 21, dpl = 2,
    me = 0.05, ml = 0.15, mup = 0.05,
    b0 = 0.55,
    d1 = 0.160527,
    enable_heterogeneity = TRUE,
    sigma_squared = 1.67,
    n_heterogeneity_groups = 5L,
    incidence_rendering_min_ages = cfg$age_min_days,
    incidence_rendering_max_ages = cfg$age_max_days,
    clinical_incidence_rendering_min_ages = cfg$age_min_days,
    clinical_incidence_rendering_max_ages = cfg$age_max_days
  ))
  p$mum[] <- 0.132
  p$beta[] <- 16
  # Native stage controls are attached after get_parameters(), as in Godrick.
  p$native_mosquito_nEIP <- 6L
  p$native_mosquito_nE <- 2L
  p$native_mosquito_nL <- 3L
  p$native_mosquito_nP <- 2L
  p$native_mosquito_nu <- 6
  p$mosquito_tau_step <- cfg$mosquito_tau_step
  # Godrick routine treatment coverage is zero (the package default).
  stopifnot(sum(unlist(p$clinical_treatment_coverages)) == 0)
  p <- malariasimulationGD::set_bednets(
    p, timesteps = 1L, coverages = 0, retention = 1e12,
    dn0 = matrix(0.437, 1L, length(p$species)),
    rn = matrix(0.555, 1L, length(p$species)),
    rnm = matrix(0.24, 1L, length(p$species)), gamman = 1e12
  )
  p$initial_bednet_coverage <- 0.684
  p$initial_bednet_time <- 1L
  p$bednet_skip_timesteps <- 1L
  p$cube <- four_node_cube()
  p$move_probs <- movement$probs
  p$move_rates <- movement$rates
  p <- malariasimulationGD::set_equilibrium(p, init_EIR = cfg$init_eir,
                                           native_total_M = TRUE)
  stopifnot(isTRUE(p$native_total_M_equilibrium), !p$model_seasonality,
            !p$carrying_capacity)
  parameters <- lapply(cfg$sites, function(site) {
    q <- p
    if (site %in% cfg$release_sites) {
      q <- malariasimulationGD::set_releases(q, list(
        releasesStart = cfg$burnin_days + cfg$release_day, releasesNumber = 1L,
        releaseCount = cfg$release_count, releaseSex = "M",
        releaseGenotype = "HH", releasesInterval = 0L
      ))
    }
    q
  })
  names(parameters) <- cfg$sites
  four_node_assert_equivalent(parameters, cfg)
  parameters
}

four_node_assert_equivalent <- function(parameters, cfg) {
  if (!identical(names(parameters), cfg$sites)) {
    stop("Parameter lists must be ordered A, B, C, D.", call. = FALSE)
  }
  baseline <- lapply(parameters, function(p) {
    p$releases <- NULL
    p$releases_schedule <- NULL
    p
  })
  if (any(!vapply(baseline, identical, logical(1), baseline[[1L]]))) {
    stop("All site parameters must be identical apart from release schedules.", call. = FALSE)
  }
  for (p in baseline) {
    stopifnot(p$human_population == cfg$population,
              p$contact_multiplier == 1, p$human_slot_contact_multiplier == 1,
              !p$human_mobility_enabled, !p$model_seasonality, !p$carrying_capacity)
  }
  invisible(TRUE)
}

four_node_extract <- function(data, cfg) {
  if (!is.list(data) || length(data) != 4L) {
    stop("Expected four node output tables in A, B, C, D order.", call. = FALSE)
  }
  age_suffix <- paste(cfg$age_min_days, cfg$age_max_days, sep = "_")
  columns <- c(paste0("n_age_", age_suffix), paste0("n_inc_", age_suffix),
               paste0("n_inc_clinical_", age_suffix))
  total_days <- cfg$burnin_days + cfg$days
  absolute_release_day <- cfg$burnin_days + cfg$release_day
  daily <- carriers <- schedules <- vector("list", 4L)
  for (i in seq_len(4L)) {
    d <- data[[i]]
    if (!all(columns %in% names(d)) || nrow(d) != total_days ||
        !identical(as.integer(d$timestep), seq_len(total_days))) {
      stop("Missing age-specific incidence output or unexpected timesteps at ",
           cfg$sites[i], ".", call. = FALSE)
    }
    f <- attr(d, "mosquito_genotype_counts_female")
    m <- attr(d, "mosquito_genotype_counts_male")
    if (!is.matrix(f) || !is.matrix(m) || !identical(dim(f), dim(m)) ||
        nrow(f) != total_days || !identical(colnames(f), colnames(m)) ||
        any(!is.finite(c(f, m))) || any(c(f, m) < 0)) {
      stop("Missing or invalid adult genotype histories at ", cfg$sites[i], ".", call. = FALSE)
    }
    h <- grepl("H", colnames(f), fixed = TRUE)
    if (!any(h)) stop("No H-carrier genotype columns.", call. = FALSE)
    n_h <- rowSums(f[, h, drop = FALSE]) + rowSums(m[, h, drop = FALSE])
    n_adult <- rowSums(f) + rowSums(m)
    cf <- ifelse(n_adult > 0, n_h / n_adult, NA_real_)
    # No drive before the consumed release date, anywhere in the landscape.
    stopifnot(all(n_h[d$timestep < absolute_release_day] == 0),
              all(is.na(cf) | (cf >= 0 & cf <= 1)))
    carriers[[i]] <- data.frame(
      replicate = 1L, site = cfg$sites[i],
      timestep = d$timestep - cfg$burnin_days, simulation_day = d$timestep,
      gd_carriers = n_h, total_adults = n_adult, carrier_frequency = cf
    )
    daily[[i]] <- data.frame(
      replicate = 1L, site = cfg$sites[i],
      timestep = d$timestep - cfg$burnin_days, simulation_day = d$timestep,
      children = d[[columns[1L]]], infections = d[[columns[2L]]],
      clinical_cases = d[[columns[3L]]]
    )
    rs <- attr(d, "mosquito_release_schedule")
    if (!is.null(rs) && nrow(rs) > 0L) {
      rs$site <- cfg$sites[i]
      schedules[[i]] <- rs
    }
    if (cfg$sites[i] %in% cfg$release_sites) {
      # Histories are end-of-day counts: mortality/movement can reduce the
      # scheduled pulse before rendering. Its size is checked in the schedule.
      stopifnot("HH" %in% colnames(m), m[absolute_release_day, "HH"] > 0)
    }
  }
  schedule <- do.call(rbind, schedules)
  if (is.null(schedule) || nrow(schedule) != 2L ||
      !setequal(schedule$site, cfg$release_sites) ||
      any(schedule$timestep != absolute_release_day) ||
      any(schedule$count != cfg$release_count) ||
      any(schedule$sex != "M") || any(schedule$genotype != "HH")) {
    stop("Consumed release schedule does not match the two planned releases.", call. = FALSE)
  }
  schedule$simulation_day <- schedule$timestep
  schedule$timestep <- schedule$timestep - cfg$burnin_days
  carriers <- do.call(rbind, carriers)
  daily <- do.call(rbind, daily)
  list(carriers = carriers[carriers$timestep > 0L, ],
       daily = daily[daily$timestep > 0L, ], schedule = schedule,
       burnin_carriers = carriers[carriers$timestep <= 0L, ],
       burnin_daily = daily[daily$timestep <= 0L, ])
}

four_node_baseline_diagnostics <- function(outputs, cfg) {
  daily <- rbind(outputs$burnin_daily, outputs$daily)
  carriers <- rbind(outputs$burnin_carriers, outputs$carriers)
  daily <- daily[daily$simulation_day < cfg$burnin_days + cfg$release_day, ]
  daily$period <- ifelse(daily$timestep <= 0L,
    paste0("burnin_year_", ceiling(daily$simulation_day / 365)), "observed_before_release")
  do.call(rbind, lapply(split(daily, interaction(daily$site, daily$period, drop = TRUE)), function(d) {
    c <- carriers[carriers$site == d$site[1L] & carriers$simulation_day %in% d$simulation_day, ]
    stopifnot(all(c$gd_carriers == 0))
    data.frame(site = d$site[1L], period = d$period[1L],
      simulation_start_day = min(d$simulation_day), simulation_end_day = max(d$simulation_day),
      days = nrow(d), mean_children = mean(d$children), mean_adults = mean(c$total_adults),
      infections = sum(d$infections), person_years = sum(d$children) / 365,
      infection_incidence_ppy = if (sum(d$children) > 0) 365 * sum(d$infections) / sum(d$children) else NA_real_)
  }))
}

four_node_incidence <- function(daily, cfg) {
  daily$window <- (daily$timestep - 1L) %/% cfg$incidence_window_days
  groups <- split(daily, interaction(daily$site, daily$window, drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(x) {
    person_days <- sum(x$children)
    data.frame(
      replicate = 1L, site = x$site[1L],
      start_day = min(x$timestep), end_day = max(x$timestep),
      midpoint_day = mean(range(x$timestep)), person_days = person_days,
      infections = sum(x$infections), clinical_cases = sum(x$clinical_cases),
      infection_incidence_ppy = if (person_days > 0) 365 * sum(x$infections) / person_days else NA_real_,
      clinical_incidence_ppy = if (person_days > 0) 365 * sum(x$clinical_cases) / person_days else NA_real_
    )
  }))
  rownames(out) <- NULL
  out[order(match(out$site, cfg$sites), out$start_day), ]
}

four_node_plots <- function(out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install suggested package ggplot2 to create the plots.", call. = FALSE)
  }
  context <- readRDS(file.path(out_dir, "context.rds"))
  cfg <- context$config
  carriers <- utils::read.csv(file.path(out_dir, "carrier_frequency.csv"))
  incidence <- utils::read.csv(file.path(out_dir, "incidence_windows.csv"))
  colours <- c(A = "#0072B2", B = "#009E73", C = "#D55E00", D = "#CC79A7")
  make_plot <- function(sites, kind) {
    if (kind == "carrier") {
      d <- carriers[carriers$site %in% sites, ]
      d$year <- d$timestep / 365
      d$value <- d$carrier_frequency
      ylab <- "GD carrier frequency (all adult mosquitoes)"
      title <- "Adult mosquito GD carriers"
    } else {
      d <- incidence[incidence$site %in% sites, ]
      d$year <- d$midpoint_day / 365
      d$value <- d$infection_incidence_ppy
      ylab <- "New infections per child-year"
      title <- "Infection incidence: 6 months to <10 years"
    }
    p <- ggplot2::ggplot(d, ggplot2::aes(x = year, y = value, colour = site)) +
      ggplot2::geom_vline(xintercept = cfg$release_day / 365, linetype = 2,
                          colour = "grey45", linewidth = 0.5) +
      ggplot2::scale_colour_manual(values = colours, breaks = sites) +
      ggplot2::scale_x_continuous(limits = c(0, cfg$days / 365),
                                breaks = seq(0, cfg$days / 365, by = 1)) +
      ggplot2::labs(
        title = title,
        subtitle = paste0("Sites ", paste(sites, collapse = " & "),
                          if (identical(sites, c("A", "C"))) " (release sites)" else " (no direct releases)"),
        x = "Years from observation start (after burn-in)", y = ylab, colour = "Site",
        caption = if (kind == "carrier") {
          "H-carrier genotypes: WH, HH, HR, HB; females and males combined."
        } else {
          sprintf("%d-day windows; cases / observed child person-time; final window may be shorter.",
                  cfg$incidence_window_days)
        }
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"),
                     plot.caption = ggplot2::element_text(size = 8))
    if (all(table(d$site) > 1L)) p <- p + ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE)
    else p <- p + ggplot2::geom_point(size = 2, na.rm = TRUE)
    if (kind == "carrier") p <- p + ggplot2::scale_y_continuous(limits = c(0, 1))
    else {
      ymax <- max(incidence$infection_incidence_ppy, na.rm = TRUE)
      if (!is.finite(ymax) || ymax <= 0) ymax <- 1
      p <- p + ggplot2::scale_y_continuous(limits = c(0, ymax * 1.05))
    }
    p
  }
  plots <- list(
    carrier_B_D = make_plot(c("B", "D"), "carrier"),
    carrier_A_C = make_plot(c("A", "C"), "carrier"),
    incidence_B_D = make_plot(c("B", "D"), "incidence"),
    incidence_A_C = make_plot(c("A", "C"), "incidence")
  )
  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(plots)) {
    for (ext in c("png", "pdf")) {
      ggplot2::ggsave(file.path(plot_dir, paste0(name, ".", ext)), plots[[name]],
                      width = 8, height = 4.8, units = "in", dpi = 160)
    }
  }
  grDevices::png(file.path(plot_dir, "overview.png"), width = 2000, height = 1300, res = 140)
  tryCatch({
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      3L, 2L, heights = grid::unit(c(0.07, 0.465, 0.465), "npc"))))
    grid::grid.text(
      sprintf("A - B - C - D | %g-year burn-in | %d humans/site | release at A and C: observed day %d (dashed)",
              cfg$burnin_days / 365, cfg$population, cfg$release_day),
      vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1:2),
      gp = grid::gpar(fontsize = 14, fontface = "bold"))
    for (i in seq_along(plots)) print(plots[[i]], newpage = FALSE,
      vp = grid::viewport(layout.pos.row = 2L + (i - 1L) %/% 2L,
                          layout.pos.col = 1L + (i - 1L) %% 2L))
  }, finally = grDevices::dev.off())
  invisible(plots)
}
