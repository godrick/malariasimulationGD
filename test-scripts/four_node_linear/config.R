# Godrick/P20 settings, copied explicitly so this example does not depend on
# the untracked Godrick_calibration bundle or its machine-specific paths.
four_node_config <- function() {
  list(
    sites = c("A", "B", "C", "D"),
    x_km = c(0, 2, 4, 6),
    population = 10000L,
    burnin_days = 5L * 365L,
    days = 5L * 365L, # observed days, following the burn-in
    release_day = 183L, # observed day; absolute release day is 1825 + 183
    release_sites = c("A", "C"),
    release_count = 10000L, # one male HH release per release site
    seed = 20261029L,
    init_eir = 18.4730728587,
    movement_decay_km = 2,
    p_move = 0.02, # Godrick competing movement-versus-mortality probability
    incidence_window_days = 30L,
    age_min_days = 182L,
    age_max_days = 3649L, # six months to strictly under ten years
    mosquito_tau_step = 0.1
  )
}

validate_four_node_config <- function(cfg) {
  for (key in c("population", "days", "release_day", "release_count", "seed",
                "incidence_window_days", "age_min_days", "age_max_days")) {
    x <- cfg[[key]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || x > .Machine$integer.max || x != as.integer(x)) {
      stop(key, " must be a positive integer.", call. = FALSE)
    }
  }
  if (!is.numeric(cfg$burnin_days) || length(cfg$burnin_days) != 1L ||
      !is.finite(cfg$burnin_days) || cfg$burnin_days < 0 ||
      cfg$burnin_days > .Machine$integer.max - cfg$days ||
      cfg$burnin_days != floor(cfg$burnin_days)) {
    stop("burnin_days must be a nonnegative integer with burnin_days + days within integer range.",
         call. = FALSE)
  }
  if (cfg$release_day > cfg$days) {
    stop("release_day must not exceed days.", call. = FALSE)
  }
  if (!identical(cfg$sites, c("A", "B", "C", "D")) ||
      !identical(as.numeric(cfg$x_km), c(0, 2, 4, 6))) {
    stop("Site order and positions must be A=0, B=2, C=4, D=6 km.", call. = FALSE)
  }
  if (!identical(cfg$release_sites, c("A", "C"))) {
    stop("Release sites must be A and C.", call. = FALSE)
  }
  if (cfg$age_min_days >= cfg$age_max_days) {
    stop("Age bounds must be increasing.", call. = FALSE)
  }
  for (key in c("init_eir", "movement_decay_km", "mosquito_tau_step", "p_move")) {
    x <- cfg[[key]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop(key, " must be a positive finite number.", call. = FALSE)
    }
  }
  if (cfg$p_move >= 1 || cfg$mosquito_tau_step > 1) {
    stop("p_move must be below 1 and mosquito_tau_step must be at most 1.", call. = FALSE)
  }
  invisible(cfg)
}
