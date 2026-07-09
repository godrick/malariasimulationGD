# ------------------------------------------------------------------------------
# busia_landscape.R
# ------------------------------------------------------------------------------
# Busia 58-cluster landscape and demographic helpers.
#
# The coordinate file contains cluster ids and latitude/longitude only. Because
# no cluster population column is available, the default workflow uses a
# deterministic synthetic village-population vector with the requested Busia
# regional total and cross-village dispersion. Set BUSIA_TOTAL_POPULATION or
# BUSIA_POPULATION_SD to override those simulation inputs.
# ------------------------------------------------------------------------------

.busia_config_file <- local({
  frames <- sys.frames()
  ofiles <- vapply(frames, function(env) {
    val <- env$ofile
    if (is.null(val)) NA_character_ else as.character(val)
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) normalizePath(tail(ofiles, 1L), mustWork = FALSE) else NA_character_
})

.busia_example_dir <- function() {
  if (!is.na(.busia_config_file)) {
    return(normalizePath(file.path(dirname(.busia_config_file), ".."), mustWork = TRUE))
  }
  normalizePath(file.path("test-scripts", "Busia_landscape"), mustWork = TRUE)
}

.busia_data_path <- function(filename) {
  file.path(.busia_example_dir(), "Busia_data", filename)
}

busia_default_total_population <- function() {
  val <- Sys.getenv("BUSIA_TOTAL_POPULATION", unset = NA_character_)
  if (!is.na(val) && nzchar(val)) {
    out <- suppressWarnings(as.integer(val))
    if (is.na(out) || out < 58L) {
      stop("BUSIA_TOTAL_POPULATION must be an integer >= 58.", call. = FALSE)
    }
    return(out)
  }
  44875L
}

busia_default_population_sd <- function() {
  val <- Sys.getenv("BUSIA_POPULATION_SD", unset = NA_character_)
  if (!is.na(val) && nzchar(val)) {
    out <- suppressWarnings(as.numeric(val))
    if (!is.finite(out) || out < 0) {
      stop("BUSIA_POPULATION_SD must be a finite non-negative number.",
           call. = FALSE)
    }
    return(out)
  }
  60.3
}

busia_default_target_pfpr <- function() {
  val <- Sys.getenv("BUSIA_TARGET_PFPR", unset = NA_character_)
  if (!is.na(val) && nzchar(val)) {
    out <- suppressWarnings(as.numeric(val))
    if (!is.finite(out) || out <= 0 || out >= 1) {
      stop("BUSIA_TARGET_PFPR must be a finite number in (0, 1).",
           call. = FALSE)
    }
    return(out)
  }
  0.39
}

busia_largest_remainder_counts <- function(counts, sample_size) {
  counts <- as.numeric(counts)
  sample_size <- as.integer(sample_size)
  if (length(counts) < 1L || any(!is.finite(counts)) || any(counts < 0)) {
    stop("counts must be a non-negative numeric vector.", call. = FALSE)
  }
  if (sample_size < 1L || sum(counts) <= 0) {
    stop("sample_size must be positive and counts must have positive total.",
         call. = FALSE)
  }
  raw <- sample_size * counts / sum(counts)
  out <- floor(raw)
  remainder <- sample_size - sum(out)
  if (remainder > 0L) {
    add <- order(raw - out, decreasing = TRUE)[seq_len(remainder)]
    out[add] <- out[add] + 1L
  }
  as.integer(out)
}

busia_deathrates_for_target_age_counts <- function(target_counts, age_width_days) {
  n_age <- length(target_counts)
  if (length(age_width_days) != n_age) {
    stop("age_width_days must have the same length as target_counts.",
         call. = FALSE)
  }
  if (any(target_counts <= 0L)) {
    stop("All target age groups must have positive counts.", call. = FALSE)
  }
  aging_rate <- 1 / age_width_days
  aging_rate[[n_age]] <- 0

  deathrates <- rep(1e-5, n_age)
  for (i in 2:n_age) {
    ratio <- target_counts[[i]] / target_counts[[i - 1L]]
    deathrates[[i]] <- aging_rate[[i - 1L]] / ratio - aging_rate[[i]]
  }
  if (any(!is.finite(deathrates)) || any(deathrates <= 0) || any(deathrates >= 1)) {
    stop("Could not derive valid daily death rates for the requested Busia age structure.",
         call. = FALSE)
  }
  deathrates
}

read_busia_age_structure <- function(
    path = .busia_data_path("busia_age_structure.csv")) {
  if (!file.exists(path)) {
    stop(sprintf("Busia age-structure file not found: %s", path), call. = FALSE)
  }
  age <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("age_group", "age_min_years", "age_max_years", "population")
  missing <- setdiff(required, names(age))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Busia age-structure file is missing required columns: %s",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  age$age_min_years <- as.numeric(age$age_min_years)
  age$age_max_years <- as.numeric(age$age_max_years)
  age$population <- as.integer(age$population)
  age <- age[order(age$age_min_years), , drop = FALSE]
  if (nrow(age) != 9L ||
      any(!is.finite(age$age_min_years)) ||
      any(!is.finite(age$age_max_years)) ||
      any(age$age_max_years <= age$age_min_years) ||
      any(is.na(age$population)) ||
      any(age$population <= 0L)) {
    stop("Busia age-structure data must contain 9 positive, increasing age groups.",
         call. = FALSE)
  }
  age
}

busia_demography_spec <- function(total_population = busia_default_total_population()) {
  age <- read_busia_age_structure()
  total_population <- as.integer(total_population)
  if (is.na(total_population) || total_population < nrow(age)) {
    stop("total_population must be an integer at least as large as the number of age groups.",
         call. = FALSE)
  }
  year <- 365L
  age_min_days <- as.integer(age$age_min_years * year)
  age_max_days <- as.integer(age$age_max_years * year)
  age_width_days <- age_max_days - age_min_days
  scaled_counts <- busia_largest_remainder_counts(age$population, total_population)
  deathrates <- busia_deathrates_for_target_age_counts(
    target_counts = scaled_counts,
    age_width_days = age_width_days
  )
  list(
    age_structure = age,
    total_population = total_population,
    age_min_days = age_min_days,
    age_max_days = age_max_days,
    age_width_days = age_width_days,
    scaled_counts = scaled_counts,
    deathrates = deathrates
  )
}

apply_busia_demography <- function(parameters,
                                   total_population = parameters$human_population,
                                   render_age_groups = TRUE) {
  spec <- busia_demography_spec(total_population = total_population)
  parameters <- malariasimulationGD::set_demography(
    parameters,
    agegroups = spec$age_max_days,
    timesteps = 0L,
    deathrates = matrix(spec$deathrates, nrow = 1L)
  )
  if (isTRUE(render_age_groups)) {
    parameters$age_group_rendering_min_ages <- spec$age_min_days
    parameters$age_group_rendering_max_ages <- spec$age_max_days
  }
  parameters
}

read_busia_coordinates <- function(
    path = .busia_data_path("busia_58_cluster_coordinates.csv")) {
  if (!file.exists(path)) {
    stop(sprintf("Busia coordinate file not found: %s", path), call. = FALSE)
  }
  coords <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("cluster_id", "latitude", "longitude")
  missing <- setdiff(required, names(coords))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Busia coordinate file is missing required columns: %s",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  coords$cluster_id <- as.integer(coords$cluster_id)
  coords$latitude <- as.numeric(coords$latitude)
  coords$longitude <- as.numeric(coords$longitude)
  if (nrow(coords) != 58L) {
    stop(sprintf("Expected 58 Busia coordinate rows, found %d.", nrow(coords)),
         call. = FALSE)
  }
  if (anyNA(coords$cluster_id) || anyDuplicated(coords$cluster_id)) {
    stop("Busia coordinate cluster_id values must be unique non-missing integers.",
         call. = FALSE)
  }
  if (any(!is.finite(coords$latitude)) || any(!is.finite(coords$longitude))) {
    stop("Busia latitude/longitude values must be finite.", call. = FALSE)
  }
  coords
}

busia_adjust_integer_sum <- function(values, target_sum) {
  values <- as.integer(values)
  target_sum <- as.integer(target_sum)
  delta <- target_sum - sum(values)
  if (delta == 0L) {
    return(values)
  }
  center <- target_sum / length(values)
  ord <- order(abs(values - center))
  step <- if (delta > 0L) 1L else -1L
  for (i in seq_len(abs(delta))) {
    values[ord[((i - 1L) %% length(ord)) + 1L]] <-
      values[ord[((i - 1L) %% length(ord)) + 1L]] + step
  }
  values
}

busia_allocate_population <- function(n_nodes,
                                      total_population = busia_default_total_population(),
                                      population_sd = busia_default_population_sd()) {
  if (n_nodes < 1L) {
    stop("n_nodes must be positive.", call. = FALSE)
  }
  total_population <- as.integer(total_population)
  population_sd <- as.numeric(population_sd)
  if (is.na(total_population) || total_population < n_nodes) {
    stop("total_population must be an integer at least as large as n_nodes.",
         call. = FALSE)
  }
  if (!is.finite(population_sd) || population_sd < 0) {
    stop("population_sd must be a finite non-negative number.", call. = FALSE)
  }
  if (population_sd == 0 || n_nodes == 1L) {
    return(busia_largest_remainder_counts(rep(1, n_nodes), total_population))
  }

  # Deterministic quantiles give a stable heterogeneous population input while
  # preserving the requested total after integer rounding.
  z <- stats::qnorm((seq_len(n_nodes) - 0.5) / n_nodes)
  z <- as.numeric(scale(z))
  out <- round(total_population / n_nodes + population_sd * z)
  out <- busia_adjust_integer_sum(out, total_population)
  if (any(out <= 0L)) {
    stop(
      "Requested Busia population SD is too large for the total population; ",
      "one or more node populations would be non-positive.",
      call. = FALSE
    )
  }
  perm <- order((seq_len(n_nodes) * 0.6180339887498949) %% 1)
  as.integer(out[perm])
}

build_busia_landscape <- function(total_population = busia_default_total_population()) {
  coords <- read_busia_coordinates()
  lat0 <- mean(coords$latitude)
  lon0 <- mean(coords$longitude)
  x <- (coords$longitude - lon0) * 111.320 * cos(lat0 * pi / 180)
  y <- (coords$latitude - lat0) * 110.574

  nodes <- data.frame(
    node = seq_len(nrow(coords)),
    cluster_id = coords$cluster_id,
    latitude = coords$latitude,
    longitude = coords$longitude,
    x = x,
    y = y,
    NH_per_node = busia_allocate_population(nrow(coords), total_population),
    stringsAsFactors = FALSE
  )
  D <- as.matrix(stats::dist(nodes[, c("x", "y")]))
  validate_busia_landscape(list(nodes = nodes, D = D, n_nodes = nrow(nodes)))
  list(nodes = nodes, D = D, n_nodes = nrow(nodes))
}

validate_busia_landscape <- function(land) {
  if (!is.list(land) || is.null(land$nodes) || is.null(land$D)) {
    stop("Busia landscape must be a list with `nodes` and `D`.", call. = FALSE)
  }
  nodes <- land$nodes
  D <- as.matrix(land$D)
  required <- c("node", "cluster_id", "latitude", "longitude", "x", "y", "NH_per_node")
  missing <- setdiff(required, names(nodes))
  if (length(missing) > 0L) {
    stop(sprintf("Busia nodes are missing columns: %s", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  if (nrow(nodes) != 58L || !identical(as.integer(land$n_nodes), 58L)) {
    stop("Busia landscape must contain exactly 58 nodes.", call. = FALSE)
  }
  if (!all(dim(D) == c(58L, 58L))) {
    stop("Busia distance matrix must be 58 x 58.", call. = FALSE)
  }
  if (any(!is.finite(D)) || any(D < 0) ||
      max(abs(D - t(D))) > 1e-8 ||
      any(abs(diag(D)) > 1e-8)) {
    stop("Busia distance matrix must be finite, symmetric, non-negative, with zero diagonal.",
         call. = FALSE)
  }
  if (any(is.na(nodes$NH_per_node)) || any(nodes$NH_per_node <= 0L)) {
    stop("Busia NH_per_node must be positive for every node.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_busia_demography_for_NH <- function(NH) {
  NH <- as.integer(NH)
  if (length(NH) != 58L || anyNA(NH) || any(NH <= 0L)) {
    stop("Busia NH must be a positive integer vector of length 58.", call. = FALSE)
  }
  spec <- busia_demography_spec(total_population = sum(NH))
  if (sum(spec$scaled_counts) != sum(NH)) {
    stop("Internal Busia demography scaling failed to preserve total population.",
         call. = FALSE)
  }
  invisible(TRUE)
}

busia_allowed_matrix <- function(n_nodes = 58L) {
  allowed <- matrix(TRUE, n_nodes, n_nodes)
  diag(allowed) <- FALSE
  allowed
}

validate_busia_simulation_output <- function(res, n_nodes = 58L,
                                             required_columns = character(),
                                             require_genotypes = FALSE,
                                             label = "simulation output") {
  data_list <- if (is.list(res) && all(vapply(res, is.data.frame, logical(1)))) {
    res
  } else if (!is.null(res$data)) {
    res$data
  } else if (!is.null(res$sim$data)) {
    res$sim$data
  } else {
    stop(sprintf("Could not locate per-node data frames in %s.", label), call. = FALSE)
  }
  if (length(data_list) != n_nodes) {
    stop(sprintf("%s has %d node outputs; expected %d.",
                 label, length(data_list), n_nodes), call. = FALSE)
  }
  row_counts <- vapply(data_list, nrow, integer(1))
  if (any(row_counts < 1L) || length(unique(row_counts)) != 1L) {
    stop(sprintf("%s node outputs must all have the same positive number of rows.", label),
         call. = FALSE)
  }
  for (i in seq_along(data_list)) {
    missing <- setdiff(required_columns, names(data_list[[i]]))
    if (length(missing) > 0L) {
      stop(sprintf("%s node %d is missing columns: %s",
                   label, i, paste(missing, collapse = ", ")), call. = FALSE)
    }
    if (isTRUE(require_genotypes)) {
      female <- attr(data_list[[i]], "mosquito_genotype_counts_female")
      male <- attr(data_list[[i]], "mosquito_genotype_counts_male")
      if (is.null(female) || is.null(male) ||
          nrow(female) != row_counts[[i]] ||
          nrow(male) != row_counts[[i]]) {
        stop(sprintf(
          "%s node %d is missing genotype count attributes with matching rows.",
          label, i
        ), call. = FALSE)
      }
    }
  }
  invisible(data_list)
}

# Backward-compatible aliases for any Busia-local code that still uses the
# seven-node example function names.
build_seven_node_landscape <- build_busia_landscape
