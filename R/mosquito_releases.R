#' @title Configure scheduled mosquito genotype releases (hybrid mode)
#'
#' @description
#' Expand and validate a simple MGDrivE-like release specification for adult
#' mosquito genotype releases in hybrid mosquito mode
#' (`individual_mosquitoes = TRUE`). The expanded schedule is stored in
#' `parameters$releases_schedule` and the original configuration is preserved in
#' `parameters$releases`.
#'
#' The release schedule is a data.frame with columns:
#' `timestep`, `species`, `sex`, `genotype`, `count`.
#'
#' `releaseSex = "M"` updates the internal male genotype pool (males are not
#' explicit individual agents in the current hybrid implementation).
#' `releaseSex = "F"` releases susceptible adult females (`Sm`) as individual
#' mosquitoes.
#'
#' @param parameters model parameters
#' @param releases a release configuration list, or `NULL` to clear releases
#' @return modified parameter list
#' @export
set_releases <- function(parameters, releases) {
  if (!is.list(parameters)) {
    stop("parameters must be a list")
  }

  if (is.null(releases)) {
    parameters$releases <- NULL
    parameters$releases_schedule <- NULL
    return(parameters)
  }
  if (!is.list(releases)) {
    stop("releases must be a list or NULL")
  }

  parameters$releases <- releases
  parameters$releases_schedule <- expand_mosquito_release_schedule(parameters, releases)
  parameters
}

genotype_debug_focus_timesteps <- function(parameters) {
  tt <- parameters$debug_genotype_timesteps
  if (!is.null(tt)) {
    return(as.integer(tt))
  }
  sch <- parameters$releases_schedule
  if (!is.null(sch) && nrow(sch) > 0L) {
    return(sort(unique(c(sch$timestep - 1L, sch$timestep, sch$timestep + 1L))))
  }
  NULL
}

genotype_debug_enabled <- function(parameters, timestep = NULL) {
  if (!isTRUE(parameters$debug_genotypes)) {
    return(FALSE)
  }
  if (is.null(timestep)) {
    return(TRUE)
  }
  focus <- genotype_debug_focus_timesteps(parameters)
  if (is.null(focus)) {
    return(TRUE)
  }
  timestep %in% focus
}

genotype_debug_fmt_counts <- function(x, names_hint = NULL) {
  if (length(x) == 0L) {
    return("[]")
  }
  if (is.null(names(x)) && !is.null(names_hint) && length(names_hint) == length(x)) {
    names(x) <- names_hint
  }
  if (is.null(names(x))) {
    return(paste0("[", paste(x, collapse = ","), "]"))
  }
  paste(paste0(names(x), "=", x), collapse = ",")
}

genotype_debug_species_counts <- function(variables, models, parameters, species_i) {
  species_name <- parameters$species[[species_i]]
  male_counts <- NULL
  if (!is.null(models[[species_i]]$genotype_state)) {
    male_counts <- models[[species_i]]$genotype_state$male_counts
  }
  if (is.null(models[[species_i]]$cube)) {
    female_counts <- c(WT = variables$species$get_index_of(species_name)$and(
      variables$mosquito_state$get_index_of("NonExistent")$not(TRUE)
    )$size())
    if (is.null(male_counts)) {
      male_counts <- c(WT = 0L)
    }
    return(list(female = female_counts, male = male_counts))
  }
  cube_info <- cube_genotype_info(models[[species_i]]$cube)
  female_counts <- adult_female_genotype_counts_by_species(variables, species_name, cube_info$G)
  names(female_counts) <- cube_info$genotypesID
  if (is.null(male_counts)) {
    male_counts <- rep.int(0L, cube_info$G)
    names(male_counts) <- cube_info$genotypesID
  } else if (is.null(names(male_counts))) {
    names(male_counts) <- cube_info$genotypesID
  }
  list(female = female_counts, male = male_counts)
}

genotype_debug_log_counts <- function(parameters, timestep, tag, species_name, counts, extra = NULL) {
  if (!genotype_debug_enabled(parameters, timestep)) {
    return(invisible(NULL))
  }
  line <- sprintf(
    "[geno-debug t=%d %s %s] F{%s} M{%s}",
    timestep,
    tag,
    species_name,
    genotype_debug_fmt_counts(counts$female),
    genotype_debug_fmt_counts(counts$male)
  )
  if (!is.null(extra) && nzchar(extra)) {
    line <- paste0(line, " ", extra)
  }
  cat(line, "\n")
}

genotype_debug_log <- function(parameters, timestep, tag, species_name = NULL, msg) {
  if (!genotype_debug_enabled(parameters, timestep)) {
    return(invisible(NULL))
  }
  prefix <- if (is.null(species_name)) {
    sprintf("[geno-debug t=%d %s]", timestep, tag)
  } else {
    sprintf("[geno-debug t=%d %s %s]", timestep, tag, species_name)
  }
  cat(prefix, msg, "\n")
}

empty_mosquito_release_schedule <- function() {
  data.frame(
    timestep = integer(0),
    species = character(0),
    sex = character(0),
    genotype = character(0),
    count = integer(0),
    stringsAsFactors = FALSE
  )
}

coerce_scalar_integer <- function(x, name) {
  if (length(x) != 1 || is.na(x) || !is.numeric(x)) {
    stop(sprintf("%s must be a single numeric value", name))
  }
  xi <- as.integer(x)
  if (!isTRUE(all.equal(as.numeric(xi), as.numeric(x)))) {
    stop(sprintf("%s must be an integer value", name))
  }
  xi
}

resolve_release_species <- function(parameters, releases) {
  release_species <- releases$releaseSpecies
  if (is.null(release_species)) {
    if (!is.null(parameters$species)) {
      if (length(parameters$species) != 1L) {
        stop("releaseSpecies must be provided when multiple mosquito species are configured")
      }
      release_species <- as.character(parameters$species[[1]])
    } else {
      release_species <- "gamb"
    }
  }
  if (length(release_species) != 1L) {
    stop("releaseSpecies must be a single character value")
  }
  release_species <- as.character(release_species[[1]])
  if (!is.null(parameters$species) && !(release_species %in% parameters$species)) {
    stop(sprintf(
      "releaseSpecies '%s' is not in parameters$species (%s)",
      release_species,
      paste(parameters$species, collapse = ", ")
    ))
  }
  release_species
}

resolve_release_genotype <- function(parameters, releases, cube_info) {
  release_genotype <- releases$releaseGenotype
  if (is.null(release_genotype)) {
    if (is.null(parameters$cube) || is.null(parameters$cube$releaseType)) {
      stop(sprintf(
        "releaseGenotype is required when cube$releaseType is not set. Allowed genotypes: %s",
        paste(cube_info$genotypesID, collapse = ", ")
      ))
    }
    release_genotype <- parameters$cube$releaseType
  }

  if (length(release_genotype) != 1L) {
    stop("releaseGenotype must be a single character value")
  }
  release_genotype <- as.character(release_genotype[[1]])
  if (!(release_genotype %in% cube_info$genotypesID)) {
    stop(sprintf(
      "releaseGenotype '%s' is not in cube$genotypesID. Allowed genotypes: %s",
      release_genotype,
      paste(cube_info$genotypesID, collapse = ", ")
    ))
  }
  release_genotype
}

expand_mosquito_release_schedule <- function(parameters, releases) {
  releases_start <- releases$releasesStart
  release_count <- releases$releaseCount
  if (is.null(releases_start)) {
    stop("releases$releasesStart is required")
  }
  if (is.null(release_count)) {
    stop("releases$releaseCount is required")
  }

  releases_start <- coerce_scalar_integer(releases_start, "releases$releasesStart")
  release_count <- coerce_scalar_integer(release_count, "releases$releaseCount")
  releases_number <- releases$releasesNumber
  if (is.null(releases_number)) {
    releases_number <- 1L
  } else {
    releases_number <- coerce_scalar_integer(releases_number, "releases$releasesNumber")
  }
  releases_interval <- releases$releasesInterval
  if (is.null(releases_interval)) {
    releases_interval <- 0L
  } else {
    releases_interval <- coerce_scalar_integer(releases_interval, "releases$releasesInterval")
  }

  if (releases_start < 1L) {
    stop("releases$releasesStart must be >= 1")
  }
  if (releases_number < 0L) {
    stop("releases$releasesNumber must be >= 0")
  }
  if (releases_interval < 0L) {
    stop("releases$releasesInterval must be >= 0")
  }
  if (release_count < 0L) {
    stop("releases$releaseCount must be >= 0")
  }

  if (releases_number == 0L || release_count == 0L) {
    return(empty_mosquito_release_schedule())
  }

  if (!isTRUE(parameters$individual_mosquitoes) && !native_mosquito_backend_enabled(parameters)) {
    stop("Mosquito genotype releases require stochastic mosquitoes or the native mosquito backend.")
  }
  if (is.null(parameters$cube)) {
    stop("Mosquito genotype releases require parameters$cube (cannot release non-WT genotypes without a cube)")
  }

  cube_info <- cube_genotype_info(parameters$cube)
  release_species <- resolve_release_species(parameters, releases)

  release_sex <- releases$releaseSex
  if (is.null(release_sex)) {
    release_sex <- "M"
  }
  if (length(release_sex) != 1L) {
    stop("releaseSex must be one of 'M', 'F', or 'both'")
  }
  release_sex <- tolower(as.character(release_sex[[1]]))
  if (!(release_sex %in% c("m", "f", "both"))) {
    stop("releaseSex must be one of 'M', 'F', or 'both'")
  }

  release_genotype <- resolve_release_genotype(parameters, releases, cube_info)

  spacing <- max(1L, releases_interval)
  event_timesteps <- releases_start + seq.int(0L, releases_number - 1L) * spacing

  rows <- vector("list", length(event_timesteps) * 2L)
  r_i <- 0L
  for (tt in event_timesteps) {
    if (release_sex == "both") {
      female_n <- floor(release_count / 2)
      male_n <- release_count - female_n
      if (female_n > 0L) {
        r_i <- r_i + 1L
        rows[[r_i]] <- list(
          timestep = as.integer(tt),
          species = release_species,
          sex = "F",
          genotype = release_genotype,
          count = as.integer(female_n)
        )
      }
      if (male_n > 0L) {
        r_i <- r_i + 1L
        rows[[r_i]] <- list(
          timestep = as.integer(tt),
          species = release_species,
          sex = "M",
          genotype = release_genotype,
          count = as.integer(male_n)
        )
      }
    } else {
      r_i <- r_i + 1L
      rows[[r_i]] <- list(
        timestep = as.integer(tt),
        species = release_species,
        sex = toupper(release_sex),
        genotype = release_genotype,
        count = as.integer(release_count)
      )
    }
  }

  rows <- rows[seq_len(r_i)]
  if (length(rows) == 0L) {
    return(empty_mosquito_release_schedule())
  }

  out <- do.call(
    rbind,
    lapply(rows, function(x) {
      data.frame(
        timestep = x$timestep,
        species = x$species,
        sex = x$sex,
        genotype = x$genotype,
        count = x$count,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out
}

ensure_mosquito_release_capacity <- function(
    variables,
    events,
    additional,
    default_species,
    default_genotype
  ) {
  if (additional <= 0L) {
    return(invisible(NULL))
  }

  variables$species$queue_extend(rep.int(default_species, additional))
  variables$mosquito_state$queue_extend(rep.int("NonExistent", additional))
  variables$geno_id$queue_extend(rep.int(as.integer(default_genotype), additional))
  events$mosquito_infection$queue_extend(additional)
  events$mosquito_death$queue_extend(additional)

  variables$species$.resize()
  variables$mosquito_state$.resize()
  variables$geno_id$.resize()
  events$mosquito_infection$.resize()
  events$mosquito_death$.resize()
}

create_mosquito_release_process <- function(
    solvers,
    models,
    variables,
    events,
    parameters
  ) {
  schedule <- parameters$releases_schedule
  if (is.null(schedule) || nrow(schedule) == 0L) {
    return(function(timestep) {})
  }
  if (is.null(parameters$cube)) {
    stop("Mosquito genotype releases require parameters$cube")
  }

  cube_info <- cube_genotype_info(parameters$cube)
  genotype_idx <- match(schedule$genotype, cube_info$genotypesID)
  if (any(is.na(genotype_idx))) {
    stop("Internal error: release schedule contains an unknown genotype")
  }

  species_lookup <- setNames(seq_along(parameters$species), parameters$species)
  if (!all(schedule$species %in% names(species_lookup))) {
    stop("Internal error: release schedule contains an unknown species")
  }

  rows_by_timestep <- split(seq_len(nrow(schedule)), schedule$timestep)

  if (native_mosquito_backend_enabled(parameters)) {
    return(function(timestep) {
      rows <- rows_by_timestep[[as.character(timestep)]]
      if (is.null(rows) || length(rows) == 0L) {
        return(invisible(NULL))
      }

      for (row_i in rows) {
        n_release <- as.integer(schedule$count[[row_i]])
        if (n_release <= 0L) {
          next
        }
        species_i <- species_lookup[[schedule$species[[row_i]]]]
        solvers[[species_i]]$apply_release(
          sex = schedule$sex[[row_i]],
          genotype_idx = genotype_idx[[row_i]],
          count = n_release,
          timestep = timestep
        )
      }
      invisible(NULL)
    })
  }

  if (!isTRUE(parameters$individual_mosquitoes)) {
    stop("Mosquito genotype releases are currently supported only when parameters$individual_mosquitoes = TRUE")
  }
  default_species <- as.character(parameters$species[[1]])
  default_genotype <- cube_wild_type_index(parameters$cube)
  emergence_rate <- .5 * 1 / parameters$dpl

  function(timestep) {
    rows <- rows_by_timestep[[as.character(timestep)]]
    if (is.null(rows) || length(rows) == 0L) {
      return(invisible(NULL))
    }

    debug_by_species <- list()
    if (genotype_debug_enabled(parameters, timestep)) {
      for (s_i in seq_along(parameters$species)) {
        debug_by_species[[s_i]] <- genotype_debug_species_counts(variables, models, parameters, s_i)
        genotype_debug_log_counts(
          parameters,
          timestep,
          "BEFORE_RELEASE",
          parameters$species[[s_i]],
          debug_by_species[[s_i]]
        )
      }
    }

    emergence_reserve <- 0L
    if (length(solvers) > 0) {
      p_counts <- vnapply(
        seq_along(solvers),
        function(i) {
          aquatic_stage_total(solvers[[i]]$get_states(), models[[i]]$cube, "P")
        }
      )
      emergence_reserve <- ceiling(sum(p_counts) * emergence_rate)
    }

    next_non_existent <- 1L
    queued_female_release <- FALSE
    for (row_i in rows) {
      n_release <- as.integer(schedule$count[[row_i]])
      if (n_release <= 0L) {
        next
      }

      species_name <- schedule$species[[row_i]]
      sex <- schedule$sex[[row_i]]
      g_idx <- genotype_idx[[row_i]]
      species_i <- species_lookup[[species_name]]

      if (sex == "M") {
        if (is.null(models[[species_i]]$genotype_state)) {
          stop("Male genotype releases require genotype tracking state (set parameters$cube in hybrid mode)")
        }
        male_counts <- models[[species_i]]$genotype_state$male_counts
        male_counts[[g_idx]] <- male_counts[[g_idx]] + n_release
        models[[species_i]]$genotype_state$male_counts <- male_counts
        genotype_debug_log(
          parameters,
          timestep,
          "RELEASE",
          species_name,
          sprintf("sex=M genotype=%s count=%d (male pool updated immediately)", schedule$genotype[[row_i]], n_release)
        )
        next
      }

      needed_index <- next_non_existent + n_release - 1L
      needed_capacity <- needed_index + emergence_reserve
      available <- variables$mosquito_state$get_size_of("NonExistent")
      if (available < needed_capacity) {
        ensure_mosquito_release_capacity(
          variables = variables,
          events = events,
          additional = needed_capacity - available,
          default_species = default_species,
          default_genotype = default_genotype
        )
      }

      non_existent <- variables$mosquito_state$get_index_of("NonExistent")
      target <- bitset_at(non_existent, seq.int(next_non_existent, needed_index))
      variables$geno_id$queue_update(g_idx, target)
      variables$mosquito_state$queue_update("Sm", target)
      variables$species$queue_update(species_name, target)
      next_non_existent <- needed_index + 1L
      queued_female_release <- TRUE
      genotype_debug_log(
        parameters,
        timestep,
        "RELEASE",
        species_name,
        sprintf("sex=F genotype=%s count=%d (female individuals queued+committed for same day)", schedule$genotype[[row_i]], n_release)
      )
    }

    if (queued_female_release) {
      variables$geno_id$.update()
      variables$mosquito_state$.update()
      variables$species$.update()
    }

    if (genotype_debug_enabled(parameters, timestep)) {
      for (s_i in seq_along(parameters$species)) {
        genotype_debug_log_counts(
          parameters,
          timestep,
          "AFTER_RELEASE",
          parameters$species[[s_i]],
          genotype_debug_species_counts(variables, models, parameters, s_i)
        )
      }
    }
  }
}
