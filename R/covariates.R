#' @title Create An Empty Node Contact Effect
#' @description Construct a no-op node contact effect specification.
#' @param note optional note recorded in the effect metadata.
#' @param label optional human-readable label.
#' @param source optional source/lineage metadata.
#' @return A validated node contact effect specification.
#' @export
node_contact_effect_none <- function(note = NULL, label = NULL, source = NULL) {
  list(
    type = "node_contact_effect",
    effect_type = "none",
    mechanism = "contact",
    normalization = "none",
    note = if (is.null(note)) NULL else as.character(note)[[1L]],
    label = if (is.null(label)) NULL else as.character(label)[[1L]],
    source = source
  )
}


#' @title Create A Binary Node Contact Effect
#' @description Construct a single-covariate binary contact effect specification.
#' @param covariate node-data column name.
#' @param reference_level reference level mapped to raw effect `1`.
#' @param exposed_level exposed level mapped to `exposed_multiplier`.
#' @param exposed_multiplier positive relative contact multiplier for the exposed
#'   level before any optional normalization.
#' @param normalization either `"mean_one"` or `"none"`.
#' @param label optional human-readable label.
#' @param source optional source/lineage metadata.
#' @return A validated node contact effect specification.
#' @export
node_contact_effect_binary <- function(
    covariate,
    reference_level,
    exposed_level,
    exposed_multiplier,
    normalization = c("mean_one", "none"),
    label = NULL,
    source = NULL
) {
  normalization <- .msimGD_contact_effect_normalization(normalization)
  covariate <- .msimGD_scalar_string(covariate, "covariate")
  reference_level <- .msimGD_effect_level(reference_level, "reference_level")
  exposed_level <- .msimGD_effect_level(exposed_level, "exposed_level")
  if (identical(reference_level, exposed_level)) {
    stop("reference_level and exposed_level must differ.", call. = FALSE)
  }
  exposed_multiplier <- .msimGD_positive_scalar(
    exposed_multiplier,
    "exposed_multiplier"
  )

  list(
    type = "node_contact_effect",
    effect_type = "binary",
    mechanism = "contact",
    covariate = covariate,
    reference_level = reference_level,
    exposed_level = exposed_level,
    exposed_multiplier = exposed_multiplier,
    normalization = normalization,
    label = if (is.null(label)) NULL else as.character(label)[[1L]],
    source = source
  )
}


#' @title Create A Categorical Node Contact Effect
#' @description Construct a single-covariate categorical contact effect
#'   specification using an explicit level-to-multiplier lookup.
#' @param covariate node-data column name.
#' @param level_multipliers named positive numeric vector mapping observed levels
#'   to raw relative contact effects before any optional normalization.
#' @param normalization either `"mean_one"` or `"none"`.
#' @param label optional human-readable label.
#' @param source optional source/lineage metadata.
#' @return A validated node contact effect specification.
#' @export
node_contact_effect_categorical <- function(
    covariate,
    level_multipliers,
    normalization = c("mean_one", "none"),
    label = NULL,
    source = NULL
) {
  normalization <- .msimGD_contact_effect_normalization(normalization)
  covariate <- .msimGD_scalar_string(covariate, "covariate")
  level_multipliers <- .msimGD_named_positive_vector(
    level_multipliers,
    "level_multipliers"
  )

  list(
    type = "node_contact_effect",
    effect_type = "categorical",
    mechanism = "contact",
    covariate = covariate,
    level_multipliers = level_multipliers,
    normalization = normalization,
    label = if (is.null(label)) NULL else as.character(label)[[1L]],
    source = source
  )
}


#' @title Create A Numeric Ratio Node Contact Effect
#' @description Construct a numeric contact effect specification where each
#'   `unit` increase relative to `reference_value` multiplies contact by
#'   `multiplier_per_unit`.
#' @param covariate node-data column name.
#' @param reference_value numeric reference value mapped to raw effect `1`.
#' @param multiplier_per_unit positive relative multiplier for each `unit`
#'   increase in the covariate.
#' @param unit positive numeric unit step used in the ratio rule.
#' @param normalization either `"mean_one"` or `"none"`.
#' @param multiplier_bounds optional positive numeric vector `c(lower, upper)` to
#'   clip raw effects after applying the ratio rule.
#' @param label optional human-readable label.
#' @param source optional source/lineage metadata.
#' @return A validated node contact effect specification.
#' @export
node_contact_effect_numeric_ratio <- function(
    covariate,
    reference_value,
    multiplier_per_unit,
    unit = 1,
    normalization = c("mean_one", "none"),
    multiplier_bounds = NULL,
    label = NULL,
    source = NULL
) {
  normalization <- .msimGD_contact_effect_normalization(normalization)
  covariate <- .msimGD_scalar_string(covariate, "covariate")
  reference_value <- .msimGD_finite_scalar(reference_value, "reference_value")
  multiplier_per_unit <- .msimGD_positive_scalar(
    multiplier_per_unit,
    "multiplier_per_unit"
  )
  unit <- .msimGD_positive_scalar(unit, "unit")
  multiplier_bounds <- .msimGD_contact_multiplier_bounds(multiplier_bounds)

  list(
    type = "node_contact_effect",
    effect_type = "numeric_ratio",
    mechanism = "contact",
    covariate = covariate,
    reference_value = reference_value,
    multiplier_per_unit = multiplier_per_unit,
    unit = unit,
    normalization = normalization,
    multiplier_bounds = multiplier_bounds,
    label = if (is.null(label)) NULL else as.character(label)[[1L]],
    source = source
  )
}


#' @title Create A Numeric Binned Node Contact Effect
#' @description Construct a numeric contact effect specification using explicit
#'   bins and per-bin positive multipliers.
#' @param covariate node-data column name.
#' @param breaks strictly increasing numeric breakpoints passed to `cut()`.
#' @param bin_multipliers positive numeric vector of length `length(breaks) - 1`.
#' @param normalization either `"mean_one"` or `"none"`.
#' @param include_lowest logical passed to `cut()`.
#' @param right logical passed to `cut()`.
#' @param label optional human-readable label.
#' @param source optional source/lineage metadata.
#' @return A validated node contact effect specification.
#' @export
node_contact_effect_numeric_bins <- function(
    covariate,
    breaks,
    bin_multipliers,
    normalization = c("mean_one", "none"),
    include_lowest = TRUE,
    right = TRUE,
    label = NULL,
    source = NULL
) {
  normalization <- .msimGD_contact_effect_normalization(normalization)
  covariate <- .msimGD_scalar_string(covariate, "covariate")
  breaks <- as.numeric(breaks)
  if (length(breaks) < 2L || anyNA(breaks) ||
      any(!(is.finite(breaks) | is.infinite(breaks))) ||
      any(diff(breaks) <= 0)) {
    stop(
      "breaks must be a strictly increasing numeric vector of length at least 2.",
      call. = FALSE
    )
  }
  bin_multipliers <- as.numeric(bin_multipliers)
  if (length(bin_multipliers) != length(breaks) - 1L ||
      any(!is.finite(bin_multipliers)) || any(bin_multipliers <= 0)) {
    stop(
      "bin_multipliers must be a positive numeric vector of length length(breaks) - 1.",
      call. = FALSE
    )
  }

  list(
    type = "node_contact_effect",
    effect_type = "numeric_bins",
    mechanism = "contact",
    covariate = covariate,
    breaks = breaks,
    bin_multipliers = bin_multipliers,
    normalization = normalization,
    include_lowest = isTRUE(include_lowest),
    right = isTRUE(right),
    label = if (is.null(label)) NULL else as.character(label)[[1L]],
    source = source
  )
}


#' @title Create A Multi-Covariate Node Contact Effect Bundle
#' @description Construct a multi-covariate contact effect specification by
#'   combining typed single-covariate contact effects with a product combiner.
#'   This is only valid when the component effects are mechanistically
#'   separable positive modifiers of the same contact process. This is not a
#'   claim of statistical independence of the observed covariates.
#' @param effects list of typed single-covariate node contact effects created by
#'   `node_contact_effect_*()`. Nested bundles are not supported in this first
#'   version.
#' @param normalization either `"mean_one"` or `"none"`, applied once after
#'   multiplying the unnormalized component effects.
#' @param label optional human-readable label.
#' @param source optional source/lineage metadata.
#' @return A validated node contact effect bundle specification.
#' @export
node_contact_effect_bundle <- function(
    effects,
    normalization = c("mean_one", "none"),
    label = NULL,
    source = NULL
) {
  normalization <- .msimGD_contact_effect_normalization(normalization)
  if (!is.list(effects) || length(effects) < 1L) {
    stop("effects must be a non-empty list of typed node contact effects.", call. = FALSE)
  }

  component_effects <- lapply(seq_along(effects), function(i) {
    effect <- effects[[i]]
    if (is.null(effect)) {
      stop(sprintf("effects[[%d]] is NULL.", i), call. = FALSE)
    }
    effect <- .msimGD_normalise_single_node_contact_effect(effect)
    if (!identical(effect$mechanism, "contact")) {
      stop(
        sprintf(
          "effects[[%d]] does not declare mechanism = \"contact\".",
          i
        ),
        call. = FALSE
      )
    }
    if (!identical(effect$normalization, "none")) {
      stop(
        paste(
          "Bundle components must use normalization = \"none\" because",
          "bundle normalization is applied once after the product combiner."
        ),
        call. = FALSE
      )
    }
    if (identical(effect$effect_type, "none")) {
      return(NULL)
    }
    effect
  })

  component_effects <- Filter(Negate(is.null), component_effects)
  if (length(component_effects) < 1L) {
    return(node_contact_effect_none(
      note = "bundle with only no-op component effects",
      label = label,
      source = source
    ))
  }

  covariates <- vapply(component_effects, `[[`, character(1), "covariate")
  duplicated_covariates <- unique(covariates[duplicated(covariates)])
  if (length(duplicated_covariates) > 0L) {
    stop(
      paste(
        "Duplicate covariate names are not allowed within one contact bundle:",
        paste(duplicated_covariates, collapse = ", "),
        "This first version assumes the listed covariates are mechanistically",
        "separable contact modifiers; obvious duplicates are rejected, but",
        "semantic overlap cannot be detected automatically."
      ),
      call. = FALSE
    )
  }

  list(
    type = "node_contact_effect",
    effect_type = "bundle",
    mechanism = "contact",
    combiner = "product",
    combination_assumption = .msimGD_contact_bundle_assumption_text(),
    covariates = covariates,
    component_effects = component_effects,
    normalization = normalization,
    label = if (is.null(label)) NULL else as.character(label)[[1L]],
    source = source
  )
}


#' @title Resolve A Node Contact Effect
#' @description Resolve a typed single-covariate contact effect specification,
#'   or a strict multi-covariate contact bundle, against node data into
#'   explicit node contact multipliers.
#' @param node_data data frame with one row per node.
#' @param effect node contact effect created by one of the
#'   `node_contact_effect_*()` constructors, or by
#'   `node_contact_effect_bundle()`.
#' @param node_index_col name of the positive integer node-index column.
#' @return A resolved node contact effect object suitable for
#'   `apply_node_contact_effect()`.
#' @export
resolve_node_contact_effect <- function(
    node_data,
    effect,
    node_index_col = "node_index"
) {
  effect <- .msimGD_normalise_node_contact_effect(effect)
  node_data <- as.data.frame(node_data, stringsAsFactors = FALSE)
  if (nrow(node_data) < 1L) {
    stop("node_data must contain at least one row.", call. = FALSE)
  }

  node_index_col <- .msimGD_scalar_string(node_index_col, "node_index_col")
  if (!(node_index_col %in% names(node_data))) {
    stop(sprintf("node_data is missing `%s`.", node_index_col), call. = FALSE)
  }
  node_index <- as.integer(node_data[[node_index_col]])
  if (length(node_index) != nrow(node_data) || anyNA(node_index) ||
      any(node_index <= 0L) || anyDuplicated(node_index)) {
    stop(
      paste(
        "node_data[[node_index_col]] must be a positive integer vector with one",
        "unique entry per node."
      ),
      call. = FALSE
    )
  }

  order_idx <- order(node_index)
  node_index <- node_index[order_idx]
  node_data <- node_data[order_idx, , drop = FALSE]
  rownames(node_data) <- NULL

  if (identical(effect$effect_type, "none")) {
    return(normalise_node_contact_surface(list(
      type = "none",
      note = effect$note
    )))
  }

  if (identical(effect$effect_type, "bundle")) {
    return(.msimGD_resolve_contact_effect_bundle(
      node_data = node_data,
      node_index = node_index,
      effect = effect
    ))
  }

  covariate <- effect$covariate
  if (!(covariate %in% names(node_data))) {
    stop(sprintf("node_data is missing covariate column `%s`.", covariate), call. = FALSE)
  }

  resolved <- .msimGD_resolve_single_contact_effect(
    node_data = node_data,
    node_index = node_index,
    effect = effect
  )
  contact_multiplier <- .msimGD_apply_contact_normalization(
    resolved$raw_contact_effect,
    effect$normalization
  )

  resolved_table <- cbind(
    data.frame(node_index = node_index, stringsAsFactors = FALSE),
    resolved$node_covariates,
    data.frame(
      raw_contact_effect = as.numeric(resolved$raw_contact_effect),
      contact_multiplier = as.numeric(contact_multiplier),
      stringsAsFactors = FALSE
    )
  )

  normalise_node_contact_surface(list(
    type = "contact_surface",
    hook = "human_blood_meal_rate",
    label = effect$label,
    source = effect$source,
    normalize = identical(effect$normalization, "mean_one"),
    node_index = node_index,
    node_covariates = resolved_table,
    raw_contact_effect = stats::setNames(
      as.numeric(resolved$raw_contact_effect),
      as.character(node_index)
    ),
    contact_multiplier = stats::setNames(
      as.numeric(contact_multiplier),
      as.character(node_index)
    ),
    effect_spec = effect
  ))
}


#' @title Normalise A Node Contact Surface
#' @description Validate and normalise a resolved node contact surface object.
#'
#' @param contact_surface NULL, a numeric node vector of contact multipliers, or a
#'   contact-surface list.
#' @param n_nodes optional expected number of nodes.
#'
#' @return A validated contact-surface specification list.
#' @export
normalise_node_contact_surface <- function(contact_surface, n_nodes = NULL) {
  if (is.null(contact_surface)) {
    return(list(
      type = "none",
      note = NULL
    ))
  }

  if (is.numeric(contact_surface)) {
    values <- as.numeric(contact_surface)
    if (length(values) < 1L || any(!is.finite(values)) || any(values <= 0)) {
      stop("Numeric contact_surface must be finite and positive.", call. = FALSE)
    }
    if (!is.null(n_nodes) && length(values) != n_nodes) {
      stop(sprintf("Numeric contact_surface must have length %d.", n_nodes), call. = FALSE)
    }

    node_index <- seq_along(values)
    return(list(
      type = "contact_surface",
      hook = "human_blood_meal_rate",
      label = NULL,
      source = NULL,
      normalize = FALSE,
      node_index = node_index,
      node_covariates = data.frame(
        node_index = node_index,
        raw_contact_effect = values,
        contact_multiplier = values,
        stringsAsFactors = FALSE
      ),
      raw_contact_effect = stats::setNames(values, as.character(node_index)),
      contact_multiplier = stats::setNames(values, as.character(node_index)),
      effect_spec = list(
        type = "node_contact_effect",
        effect_type = "manual_resolved",
        normalization = "none"
      )
    ))
  }

  if (!is.list(contact_surface) || is.null(contact_surface$type)) {
    stop("contact_surface must be NULL, numeric, or a contact-surface config list.", call. = FALSE)
  }

  type <- match.arg(as.character(contact_surface$type), c("none", "contact_surface"))
  if (identical(type, "none")) {
    return(list(
      type = "none",
      note = contact_surface$note
    ))
  }

  raw_values <- contact_surface$contact_multiplier
  values <- as.numeric(raw_values)
  if (length(values) < 1L || any(!is.finite(values)) || any(values <= 0)) {
    stop(
      "contact_surface$contact_multiplier must be finite and positive.",
      call. = FALSE
    )
  }

  node_index <- contact_surface$node_index
  if (is.null(node_index)) {
    if (!is.null(names(raw_values)) && all(nzchar(names(raw_values)))) {
      node_index <- as.integer(names(raw_values))
    } else {
      node_index <- seq_along(values)
    }
  }
  node_index <- as.integer(node_index)
  if (length(node_index) != length(values) ||
      anyNA(node_index) || any(node_index <= 0L) || anyDuplicated(node_index)) {
    stop(
      paste(
        "contact_surface node_index must be positive integers with one entry",
        "per contact_multiplier value and no duplicates."
      ),
      call. = FALSE
    )
  }
  if (!is.null(n_nodes) && length(values) != n_nodes) {
    stop(sprintf("contact_surface must define exactly %d nodes.", n_nodes), call. = FALSE)
  }

  order_idx <- order(node_index)
  node_index <- node_index[order_idx]
  values <- values[order_idx]
  names(values) <- as.character(node_index)

  raw_contact_effect <- contact_surface$raw_contact_effect
  if (is.null(raw_contact_effect)) {
    raw_contact_effect <- values
  } else {
    raw_contact_effect <- as.numeric(raw_contact_effect)
    if (length(raw_contact_effect) != length(values) ||
        any(!is.finite(raw_contact_effect)) || any(raw_contact_effect <= 0)) {
      stop("raw_contact_effect must be finite, positive, and match contact_multiplier length.", call. = FALSE)
    }
    raw_contact_effect <- raw_contact_effect[order_idx]
    names(raw_contact_effect) <- as.character(node_index)
  }

  resolved_table <- contact_surface$node_covariates
  if (is.null(resolved_table)) {
    resolved_table <- data.frame(
      node_index = node_index,
      raw_contact_effect = raw_contact_effect,
      contact_multiplier = values,
      stringsAsFactors = FALSE
    )
  } else {
    resolved_table <- as.data.frame(resolved_table, stringsAsFactors = FALSE)
    if (!("node_index" %in% names(resolved_table))) {
      resolved_table$node_index <- node_index
    }
    resolved_table <- resolved_table[
      match(node_index, as.integer(resolved_table$node_index)),
      ,
      drop = FALSE
    ]
    if (!("raw_contact_effect" %in% names(resolved_table))) {
      resolved_table$raw_contact_effect <- as.numeric(raw_contact_effect)
    }
    if (!("contact_multiplier" %in% names(resolved_table))) {
      resolved_table$contact_multiplier <- values
    }
  }

  effect_spec <- contact_surface$effect_spec
  if (is.null(effect_spec)) {
    legacy_fields <- list(
      link = contact_surface$link,
      intercept = contact_surface$intercept,
      covariate_names = contact_surface$covariate_names,
      coefficients = contact_surface$coefficients,
      centers = contact_surface$centers,
      scales = contact_surface$scales
    )
    legacy_fields <- legacy_fields[!vapply(legacy_fields, is.null, logical(1))]
    if (length(legacy_fields) > 0L) {
      effect_spec <- c(
        list(type = "legacy_contact_effect"),
        legacy_fields
      )
    }
  }

  list(
    type = "contact_surface",
    hook = as.character(if (!is.null(contact_surface$hook)) contact_surface$hook else "human_blood_meal_rate")[[1L]],
    label = if (is.null(contact_surface$label)) NULL else as.character(contact_surface$label)[[1L]],
    source = contact_surface$source,
    normalize = isTRUE(contact_surface$normalize),
    node_index = node_index,
    node_covariates = resolved_table,
    raw_contact_effect = raw_contact_effect,
    contact_multiplier = values,
    effect_spec = effect_spec
  )
}


#' @title Apply A Node Contact Surface
#' @description Apply one node's resolved contact multiplier and metadata to a
#' model parameter list.
#'
#' @param parameters model parameter list for a single node.
#' @param contact_surface resolved contact surface from
#'   `normalise_node_contact_surface()`.
#' @param node_index positive integer node id.
#'
#' @return Updated parameter list.
#' @export
apply_node_contact_surface <- function(parameters, contact_surface, node_index = 1L) {
  surface <- normalise_node_contact_surface(contact_surface)
  if (identical(surface$type, "none")) {
    return(parameters)
  }

  node_index <- as.integer(node_index)
  if (length(node_index) != 1L || is.na(node_index) || node_index <= 0L) {
    stop("node_index must be a single integer > 0.", call. = FALSE)
  }

  node_name <- as.character(node_index)
  match_idx <- match(node_name, names(surface$contact_multiplier))
  if (is.na(match_idx)) {
    stop(sprintf("contact_surface has no entry for node %d.", node_index), call. = FALSE)
  }

  # Backward-compatible surface field: scalar contact_multiplier is folded into
  # human_slot_contact_multiplier when human variables are initialized.
  parameters$contact_multiplier <- as.numeric(surface$contact_multiplier[[match_idx]])
  parameters$contact_multiplier_hook <- surface$hook
  parameters$contact_multiplier_label <- surface$label
  parameters$contact_multiplier_source <- surface$source
  parameters$contact_multiplier_normalize <- isTRUE(surface$normalize)

  effect_spec <- surface$effect_spec
  if (is.null(effect_spec)) {
    effect_spec <- list(
      type = "node_contact_effect",
      effect_type = "manual_resolved",
      normalization = if (isTRUE(surface$normalize)) "mean_one" else "none"
    )
  }
  parameters$contact_multiplier_effect_spec <- effect_spec

  node_row <- surface$node_covariates[
    match(node_index, as.integer(surface$node_covariates$node_index)),
    ,
    drop = FALSE
  ]
  if (nrow(node_row) == 1L) {
    drop_cols <- c(
      "node_index",
      "contact_multiplier",
      "raw_contact_effect"
    )
    keep <- !(names(node_row) %in% drop_cols)
    if (any(keep)) {
      node_values <- lapply(node_row[1, keep, drop = FALSE], function(x) {
        if (length(x) == 1L) {
          return(x[[1L]])
        }
        x
      })
      parameters$contact_multiplier_covariates <- node_values
    }
  }

  parameters
}


#' @title Apply A Resolved Node Contact Effect
#' @description Apply one node's resolved contact effect to a model parameter
#'   list.
#' @param parameters model parameter list for a single node.
#' @param resolved_effect resolved effect from `resolve_node_contact_effect()`.
#' @param node_index positive integer node id.
#' @return Updated parameter list.
#' @export
apply_node_contact_effect <- function(parameters, resolved_effect, node_index = 1L) {
  apply_node_contact_surface(
    parameters = parameters,
    contact_surface = resolved_effect,
    node_index = node_index
  )
}


.msimGD_scalar_string <- function(x, name) {
  x <- as.character(x)
  if (length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("%s must be a single non-empty string.", name), call. = FALSE)
  }
  x[[1L]]
}


.msimGD_effect_level <- function(x, name) {
  if (length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be a single non-missing level.", name), call. = FALSE)
  }
  as.character(x)[[1L]]
}


.msimGD_finite_scalar <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x)) {
    stop(sprintf("%s must be a single finite number.", name), call. = FALSE)
  }
  as.numeric(x)
}


.msimGD_positive_scalar <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("%s must be a single finite number > 0.", name), call. = FALSE)
  }
  as.numeric(x)
}


.msimGD_named_positive_vector <- function(x, name) {
  x_names <- names(x)
  x <- as.numeric(x)
  if (length(x) < 1L || is.null(x_names) || anyNA(x_names) || any(!nzchar(x_names)) ||
      anyDuplicated(x_names) || any(!is.finite(x)) || any(x <= 0)) {
    stop(
      sprintf(
        "%s must be a named positive numeric vector with unique non-empty names.",
        name
      ),
      call. = FALSE
    )
  }
  stats::setNames(as.numeric(x), as.character(x_names))
}


.msimGD_contact_effect_normalization <- function(normalization) {
  match.arg(as.character(normalization), c("mean_one", "none"))
}


.msimGD_contact_multiplier_bounds <- function(multiplier_bounds) {
  if (is.null(multiplier_bounds)) {
    return(NULL)
  }
  multiplier_bounds <- as.numeric(multiplier_bounds)
  if (length(multiplier_bounds) != 2L ||
      any(!is.finite(multiplier_bounds)) ||
      multiplier_bounds[[1L]] <= 0 ||
      multiplier_bounds[[2L]] < multiplier_bounds[[1L]]) {
    stop(
      "multiplier_bounds must be NULL or a positive numeric vector c(lower, upper) with lower <= upper.",
      call. = FALSE
    )
  }
  multiplier_bounds
}


.msimGD_contact_bundle_assumption_text <- function() {
  paste(
    "Product combination is valid only when component effects are",
    "mechanistically separable positive modifiers of the same contact",
    "process. This is not a claim of statistical independence of the",
    "observed covariates."
  )
}


.msimGD_normalise_node_contact_effect <- function(effect) {
  if (is.null(effect)) {
    return(node_contact_effect_none())
  }
  if (!is.list(effect) || !identical(effect$type, "node_contact_effect")) {
    stop("effect must be a node contact effect created by node_contact_effect_*().", call. = FALSE)
  }

  effect_type <- match.arg(
    as.character(effect$effect_type),
    c("none", "binary", "categorical", "numeric_ratio", "numeric_bins", "bundle")
  )

  if (identical(effect_type, "bundle")) {
    return(node_contact_effect_bundle(
      effects = effect$component_effects,
      normalization = effect$normalization,
      label = effect$label,
      source = effect$source
    ))
  }

  .msimGD_normalise_single_node_contact_effect(effect)
}


.msimGD_normalise_single_node_contact_effect <- function(effect) {
  if (is.null(effect)) {
    return(node_contact_effect_none())
  }
  if (!is.list(effect) || !identical(effect$type, "node_contact_effect")) {
    stop("effect must be a node contact effect created by node_contact_effect_*().", call. = FALSE)
  }

  effect_type <- match.arg(
    as.character(effect$effect_type),
    c("none", "binary", "categorical", "numeric_ratio", "numeric_bins")
  )

  switch(
    effect_type,
    none = node_contact_effect_none(
      note = effect$note,
      label = effect$label,
      source = effect$source
    ),
    binary = node_contact_effect_binary(
      covariate = effect$covariate,
      reference_level = effect$reference_level,
      exposed_level = effect$exposed_level,
      exposed_multiplier = effect$exposed_multiplier,
      normalization = effect$normalization,
      label = effect$label,
      source = effect$source
    ),
    categorical = node_contact_effect_categorical(
      covariate = effect$covariate,
      level_multipliers = effect$level_multipliers,
      normalization = effect$normalization,
      label = effect$label,
      source = effect$source
    ),
    numeric_ratio = node_contact_effect_numeric_ratio(
      covariate = effect$covariate,
      reference_value = effect$reference_value,
      multiplier_per_unit = effect$multiplier_per_unit,
      unit = effect$unit,
      normalization = effect$normalization,
      multiplier_bounds = effect$multiplier_bounds,
      label = effect$label,
      source = effect$source
    ),
    numeric_bins = node_contact_effect_numeric_bins(
      covariate = effect$covariate,
      breaks = effect$breaks,
      bin_multipliers = effect$bin_multipliers,
      normalization = effect$normalization,
      include_lowest = effect$include_lowest,
      right = effect$right,
      label = effect$label,
      source = effect$source
    )
  )
}


.msimGD_apply_contact_normalization <- function(raw_contact_effect, normalization) {
  normalization <- .msimGD_contact_effect_normalization(normalization)
  raw_contact_effect <- as.numeric(raw_contact_effect)
  if (length(raw_contact_effect) < 1L ||
      any(!is.finite(raw_contact_effect)) ||
      any(raw_contact_effect <= 0)) {
    stop("raw_contact_effect must be finite and positive.", call. = FALSE)
  }

  if (identical(normalization, "none")) {
    return(raw_contact_effect)
  }

  raw_contact_effect / mean(raw_contact_effect)
}


.msimGD_resolve_single_contact_effect <- function(node_data, node_index, effect) {
  covariate <- effect$covariate

  switch(
    effect$effect_type,
    binary = {
      values <- node_data[[covariate]]
      if (anyNA(values)) {
        stop(sprintf("Binary covariate `%s` contains missing values.", covariate), call. = FALSE)
      }
      values_chr <- as.character(values)
      allowed <- c(effect$reference_level, effect$exposed_level)
      bad_levels <- setdiff(unique(values_chr), allowed)
      if (length(bad_levels) > 0L) {
        stop(
          sprintf(
            "Binary covariate `%s` contains unexpected level(s): %s",
            covariate,
            paste(bad_levels, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      raw_effect <- ifelse(
        values_chr == effect$reference_level,
        1,
        effect$exposed_multiplier
      )
      list(
        node_covariates = data.frame(
          setNames(list(values_chr), covariate),
          stringsAsFactors = FALSE
        ),
        raw_contact_effect = as.numeric(raw_effect)
      )
    },
    categorical = {
      values <- node_data[[covariate]]
      if (anyNA(values)) {
        stop(sprintf("Categorical covariate `%s` contains missing values.", covariate), call. = FALSE)
      }
      values_chr <- as.character(values)
      bad_levels <- setdiff(unique(values_chr), names(effect$level_multipliers))
      if (length(bad_levels) > 0L) {
        stop(
          sprintf(
            "Categorical covariate `%s` contains unconfigured level(s): %s",
            covariate,
            paste(bad_levels, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      raw_effect <- unname(effect$level_multipliers[values_chr])
      list(
        node_covariates = data.frame(
          setNames(list(values_chr), covariate),
          stringsAsFactors = FALSE
        ),
        raw_contact_effect = as.numeric(raw_effect)
      )
    },
    numeric_ratio = {
      values <- as.numeric(node_data[[covariate]])
      if (length(values) != nrow(node_data) || any(!is.finite(values))) {
        stop(sprintf("Numeric covariate `%s` must be finite numeric.", covariate), call. = FALSE)
      }
      scaled_distance <- (values - effect$reference_value) / effect$unit
      raw_effect <- exp(log(effect$multiplier_per_unit) * scaled_distance)
      if (!is.null(effect$multiplier_bounds)) {
        raw_effect <- pmax(raw_effect, effect$multiplier_bounds[[1L]])
        raw_effect <- pmin(raw_effect, effect$multiplier_bounds[[2L]])
      }
      list(
        node_covariates = data.frame(
          setNames(list(values), covariate),
          scaled_distance = scaled_distance,
          stringsAsFactors = FALSE
        ),
        raw_contact_effect = as.numeric(raw_effect)
      )
    },
    numeric_bins = {
      values <- as.numeric(node_data[[covariate]])
      if (length(values) != nrow(node_data) || any(!is.finite(values))) {
        stop(sprintf("Numeric covariate `%s` must be finite numeric.", covariate), call. = FALSE)
      }
      bin_index <- cut(
        values,
        breaks = effect$breaks,
        labels = FALSE,
        include_lowest = isTRUE(effect$include_lowest),
        right = isTRUE(effect$right)
      )
      if (anyNA(bin_index)) {
        stop(
          sprintf(
            "Numeric covariate `%s` could not be assigned to the supplied bins.",
            covariate
          ),
          call. = FALSE
        )
      }
      bin_label <- as.character(cut(
        values,
        breaks = effect$breaks,
        include_lowest = isTRUE(effect$include_lowest),
        right = isTRUE(effect$right)
      ))
      raw_effect <- effect$bin_multipliers[bin_index]
      list(
        node_covariates = data.frame(
          setNames(list(values), covariate),
          bin_index = as.integer(bin_index),
          bin_label = bin_label,
          stringsAsFactors = FALSE
        ),
        raw_contact_effect = as.numeric(raw_effect)
      )
    },
    stop(sprintf("Unsupported contact effect type `%s`.", effect$effect_type), call. = FALSE)
  )
}


.msimGD_resolve_contact_effect_bundle <- function(node_data, node_index, effect) {
  component_effects <- effect$component_effects
  component_resolved <- lapply(component_effects, function(component) {
    if (!(component$covariate %in% names(node_data))) {
      stop(
        sprintf(
          "node_data is missing covariate column `%s` required by the contact bundle.",
          component$covariate
        ),
        call. = FALSE
      )
    }
    .msimGD_resolve_single_contact_effect(
      node_data = node_data,
      node_index = node_index,
      effect = component
    )
  })

  raw_component_matrix <- vapply(
    component_resolved,
    function(component) as.numeric(component$raw_contact_effect),
    numeric(length(node_index))
  )
  if (is.null(dim(raw_component_matrix))) {
    raw_component_matrix <- matrix(
      raw_component_matrix,
      ncol = 1L,
      dimnames = list(NULL, component_effects[[1L]]$covariate)
    )
  }

  # Product combination is only valid under mechanistic separability:
  # each component effect is treated as a distinct positive modifier of the
  # same contact process. This is not a claim of statistical independence of
  # the observed covariates.
  raw_contact_effect <- apply(raw_component_matrix, 1L, prod)
  contact_multiplier <- .msimGD_apply_contact_normalization(
    raw_contact_effect,
    effect$normalization
  )

  resolved_table <- data.frame(node_index = node_index, stringsAsFactors = FALSE)
  for (i in seq_along(component_effects)) {
    component <- component_effects[[i]]
    component_covariates <- component_resolved[[i]]$node_covariates
    covariate <- component$covariate

    resolved_table[[covariate]] <- component_covariates[[covariate]]
    extra_cols <- setdiff(names(component_covariates), covariate)
    if (length(extra_cols) > 0L) {
      for (extra_col in extra_cols) {
        resolved_table[[paste0(covariate, "__", extra_col)]] <- component_covariates[[extra_col]]
      }
    }
    resolved_table[[paste0(covariate, "__raw_contact_effect")]] <- as.numeric(
      component_resolved[[i]]$raw_contact_effect
    )
  }
  resolved_table$raw_contact_effect <- as.numeric(raw_contact_effect)
  resolved_table$contact_multiplier <- as.numeric(contact_multiplier)

  normalise_node_contact_surface(list(
    type = "contact_surface",
    hook = "human_blood_meal_rate",
    label = effect$label,
    source = effect$source,
    normalize = identical(effect$normalization, "mean_one"),
    node_index = node_index,
    node_covariates = resolved_table,
    raw_contact_effect = stats::setNames(
      as.numeric(raw_contact_effect),
      as.character(node_index)
    ),
    contact_multiplier = stats::setNames(
      as.numeric(contact_multiplier),
      as.character(node_index)
    ),
    effect_spec = effect
  ))
}
