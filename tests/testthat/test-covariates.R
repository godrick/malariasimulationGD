test_that("normalise_node_contact_surface accepts explicit node multipliers", {
  surface <- normalise_node_contact_surface(list(
    type = "contact_surface",
    hook = "human_blood_meal_rate",
    label = "test_surface",
    node_index = c(3L, 1L, 2L),
    contact_multiplier = c(`3` = 0.7, `1` = 1.0, `2` = 1.3)
  ))

  expect_identical(surface$type, "contact_surface")
  expect_identical(surface$hook, "human_blood_meal_rate")
  expect_identical(surface$label, "test_surface")
  expect_equal(surface$node_index, c(1L, 2L, 3L))
  expect_equal(unname(surface$contact_multiplier), c(1.0, 1.3, 0.7), tolerance = 1e-12)
  expect_true(all(surface$contact_multiplier > 0))
})


test_that("binary node contact effect resolves to mean-one multipliers", {
  effect <- node_contact_effect_binary(
    covariate = "open_eaves",
    reference_level = "no",
    exposed_level = "yes",
    exposed_multiplier = 1.2,
    normalization = "mean_one",
    label = "open_eaves_effect"
  )

  resolved <- resolve_node_contact_effect(
    node_data = data.frame(
      node_index = 1:3,
      open_eaves = c("no", "yes", "yes"),
      stringsAsFactors = FALSE
    ),
    effect = effect
  )

  raw <- c(1, 1.2, 1.2)
  expected <- raw / mean(raw)

  expect_identical(resolved$effect_spec$effect_type, "binary")
  expect_equal(unname(resolved$raw_contact_effect), raw, tolerance = 1e-12)
  expect_equal(unname(resolved$contact_multiplier), expected, tolerance = 1e-12)
  expect_equal(mean(unname(resolved$contact_multiplier)), 1, tolerance = 1e-12)
})


test_that("categorical node contact effect resolves by level lookup", {
  effect <- node_contact_effect_categorical(
    covariate = "roof_type",
    level_multipliers = c(metal = 0.95, tile = 1.0, thatch = 1.1),
    normalization = "none"
  )

  resolved <- resolve_node_contact_effect(
    node_data = data.frame(
      node_index = 1:3,
      roof_type = c("metal", "tile", "thatch"),
      stringsAsFactors = FALSE
    ),
    effect = effect
  )

  expect_identical(resolved$effect_spec$effect_type, "categorical")
  expect_equal(unname(resolved$contact_multiplier), c(0.95, 1.0, 1.1), tolerance = 1e-12)
  expect_identical(
    as.character(resolved$node_covariates$roof_type),
    c("metal", "tile", "thatch")
  )
})


test_that("numeric ratio effect handles negative values and preserves positivity", {
  effect <- node_contact_effect_numeric_ratio(
    covariate = "window_score",
    reference_value = 0,
    multiplier_per_unit = 1.2,
    unit = 1,
    normalization = "none"
  )

  resolved <- resolve_node_contact_effect(
    node_data = data.frame(
      node_index = 1:3,
      window_score = c(-1, 0, 1)
    ),
    effect = effect
  )

  expect_equal(
    unname(resolved$contact_multiplier),
    c(1 / 1.2, 1, 1.2),
    tolerance = 1e-12
  )
  expect_true(all(resolved$contact_multiplier > 0))
})


test_that("numeric binned effect resolves bin labels and multipliers", {
  effect <- node_contact_effect_numeric_bins(
    covariate = "windows_count",
    breaks = c(-Inf, 1, 3, Inf),
    bin_multipliers = c(0.9, 1.0, 1.15),
    normalization = "none"
  )

  resolved <- resolve_node_contact_effect(
    node_data = data.frame(
      node_index = 1:4,
      windows_count = c(0, 1, 2, 5)
    ),
    effect = effect
  )

  expect_equal(
    unname(resolved$contact_multiplier),
    c(0.9, 0.9, 1.0, 1.15),
    tolerance = 1e-12
  )
  expect_true(all(c("bin_index", "bin_label") %in% names(resolved$node_covariates)))
})


test_that("multi-covariate contact bundle multiplies component effects and normalizes once", {
  effect <- node_contact_effect_bundle(
    effects = list(
      node_contact_effect_binary(
        covariate = "open_eaves",
        reference_level = "no",
        exposed_level = "yes",
        exposed_multiplier = 1.2,
        normalization = "none"
      ),
      node_contact_effect_categorical(
        covariate = "roof_type",
        level_multipliers = c(metal = 0.95, tile = 1.0, thatch = 1.1),
        normalization = "none"
      ),
      node_contact_effect_numeric_ratio(
        covariate = "windows_score",
        reference_value = 0,
        multiplier_per_unit = 1.1,
        normalization = "none"
      )
    ),
    normalization = "mean_one",
    label = "housing_bundle",
    source = "unit_test"
  )

  resolved <- resolve_node_contact_effect(
    node_data = data.frame(
      node_index = 1:3,
      open_eaves = c("no", "yes", "yes"),
      roof_type = c("metal", "tile", "thatch"),
      windows_score = c(-1, 0, 1),
      stringsAsFactors = FALSE
    ),
    effect = effect
  )

  raw <- c(1.0, 1.2, 1.2) * c(0.95, 1.0, 1.1) * c(1 / 1.1, 1.0, 1.1)
  expected <- raw / mean(raw)

  expect_identical(resolved$effect_spec$effect_type, "bundle")
  expect_identical(resolved$effect_spec$combiner, "product")
  expect_match(
    resolved$effect_spec$combination_assumption,
    "mechanistically separable",
    fixed = TRUE
  )
  expect_match(
    resolved$effect_spec$combination_assumption,
    "not a claim of statistical independence",
    fixed = TRUE
  )
  expect_equal(unname(resolved$raw_contact_effect), raw, tolerance = 1e-12)
  expect_equal(unname(resolved$contact_multiplier), expected, tolerance = 1e-12)
  expect_equal(mean(unname(resolved$contact_multiplier)), 1, tolerance = 1e-12)
  expect_true(all(c(
    "open_eaves__raw_contact_effect",
    "roof_type__raw_contact_effect",
    "windows_score__scaled_distance",
    "windows_score__raw_contact_effect"
  ) %in% names(resolved$node_covariates)))
})


test_that("apply_node_contact_effect injects multiplier and typed metadata", {
  effect <- node_contact_effect_binary(
    covariate = "open_eaves",
    reference_level = "no",
    exposed_level = "yes",
    exposed_multiplier = 1.2,
    normalization = "mean_one",
    label = "open_eaves_effect",
    source = "unit_test"
  )
  resolved <- resolve_node_contact_effect(
    node_data = data.frame(
      node_index = 1:2,
      open_eaves = c("no", "yes"),
      stringsAsFactors = FALSE
    ),
    effect = effect
  )

  params <- apply_node_contact_effect(
    parameters = get_parameters(),
    resolved_effect = resolved,
    node_index = 2L
  )

  expect_equal(params$contact_multiplier, unname(resolved$contact_multiplier[["2"]]), tolerance = 1e-12)
  expect_identical(params$contact_multiplier_hook, "human_blood_meal_rate")
  expect_identical(params$contact_multiplier_label, "open_eaves_effect")
  expect_identical(params$contact_multiplier_source, "unit_test")
  expect_identical(params$contact_multiplier_effect_spec$effect_type, "binary")
  expect_true(is.list(params$contact_multiplier_covariates))
  expect_identical(params$contact_multiplier_covariates$open_eaves, "yes")
})


test_that("bundle contact effect injects combined typed metadata", {
  effect <- node_contact_effect_bundle(
    effects = list(
      node_contact_effect_binary(
        covariate = "open_eaves",
        reference_level = "no",
        exposed_level = "yes",
        exposed_multiplier = 1.2,
        normalization = "none"
      ),
      node_contact_effect_numeric_bins(
        covariate = "windows_count",
        breaks = c(-Inf, 1, Inf),
        bin_multipliers = c(0.9, 1.1),
        normalization = "none"
      )
    ),
    normalization = "mean_one",
    label = "bundle_effect",
    source = "unit_test"
  )
  resolved <- resolve_node_contact_effect(
    node_data = data.frame(
      node_index = 1:2,
      open_eaves = c("no", "yes"),
      windows_count = c(0, 3)
    ),
    effect = effect
  )

  params <- apply_node_contact_effect(
    parameters = get_parameters(),
    resolved_effect = resolved,
    node_index = 2L
  )

  expect_identical(params$contact_multiplier_effect_spec$effect_type, "bundle")
  expect_identical(params$contact_multiplier_effect_spec$combiner, "product")
  expect_identical(params$contact_multiplier_covariates$open_eaves, "yes")
  expect_identical(params$contact_multiplier_covariates$windows_count, 3)
  expect_equal(
    params$contact_multiplier_covariates$open_eaves__raw_contact_effect,
    1.2,
    tolerance = 1e-12
  )
  expect_equal(
    params$contact_multiplier_covariates$windows_count__raw_contact_effect,
    1.1,
    tolerance = 1e-12
  )
})


test_that("scalar legacy contact multiplier routes through human slot contact", {
  p_base <- get_parameters(list(
    human_population = 10,
    enable_heterogeneity = FALSE,
    contact_multiplier = 1
  ))
  p_scaled <- get_parameters(list(
    human_population = 10,
    enable_heterogeneity = FALSE,
    contact_multiplier = 2.5
  ))

  set.seed(101)
  v_base <- create_variables(p_base)
  set.seed(101)
  v_scaled <- create_variables(p_scaled)

  expect_equal(runtime_contact_multiplier(p_scaled, 1L), 1, tolerance = 1e-12)
  expect_equal(
    v_scaled$human_slot_contact_multiplier$get_values(),
    rep(2.5, p_scaled$human_population),
    tolerance = 1e-12
  )
  expect_equal(
    human_blood_meal_rate(1L, v_scaled, p_scaled, 0) /
      human_blood_meal_rate(1L, v_base, p_base, 0),
    2.5,
    tolerance = 1e-12
  )
  expect_equal(
    equilibrium_species_traits(p_base, 1L)$a,
    equilibrium_species_traits(p_scaled, 1L)$a,
    tolerance = 1e-12
  )
})


test_that("human slot contact multipliers scale contact and redistribute biting", {
  zeta <- c(1, 1)
  psi <- c(1, 1)

  expect_equal(human_pi(zeta, psi), c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(human_pi(zeta, psi, c(2, 2)), c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(human_pi(zeta, psi, c(1, 3)), c(0.25, 0.75), tolerance = 1e-12)
  expect_equal(human_slot_contact_rate_multiplier(zeta, psi), 1, tolerance = 1e-12)
  expect_equal(human_slot_contact_rate_multiplier(zeta, psi, c(2, 2)), 2, tolerance = 1e-12)
  expect_equal(human_slot_contact_rate_multiplier(zeta, psi, c(1, 3)), 2, tolerance = 1e-12)
  expect_equal(
    human_slot_contact_rate_multiplier(c(1, 3), psi, c(1, 3)),
    2.5,
    tolerance = 1e-12
  )
})


test_that("human slot contact multiplier initializes as a human variable", {
  parameters <- get_parameters(list(
    human_population = 3,
    human_slot_contact_multiplier = c(1, 2, 3)
  ))
  variables <- create_variables(parameters)

  expect_equal(
    variables$human_slot_contact_multiplier$get_values(),
    c(1, 2, 3),
    tolerance = 1e-12
  )
  expect_error(
    create_variables(get_parameters(list(
      human_population = 3,
      human_slot_contact_multiplier = c(1, 2)
    ))),
    "length 1 or match"
  )
})


test_that("legacy scalar contact multiplier combines with explicit slot contact", {
  parameters <- get_parameters(list(
    human_population = 3,
    contact_multiplier = 2,
    human_slot_contact_multiplier = c(1, 2, 3)
  ))
  variables <- create_variables(parameters)

  expect_equal(
    variables$human_slot_contact_multiplier$get_values(),
    c(2, 4, 6),
    tolerance = 1e-12
  )
  expect_equal(runtime_contact_multiplier(parameters, 1L), 1, tolerance = 1e-12)
})


test_that("legacy contact multiplier rejects vector inputs", {
  parameters <- get_parameters(list(human_population = 3))
  parameters$contact_multiplier <- c(gamb = 0.8, fun = 1.2)

  expect_error(
    create_variables(parameters),
    "single positive finite number",
    fixed = TRUE
  )
})


test_that("contact-effect helpers reject invalid inputs", {
  expect_error(
    node_contact_effect_binary(
      covariate = "open_eaves",
      reference_level = "no",
      exposed_level = "yes",
      exposed_multiplier = 0
    ),
    "> 0"
  )
  expect_error(
    resolve_node_contact_effect(
      node_data = data.frame(node_index = 1:2, open_eaves = c("no", "maybe")),
      effect = node_contact_effect_binary(
        covariate = "open_eaves",
        reference_level = "no",
        exposed_level = "yes",
        exposed_multiplier = 1.2
      )
    ),
    "unexpected level"
  )
  expect_error(
    node_contact_effect_bundle(
      effects = list(
        node_contact_effect_binary(
          covariate = "open_eaves",
          reference_level = "no",
          exposed_level = "yes",
          exposed_multiplier = 1.2,
          normalization = "mean_one"
        ),
        node_contact_effect_categorical(
          covariate = "roof_type",
          level_multipliers = c(metal = 1, tile = 1.1),
          normalization = "none"
        )
      )
    ),
    "normalization = \"none\""
  )
  expect_error(
    node_contact_effect_bundle(
      effects = list(
        node_contact_effect_binary(
          covariate = "open_eaves",
          reference_level = "no",
          exposed_level = "yes",
          exposed_multiplier = 1.2,
          normalization = "none"
        ),
        node_contact_effect_categorical(
          covariate = "open_eaves",
          level_multipliers = c(no = 1, yes = 1.1),
          normalization = "none"
        )
      )
    ),
    "Duplicate covariate names"
  )
  expect_error(
    normalise_node_contact_surface(c(1, -1)),
    "finite and positive"
  )
  expect_error(
    runtime_contact_multiplier(get_parameters(list(contact_multiplier = 0)), 1L),
    "positive finite"
  )
})
