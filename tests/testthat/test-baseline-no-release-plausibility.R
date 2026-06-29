build_no_release_baseline_cube <- function() {
  genotypes <- c("WW", "HH")
  G <- length(genotypes)
  ih <- array(
    0,
    dim = c(G, G, G),
    dimnames = list(genotypes, genotypes, genotypes)
  )
  ih["WW", "WW", "WW"] <- 1
  ih["WW", "HH", "HH"] <- 1
  ih["HH", "WW", "HH"] <- 1
  ih["HH", "HH", "HH"] <- 1

  list(
    ih = ih,
    tau = array(1, dim = c(G, G, G), dimnames = list(genotypes, genotypes, genotypes)),
    eta = matrix(1, nrow = G, ncol = G, dimnames = list(genotypes, genotypes)),
    b = setNames(rep(1, G), genotypes),
    c = setNames(rep(1, G), genotypes),
    phi = setNames(rep(0.5, G), genotypes),
    omega = setNames(rep(1, G), genotypes),
    xiF = setNames(rep(1, G), genotypes),
    xiM = setNames(rep(1, G), genotypes),
    s = setNames(rep(1, G), genotypes),
    genotypesID = genotypes,
    wildType = "WW"
  )
}

build_no_release_baseline_parameters <- function(
    human_population = 1000L,
    init_EIR = 10
) {
  parameters <- get_parameters(list(
    human_population = human_population,
    native_mosquito_backend = TRUE,
    individual_mosquitoes = FALSE,
    model_seasonality = FALSE,
    human_mobility_enabled = FALSE,
    human_move_probs = NULL,
    human_move_rates = NULL,
    move_probs = matrix(1, 1, 1),
    move_rates = 0,
    bednets = FALSE,
    spraying = FALSE,
    progress_bar = FALSE,
    cube = build_no_release_baseline_cube()
  ))

  parameters$incidence_rendering_min_ages <- 0
  parameters$incidence_rendering_max_ages <- 100 * 365 - 1
  parameters$clinical_incidence_rendering_min_ages <- c(0, 5 * 365)
  parameters$clinical_incidence_rendering_max_ages <- c(5 * 365 - 1, 100 * 365 - 1)

  set_equilibrium(parameters, init_EIR = init_EIR, native_total_M = TRUE)
}

expect_no_negative_or_nonfinite <- function(data, columns) {
  present <- intersect(columns, names(data))
  expect_true(length(present) > 0)

  values <- as.matrix(data[present])
  expect_true(all(is.finite(values)))
  expect_true(all(values >= 0))
}

test_that("high-transmission native baseline remains wildtype-only without releases", {
  parameters <- build_no_release_baseline_parameters()

  set.seed(2)
  out <- run_resumable_simulation(3 * 365, parameters = parameters)
  data <- out$data
  final_year <- utils::tail(data, 365)

  female <- out$mosquito_genotypes$female
  male <- out$mosquito_genotypes$male
  aquatic_E <- attr(data, "mosquito_aquatic_genotype_E")
  aquatic_L <- attr(data, "mosquito_aquatic_genotype_L")
  aquatic_P <- attr(data, "mosquito_aquatic_genotype_P")

  expect_null(attr(data, "mosquito_release_schedule"))
  expect_false(any(grepl("^n_released_", names(data))))
  expect_identical(colnames(female), c("WW", "HH"))
  expect_identical(colnames(male), c("WW", "HH"))
  expect_equal(sum(female[, "HH"]), 0)
  expect_equal(sum(male[, "HH"]), 0)
  expect_equal(sum(aquatic_E[, "HH"]), 0)
  expect_equal(sum(aquatic_L[, "HH"]), 0)
  expect_equal(sum(aquatic_P[, "HH"]), 0)

  rendered_adult_females <- data$Sm_gamb_count + data$Pm_gamb_count + data$Im_gamb_count
  expect_equal(rowSums(female), rendered_adult_females, tolerance = 1e-8)
  expect_equal(female[, "WW"], rendered_adult_females, tolerance = 1e-8)
  expect_true(all(is.finite(rowSums(male))))
  expect_true(all(rowSums(male) >= 0))
  expect_equal(rowSums(male), rowSums(female), tolerance = 1e-8)

  human_total <- data$S_count + data$A_count + data$D_count + data$U_count + data$Tr_count
  expect_equal(human_total, rep(parameters$human_population, nrow(data)))

  expect_no_negative_or_nonfinite(
    data,
    c(
      "E_gamb_count", "L_gamb_count", "P_gamb_count",
      "Sm_gamb_count", "Pm_gamb_count", "Im_gamb_count", "total_M_gamb",
      "EIR_gamb", "FOIM_gamb", "mu_gamb", "infectivity",
      "n_infections", "n_bitten",
      "S_count", "A_count", "D_count", "U_count", "Tr_count",
      "n_age_730_3650", "n_detect_lm_730_3650", "n_detect_pcr_730_3650",
      "n_inc_0_36499", "n_inc_clinical_0_1824"
    )
  )

  annual_eir <- sum(final_year$EIR_gamb) / parameters$human_population
  pfpr_2_10_lm <- mean(
    final_year$n_detect_lm_730_3650 /
      pmax(final_year$n_age_730_3650, 1)
  )
  annual_infection_incidence <- sum(final_year$n_inc_0_36499) /
    parameters$human_population
  under5_clinical_incidence <- sum(final_year$n_inc_clinical_0_1824) /
    mean(final_year$n_age_0_1824)
  adult_female_cv <- stats::sd(final_year$total_M_gamb) /
    mean(final_year$total_M_gamb)

  expect_gt(annual_eir, 5)
  expect_lt(annual_eir, 20)
  expect_gt(pfpr_2_10_lm, 0.25)
  expect_lt(pfpr_2_10_lm, 0.85)
  expect_gt(annual_infection_incidence, 1)
  expect_lt(annual_infection_incidence, 10)
  expect_gt(under5_clinical_incidence, 0.2)
  expect_lt(under5_clinical_incidence, 4)
  expect_lt(adult_female_cv, 1e-4)
})
