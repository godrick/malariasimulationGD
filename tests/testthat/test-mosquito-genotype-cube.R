test_that("hybrid mosquito genotype tracking is backward compatible and consistent", {
  skip_if_not_installed("MGDrivE")

  base_params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 50,
    total_M = 200,
    init_foim = 0,
    progress_bar = FALSE
  ))
  base_params <- malariasimulationGD::parameterise_total_M(base_params, base_params$total_M)

  timesteps <- 20

  set.seed(123)
  baseline <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = base_params)
  baseline_df <- baseline$data
  baseline_adults <- baseline_df$Sm_gamb_count + baseline_df$Pm_gamb_count + baseline_df$Im_gamb_count

  cube1 <- MGDrivE::cubeMendelian(gtype = c("AA"))
  params_cube1 <- base_params
  params_cube1$cube <- cube1

  set.seed(123)
  run_cube1 <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params_cube1)
  cube1_df <- run_cube1$data
  aq1_E <- attr(cube1_df, "mosquito_aquatic_genotype_E")
  aq1_L <- attr(cube1_df, "mosquito_aquatic_genotype_L")
  aq1_P <- attr(cube1_df, "mosquito_aquatic_genotype_P")

  expect_identical(names(cube1_df), names(baseline_df))
  for (nm in names(baseline_df)) {
    expect_identical(cube1_df[[nm]], baseline_df[[nm]], info = paste("column", nm))
  }
  expect_equal(
    cube1_df$Sm_gamb_count + cube1_df$Pm_gamb_count + cube1_df$Im_gamb_count,
    baseline_adults,
    tolerance = 0
  )
  expect_false(is.null(aq1_E))
  expect_false(is.null(aq1_L))
  expect_false(is.null(aq1_P))
  expect_identical(dim(aq1_E), c(as.integer(timesteps), 1L))
  expect_equal(drop(aq1_E[, 1]), cube1_df$E_gamb_count, tolerance = 0)
  expect_equal(drop(aq1_L[, 1]), cube1_df$L_gamb_count, tolerance = 0)
  expect_equal(drop(aq1_P[, 1]), cube1_df$P_gamb_count, tolerance = 0)

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  params_cube3 <- base_params
  params_cube3$cube <- cube3

  set.seed(123)
  run_cube3 <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params_cube3)
  geno <- run_cube3$mosquito_genotypes
  cube3_df <- run_cube3$data
  aq3_E <- attr(cube3_df, "mosquito_aquatic_genotype_E")
  aq3_L <- attr(cube3_df, "mosquito_aquatic_genotype_L")
  aq3_P <- attr(cube3_df, "mosquito_aquatic_genotype_P")

  expect_false(is.null(geno))
  expect_identical(colnames(geno$female), cube3$genotypesID)
  expect_identical(colnames(geno$male), cube3$genotypesID)
  expect_equal(geno$V[, 1], rep(1, nrow(geno$V)), tolerance = 0)

  female_total <- rowSums(geno$female)
  male_total <- rowSums(geno$male)
  expect_equal(female_total + male_total, geno$total_adults, tolerance = 0)
  expect_equal(
    female_total,
    cube3_df$Sm_gamb_count + cube3_df$Pm_gamb_count + cube3_df$Im_gamb_count,
    tolerance = 0
  )
  expect_false(is.null(aq3_E))
  expect_false(is.null(aq3_L))
  expect_false(is.null(aq3_P))
  expect_equal(rowSums(aq3_E), cube3_df$E_gamb_count, tolerance = 0)
  expect_equal(rowSums(aq3_L), cube3_df$L_gamb_count, tolerance = 0)
  expect_equal(rowSums(aq3_P), cube3_df$P_gamb_count, tolerance = 0)
  expect_true(all(aq3_E[, c("Aa", "aa"), drop = FALSE] == 0))
  expect_true(all(aq3_L[, c("Aa", "aa"), drop = FALSE] == 0))
  expect_true(all(aq3_P[, c("Aa", "aa"), drop = FALSE] == 0))
  expect_true(all(geno$female[, c("Aa", "aa"), drop = FALSE] == 0))
  expect_true(all(geno$male[, c("Aa", "aa"), drop = FALSE] == 0))
  expect_null(attr(cube3_df, "mosquito_release_schedule"))
})

test_that("scheduled male genotype release has delayed adult genotype effects with genotype-resolved aquatic states", {
  skip_if_not_installed("MGDrivE")

  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 50,
    total_M = 200,
    init_foim = 0,
    progress_bar = FALSE
  ))
  params <- malariasimulationGD::parameterise_total_M(params, params$total_M)

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  cube3$releaseType <- "aa"
  params$cube <- cube3
  params <- malariasimulationGD::set_releases(params, list(
    releasesStart = 20,
    releasesNumber = 1,
    releaseCount = 200,
    releaseSex = "M"
  ))

  set.seed(321)
  out <- malariasimulationGD::run_resumable_simulation(50, parameters = params)
  geno <- out$mosquito_genotypes
  rel <- attr(out$data, "mosquito_release_schedule")
  aqE <- attr(out$data, "mosquito_aquatic_genotype_E")
  aqL <- attr(out$data, "mosquito_aquatic_genotype_L")
  aqP <- attr(out$data, "mosquito_aquatic_genotype_P")

  expect_false(is.null(rel))
  expect_identical(
    rel,
    data.frame(
      timestep = 20L,
      species = "gamb",
      sex = "M",
      genotype = "aa",
      count = 200L,
      stringsAsFactors = FALSE
    )
  )
  expect_equal(out$data$n_released_gamb[20], 200)
  expect_equal(unname(geno$male[20, "aa"]), 200)
  expect_equal(unname(geno$male[20, "Aa"]), 0)
  expect_true(all(geno$female[20, c("Aa", "aa"), drop = FALSE] == 0))

  # Non-WT aquatic stages can appear immediately via mating (after the day-20 ODE
  # step), but there should be no same-day adult Aa creation from the old p_g
  # emergence shortcut.
  expect_true(any(aqE[20:nrow(aqE), c("Aa", "aa"), drop = FALSE] > 0))

  adult_Aa <- geno$female[, "Aa"] + geno$male[, "Aa"]
  first_Aa_day <- which(adult_Aa > 0)[1]
  expect_false(is.na(first_Aa_day))
  expect_gt(first_Aa_day, 20)
})

test_that("set_releases validates release genotype against cube genotypes", {
  skip_if_not_installed("MGDrivE")

  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    progress_bar = FALSE
  ))
  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  params$cube <- cube3

  expect_error(
    malariasimulationGD::set_releases(params, list(
      releasesStart = 5,
      releaseCount = 10,
      releaseGenotype = "ZZ"
    )),
    "releaseGenotype 'ZZ' is not in cube\\$genotypesID"
  )
})

test_that("female releases can expand mosquito capacity beyond mosquito_limit", {
  skip_if_not_installed("MGDrivE")

  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 20,
    total_M = 5,
    init_foim = 0,
    mosquito_limit = 20,
    progress_bar = FALSE
  ))
  params <- malariasimulationGD::parameterise_total_M(params, params$total_M)

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  cube3$releaseType <- "aa"
  params$cube <- cube3
  params <- malariasimulationGD::set_releases(params, list(
    releasesStart = 2,
    releaseCount = 200,
    releaseSex = "F"
  ))

  set.seed(456)
  out <- malariasimulationGD::run_resumable_simulation(4, parameters = params)

  expect_equal(out$data$n_released_gamb[2], 200)
  expect_equal(unname(out$mosquito_genotypes$female[2, "aa"]), 200)
  expect_gt(out$data$Sm_gamb_count[2] + out$data$Pm_gamb_count[2] + out$data$Im_gamb_count[2], 150)
})

test_that("sterile tau (all zero) yields zero egg input and aquatic collapse", {
  skip_if_not_installed("MGDrivE")

  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 50,
    total_M = 200,
    init_foim = 0,
    progress_bar = FALSE
  ))
  params <- malariasimulationGD::parameterise_total_M(params, params$total_M)

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  cube3$tau <- array(0, dim = dim(cube3$ih))
  params$cube <- cube3

  set.seed(789)
  out <- malariasimulationGD::run_resumable_simulation(80, parameters = params)
  aqE <- attr(out$data, "mosquito_aquatic_genotype_E")
  aqL <- attr(out$data, "mosquito_aquatic_genotype_L")
  aqP <- attr(out$data, "mosquito_aquatic_genotype_P")

  aquatic_total <- rowSums(aqE) + rowSums(aqL) + rowSums(aqP)
  expect_true(all(aquatic_total >= 0))
  expect_lt(tail(aquatic_total, 1), aquatic_total[1])
  expect_lt(tail(aquatic_total, 1) / aquatic_total[1], 0.01)
})

test_that("vector_infectivity_g is ignored with warning when cube is NULL", {
  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 40,
    total_M = 120,
    init_foim = 0,
    progress_bar = FALSE
  ))
  params <- malariasimulationGD::parameterise_total_M(params, params$total_M)

  timesteps <- 15

  set.seed(1201)
  baseline <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params)

  params_vi <- params
  params_vi$vector_infectivity_g <- c(AA = 0.5)

  expect_warning(
    {
      set.seed(1201)
      ignored <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params_vi)
    },
    "vector_infectivity_g ignored because .*cube=NULL"
  )

  expect_identical(ignored$data, baseline$data)
})

test_that("WT-only vector_infectivity_g=1 is a strict no-op", {
  skip_if_not_installed("MGDrivE")

  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 50,
    total_M = 200,
    init_foim = 0,
    progress_bar = FALSE
  ))
  params <- malariasimulationGD::parameterise_total_M(params, params$total_M)
  params$cube <- MGDrivE::cubeMendelian(gtype = c("AA"))

  params_vi <- params
  params_vi$vector_infectivity_g <- c(AA = 1)

  set.seed(1202)
  base <- malariasimulationGD::run_resumable_simulation(20, parameters = params)
  set.seed(1202)
  weighted <- malariasimulationGD::run_resumable_simulation(20, parameters = params_vi)

  expect_identical(weighted$data, base$data)
  expect_false("infectivity_weighted_I_gamb" %in% names(weighted$data))
  expect_false("vector_infectivity_mean_gamb" %in% names(weighted$data))
  expect_null(attr(weighted$data, "mosquito_infectious_genotype_counts"))
})

test_that("vector_infectivity_g reduces human transmission when low-infectivity genotype rises among infectious females", {
  skip_if_not_installed("MGDrivE")

  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 80,
    total_M = 300,
    init_foim = 0,
    progress_bar = FALSE
  ))
  params <- malariasimulationGD::parameterise_total_M(params, params$total_M)

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  cube3$releaseType <- "aa"
  params$cube <- cube3
  params <- malariasimulationGD::set_releases(params, list(
    releasesStart = 10,
    releasesNumber = 10,
    releasesInterval = 7,
    releaseCount = 200,
    releaseSex = "M"
  ))

  params_all_one <- params
  params_all_one$vector_infectivity_g <- c(AA = 1, Aa = 1, aa = 1)

  params_aa_zero <- params
  params_aa_zero$vector_infectivity_g <- list(
    gamb = c(AA = 1, Aa = 1, aa = 0)
  )

  timesteps <- 160

  set.seed(1203)
  out_all_one <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params_all_one)
  set.seed(1203)
  out_aa_zero <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params_aa_zero)

  df_a <- out_all_one$data
  df_b <- out_aa_zero$data

  expect_true(all(c(
    "infectivity_weighted_I_gamb",
    "vector_infectivity_mean_gamb"
  ) %in% names(df_a)))
  expect_true(all(c(
    "infectivity_weighted_I_gamb",
    "vector_infectivity_mean_gamb"
  ) %in% names(df_b)))

  expect_equal(df_a$infectivity_weighted_I_gamb, df_a$Im_gamb_count, tolerance = 0)
  infected_days_a <- df_a$Im_gamb_count > 0
  expect_true(all(df_a$vector_infectivity_mean_gamb[infected_days_a] == 1))

  tail_idx <- seq.int(max(1L, nrow(df_b) - 29L), nrow(df_b))
  expect_true(any(!is.na(df_b$vector_infectivity_mean_gamb[tail_idx])))
  expect_lt(
    mean(df_b$vector_infectivity_mean_gamb[tail_idx], na.rm = TRUE),
    mean(df_a$vector_infectivity_mean_gamb[tail_idx], na.rm = TRUE)
  )
  expect_lt(
    mean(df_b$infectivity_weighted_I_gamb[tail_idx], na.rm = TRUE),
    mean(df_b$Im_gamb_count[tail_idx], na.rm = TRUE)
  )
  expect_lt(
    sum(df_b$n_infections[tail_idx], na.rm = TRUE),
    sum(df_a$n_infections[tail_idx], na.rm = TRUE)
  )
})

test_that("vector_infectivity_g validates names, allows relative multipliers, and checks cube$b ranges", {
  skip_if_not_installed("MGDrivE")

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))

  params_missing <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    progress_bar = FALSE
  ))
  params_missing$cube <- cube3
  params_missing$vector_infectivity_g <- c(AA = 1, Aa = 1)
  expect_error(
    malariasimulationGD::run_resumable_simulation(1, parameters = params_missing),
    "missing genotype names.*aa"
  )

  params_range <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    progress_bar = FALSE
  ))
  params_range$cube <- cube3
  params_range$vector_infectivity_g <- c(AA = 1, Aa = 1.2, aa = 0)
  expect_no_error(
    malariasimulationGD::run_resumable_simulation(1, parameters = params_range)
  )

  params_cube_range <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    progress_bar = FALSE
  ))
  params_cube_range$cube <- cube3
  params_cube_range$cube$b <- c(AA = 1, Aa = 1.2, aa = 0)
  expect_error(
    malariasimulationGD::run_resumable_simulation(1, parameters = params_cube_range),
    "parameters\\$cube\\$b values must be in \\[0, 1\\]"
  )

  params_compartmental <- malariasimulationGD::get_parameters(list(
    individual_mosquitoes = FALSE,
    native_mosquito_backend = FALSE,
    progress_bar = FALSE
  ))
  params_compartmental$cube <- cube3
  params_compartmental$vector_infectivity_g <- c(AA = 1, Aa = 1, aa = 0)
  expect_error(
    malariasimulationGD::run_resumable_simulation(1, parameters = params_compartmental),
    "requires genotype-resolved mosquitoes"
  )
})

test_that("cube$omega all ones is a strict no-op (including WT-only)", {
  skip_if_not_installed("MGDrivE")

  base_params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 50,
    total_M = 200,
    init_foim = 0,
    progress_bar = FALSE
  ))
  base_params <- malariasimulationGD::parameterise_total_M(base_params, base_params$total_M)

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  params_cube3 <- base_params
  params_cube3$cube <- cube3
  params_cube3_omega <- params_cube3
  params_cube3_omega$cube$omega <- c(AA = 1, Aa = 1, aa = 1)

  set.seed(1301)
  out_cube3 <- malariasimulationGD::run_resumable_simulation(30, parameters = params_cube3)
  set.seed(1301)
  out_cube3_omega <- malariasimulationGD::run_resumable_simulation(30, parameters = params_cube3_omega)
  expect_identical(out_cube3_omega$data, out_cube3$data)
  expect_identical(out_cube3_omega$mosquito_genotypes$female, out_cube3$mosquito_genotypes$female)
  expect_identical(out_cube3_omega$mosquito_genotypes$male, out_cube3$mosquito_genotypes$male)

  cube1 <- MGDrivE::cubeMendelian(gtype = c("AA"))
  params_cube1 <- base_params
  params_cube1$cube <- cube1
  params_cube1_omega <- params_cube1
  params_cube1_omega$cube$omega <- c(AA = 1)

  set.seed(1302)
  out_cube1 <- malariasimulationGD::run_resumable_simulation(30, parameters = params_cube1)
  set.seed(1302)
  out_cube1_omega <- malariasimulationGD::run_resumable_simulation(30, parameters = params_cube1_omega)
  expect_identical(out_cube1_omega$data, out_cube1$data)
})

test_that("cube$omega > 1 for aa reduces aa adult share under repeated aa releases", {
  skip_if_not_installed("MGDrivE")

  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    human_population = 80,
    total_M = 300,
    init_foim = 0,
    progress_bar = FALSE
  ))
  params <- malariasimulationGD::parameterise_total_M(params, params$total_M)

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  cube3$releaseType <- "aa"
  params$cube <- cube3
  params <- malariasimulationGD::set_releases(params, list(
    releasesStart = 10,
    releasesNumber = 10,
    releasesInterval = 7,
    releaseCount = 200,
    releaseSex = "M"
  ))

  params_omega <- params
  params_omega$cube$omega <- c(AA = 1, Aa = 1, aa = 2)

  timesteps <- 160

  set.seed(1303)
  out_base <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params)
  set.seed(1303)
  out_omega <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params_omega)

  total_adults_base <- out_base$mosquito_genotypes$total_adults
  total_adults_omega <- out_omega$mosquito_genotypes$total_adults
  aa_total_base <- out_base$mosquito_genotypes$female[, "aa"] + out_base$mosquito_genotypes$male[, "aa"]
  aa_total_omega <- out_omega$mosquito_genotypes$female[, "aa"] + out_omega$mosquito_genotypes$male[, "aa"]

  aa_frac_base <- ifelse(total_adults_base > 0, aa_total_base / total_adults_base, NA_real_)
  aa_frac_omega <- ifelse(total_adults_omega > 0, aa_total_omega / total_adults_omega, NA_real_)
  tail_idx <- seq.int(max(1L, timesteps - 29L), timesteps)

  expect_true(any(!is.na(aa_frac_base[tail_idx])))
  expect_true(any(!is.na(aa_frac_omega[tail_idx])))
  expect_lt(
    mean(aa_frac_omega[tail_idx], na.rm = TRUE),
    mean(aa_frac_base[tail_idx], na.rm = TRUE)
  )
})

test_that("cube$omega validates genotype names and values", {
  skip_if_not_installed("MGDrivE")

  params <- malariasimulationGD::get_parameters(list(
    native_mosquito_backend = FALSE,
    individual_mosquitoes = TRUE,
    progress_bar = FALSE
  ))
  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  params$cube <- cube3

  params_missing <- params
  params_missing$cube$omega <- c(AA = 1, Aa = 1)
  expect_error(
    malariasimulationGD::run_resumable_simulation(1, parameters = params_missing),
    "missing genotype names.*aa"
  )

  params_negative <- params
  params_negative$cube$omega <- c(AA = 1, Aa = -0.1, aa = 1)
  expect_error(
    malariasimulationGD::run_resumable_simulation(1, parameters = params_negative),
    "entries must be finite and >= 0"
  )
})
