four_node_test_root <- normalizePath(
  file.path(developer_test_repo(), "test-scripts/four_node_linear"),
  mustWork = FALSE
)
if (!file.exists(file.path(four_node_test_root, "config.R"))) {
  testthat::skip("Four-node developer example sources are not available.")
}
four_node_test_env <- new.env(parent = globalenv())
sys.source(file.path(four_node_test_root, "config.R"), four_node_test_env)
sys.source(file.path(four_node_test_root, "helpers.R"), four_node_test_env)

test_that("four-node burn-in and observed timing are validated", {
  e <- four_node_test_env
  cfg <- e$four_node_config()
  expect_equal(cfg$burnin_days, 1825L)
  expect_equal(cfg$days, 1825L)
  expect_equal(cfg$burnin_days + cfg$release_day, 2008L)
  expect_gte(cfg$population, 1000L)
  expect_silent(e$validate_four_node_config(cfg))
  for (value in c(-1, 1.5, NA_real_, Inf)) {
    bad <- cfg
    bad$burnin_days <- value
    expect_error(e$validate_four_node_config(bad), "burnin_days")
  }
  cfg$burnin_days <- 0L
  expect_silent(e$validate_four_node_config(cfg))
})

test_that("four-node site equivalence rejects contact differences", {
  e <- four_node_test_env
  cfg <- e$four_node_config()
  p <- list(human_population = cfg$population, contact_multiplier = 1,
            human_slot_contact_multiplier = 1, human_mobility_enabled = FALSE,
            model_seasonality = FALSE, carrying_capacity = FALSE)
  parameters <- stats::setNames(rep(list(p), 4L), cfg$sites)
  parameters$A$releases <- list(releasesStart = 2008L)
  expect_silent(e$four_node_assert_equivalent(parameters, cfg))
  parameters$B$contact_multiplier <- 1.1
  expect_error(e$four_node_assert_equivalent(parameters, cfg), "identical")
})

test_that("burn-in is excluded while release and genotype histories retain their alignment", {
  e <- four_node_test_env
  cfg <- e$four_node_config()
  cfg$burnin_days <- 4L
  cfg$days <- 4L
  cfg$release_day <- 2L
  cfg$incidence_window_days <- 3L
  data <- lapply(cfg$sites, function(site) {
    d <- data.frame(timestep = 1:8, n_age_182_3649 = 100,
                    n_inc_182_3649 = c(rep(100L, 4L), rep(1L, 4L)),
                    n_inc_clinical_182_3649 = 0)
    f <- m <- matrix(0, 8L, 2L, dimnames = list(NULL, c("WW", "HH")))
    f[, "WW"] <- m[, "WW"] <- 10
    if (site %in% cfg$release_sites) {
      m[6:8, "HH"] <- 5
      attr(d, "mosquito_release_schedule") <- data.frame(
        timestep = 6L, species = "gamb", sex = "M", genotype = "HH",
        count = cfg$release_count
      )
    }
    attr(d, "mosquito_genotype_counts_female") <- f
    attr(d, "mosquito_genotype_counts_male") <- m
    d
  })
  x <- e$four_node_extract(data, cfg)
  expect_equal(nrow(x$daily), 16L)
  expect_equal(nrow(x$burnin_daily), 16L)
  expect_equal(x$daily$timestep, rep(1:4, 4L))
  expect_equal(x$daily$simulation_day, rep(5:8, 4L))
  expect_true(all(x$burnin_carriers$gd_carriers == 0))
  expect_equal(x$schedule$timestep, c(2L, 2L))
  expect_equal(x$schedule$simulation_day, c(6L, 6L))
  expect_equal(x$carriers$gd_carriers[x$carriers$site == "A"], c(0, 5, 5, 5))
  windows <- e$four_node_incidence(x$daily, cfg)
  expect_equal(sum(windows$infections), 16L)
  expect_equal(sum(windows$person_days), 1600)
  expect_equal(windows$infection_incidence_ppy, rep(3.65, 8L))
  diagnostics <- e$four_node_baseline_diagnostics(x, cfg)
  expect_equal(diagnostics$days[diagnostics$period == "observed_before_release"], rep(1L, 4L))
  # A stray H carrier during burn-in must be rejected, even though hidden
  # from the observed plots.
  attr(data[[1L]], "mosquito_genotype_counts_male")[4L, "HH"] <- 1
  expect_error(e$four_node_extract(data, cfg))
})
