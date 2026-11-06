test_that("the optional cache-key fallback is stable and does not draw random numbers", {
  parameters <- get_parameters(list(human_population = 8L, human_initialization_seed = 123L,
                                    human_initialization_burnin_timesteps = 2L))
  cache_key <- stationary_human_initializer_cache_key
  mockery::stub(cache_key, "requireNamespace", function(...) FALSE)
  set.seed(20261105)
  rng_before <- .Random.seed

  key <- cache_key(parameters)
  expect_identical(.Random.seed, rng_before)
  expect_match(key, "^[[:xdigit:]]{32}$")
  expect_identical(cache_key(parameters), key)
  parameters$human_initialization_seed <- 124L
  expect_false(identical(cache_key(parameters), key))
})

test_that("digest cache keys are stable when the optional dependency is installed", {
  skip_if_not_installed("digest")
  parameters <- get_parameters(list(human_population = 8L, human_initialization_seed = 123L,
                                    human_initialization_burnin_timesteps = 2L))
  set.seed(20261105)
  rng_before <- .Random.seed

  key <- stationary_human_initializer_cache_key(parameters)
  expect_identical(.Random.seed, rng_before)
  expect_match(key, "^[[:xdigit:]]{16}$")
  expect_identical(stationary_human_initializer_cache_key(parameters), key)
  parameters$human_initialization_seed <- 124L
  expect_false(identical(stationary_human_initializer_cache_key(parameters), key))
})
