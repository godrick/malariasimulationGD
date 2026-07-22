test_that("standalone mosquito ODE returns lifecycle counts without infection stages", {
  parameters <- get_mosquito_parameters(list(
    total_M = 50,
    native_mosquito_nE = 2L,
    native_mosquito_nL = 2L,
    native_mosquito_nP = 1L
  ))

  out <- run_mosquito_simulation(3L, parameters)

  expect_s3_class(parameters, "mosquito_parameters")
  expect_equal(sort(unique(out$counts$timestep)), 0:3)
  expect_setequal(
    unique(out$counts$stage),
    c("egg", "larva", "pupa", "adult_male", "adult_female")
  )
  expect_false(any(c("S", "E", "I", "Pm", "Im", "Sm") %in% names(out)))
  expect_false(any(out$state_counts$stage %in% c("S", "E", "I")))

  reconciled <- aggregate(
    count ~ timestep + node + species + stage + genotype,
    data = out$state_counts,
    FUN = sum
  )
  expect_equal(out$counts[order(out$counts$timestep, out$counts$stage), ],
               reconciled[order(reconciled$timestep, reconciled$stage), ],
               ignore_attr = TRUE)
})

test_that("standalone releases can target each mosquito lifecycle stage", {
  parameters <- get_mosquito_parameters(list(
    total_M = 50,
    native_mosquito_nE = 2L,
    native_mosquito_nL = 2L,
    native_mosquito_nP = 2L
  ))
  empty_state <- data.frame(
    node = character(0), species = character(0), stage = character(0),
    genotype = character(0), count = numeric(0)
  )
  releases <- data.frame(
    timestep = 0L,
    node = "node1",
    species = "gamb",
    stage = c("egg", "larva", "pupa", "adult_male", "adult_female"),
    genotype = "WT",
    count = c(1, 2, 3, 4, 5),
    substage = c(2L, 2L, 2L, NA, NA),
    mate_genotype = c(NA, NA, NA, NA, "WT")
  )

  out <- run_mosquito_simulation(0L, parameters, releases = releases, initial_state = empty_state)

  got <- out$counts[out$counts$timestep == 0L, c("stage", "count")]
  got <- got[match(c("egg", "larva", "pupa", "adult_male", "adult_female"), got$stage), ]
  expect_equal(got$count, c(1, 2, 3, 4, 5))

  female_detail <- out$state_counts[
    out$state_counts$stage == "adult_female" & out$state_counts$count > 0,
  ]
  expect_equal(female_detail$mating_status, "mated")
  expect_equal(female_detail$mate_genotype, "WT")
})

test_that("standalone deterministic checkpoint resumes equivalently", {
  parameters <- get_mosquito_parameters(list(
    total_M = 50,
    native_mosquito_nE = 2L,
    native_mosquito_nL = 2L,
    native_mosquito_nP = 1L
  ))

  full <- run_mosquito_simulation(5L, parameters)
  first <- run_mosquito_simulation(2L, parameters)
  resumed <- run_mosquito_simulation(5L, parameters, initial_state = first$state)

  full_tail <- full$counts[full$counts$timestep >= 2L, ]
  rownames(full_tail) <- NULL
  rownames(resumed$counts) <- NULL
  expect_equal(resumed$counts, full_tail, tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("standalone tau sampler is integer, nonnegative, and checkpoint reproducible", {
  parameters <- get_mosquito_parameters(list(
    total_M = 50,
    native_mosquito_nE = 2L,
    native_mosquito_nL = 2L,
    native_mosquito_nP = 1L,
    mosquito_tau_step = 0.25
  ))

  full <- run_mosquito_simulation(5L, parameters, sampler = "tau", seed = 123)
  first <- run_mosquito_simulation(2L, parameters, sampler = "tau", seed = 123)
  resumed <- run_mosquito_simulation(5L, parameters, sampler = "tau", initial_state = first$state)

  expect_true(all(full$state_counts$count >= 0))
  expect_true(all(abs(full$state_counts$count - round(full$state_counts$count)) < 1e-8))

  full_tail <- full$counts[full$counts$timestep >= 2L, ]
  rownames(full_tail) <- NULL
  rownames(resumed$counts) <- NULL
  expect_identical(resumed$counts, full_tail)
})

test_that("standalone mosquito network accepts movement configuration", {
  p1 <- get_mosquito_parameters(list(total_M = 50))
  p2 <- get_mosquito_parameters(list(total_M = 25))
  movement <- matrix(c(0, 1, 0, 0), nrow = 2L, byrow = TRUE)
  p1$move_probs <- movement
  p2$move_probs <- movement
  p1$move_rates <- c(0.2, 0)
  p2$move_rates <- c(0.2, 0)

  out <- run_mosquito_simulation(1L, list(origin = p1, dest = p2))

  expect_setequal(unique(out$counts$node), c("origin", "dest"))
  expect_true(any(out$counts$count[out$counts$node == "dest" & out$counts$timestep == 1L] > 0))
})
