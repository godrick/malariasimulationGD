create_transmission_mixer <- function(
  variables,
  parameters,
  lagged_eir,
  lagged_infectivity,
  mixing_tt,
  export_mixing,
  import_mixing,
  p_captured_tt,
  p_captured,
  p_success,
  lagged_transmission_eir = lagged_eir
  ) {
  function (timestep) {
    n_pops <- length(variables)
    mix_t <- match_timestep(mixing_tt, timestep)
    p_mix_export <- export_mixing[[mix_t]]
    p_mix_import <- import_mixing[[mix_t]]
    p_captured_t <- p_captured[[match_timestep(p_captured_tt, timestep)]]
    n_species <- length(parameters[[1]]$species)
    species_eir_raw <- vapply(
      seq_along(lagged_eir),
      function(i) {
        vnapply(
          lagged_eir[[i]],
          function(l) l$get(timestep - parameters[[i]]$de)
        )
      },
      numeric(n_species)
    )
    # vapply returns n_species x n_pops when n_species > 1, or a plain
    # vector of length n_pops when n_species == 1.  We need an
    # n_pops x n_species matrix in both cases.
    if (is.matrix(species_eir_raw)) {
      species_eir <- t(species_eir_raw)
    } else {
      species_eir <- matrix(species_eir_raw, ncol = 1L)
    }
    species_transmission_eir_raw <- vapply(
      seq_along(lagged_transmission_eir),
      function(i) {
        vnapply(
          lagged_transmission_eir[[i]],
          function(l) l$get(timestep - parameters[[i]]$de)
        )
      },
      numeric(n_species)
    )
    if (is.matrix(species_transmission_eir_raw)) {
      species_transmission_eir <- t(species_transmission_eir_raw)
    } else {
      species_transmission_eir <- matrix(species_transmission_eir_raw, ncol = 1L)
    }
    inf <- vnapply(
      seq_along(lagged_infectivity),
      function(i) lagged_infectivity[[i]]$get(timestep - parameters[[i]]$delay_gam)
    )

    if (identical(p_success, 0) || all(p_captured_t == 0)) {
      test_and_treat_coeff <- matrix(1, nrow = n_pops, ncol = n_pops)
    } else {
      rdt_positive <- vnapply(
        seq_len(n_pops),
        function(i) {
          rdt_detectable(variables[[i]], parameters[[i]], timestep)
        }
      )
      test_and_treat_coeff <- (1 - p_captured_t * rdt_positive * p_success)
      diag(test_and_treat_coeff) <- 1
    }

    eir <- vapply(
      seq_len(n_species),
      function(i) {
        mixed_eir <- t(species_eir[,i] * p_mix_import)
        rowSums(mixed_eir * test_and_treat_coeff)
      },
      numeric(n_pops)
    )
    if (!is.matrix(eir)) {
      eir <- matrix(eir, nrow = n_pops, ncol = n_species)
    }
    transmission_eir <- vapply(
      seq_len(n_species),
      function(i) {
        mixed_eir <- t(species_transmission_eir[,i] * p_mix_import)
        rowSums(mixed_eir * test_and_treat_coeff)
      },
      numeric(n_pops)
    )
    if (!is.matrix(transmission_eir)) {
      transmission_eir <- matrix(transmission_eir, nrow = n_pops, ncol = n_species)
    }

    inf <- rowSums(inf * p_mix_export * test_and_treat_coeff)

    list(eir = eir, transmission_eir = transmission_eir, inf = inf)
  }
}

# Estimates RDT prevalence from PCR prevalence,
# We assume PCR prevalence is close to the true proportion of infectious people
# in the population
# values take from Wu et al 2015: https://doi.org/10.1038/nature16039
rdt_detectable <- function(variables, parameters, timestep) {
  infectious_prev <- variables$state$get_size_of(
    c('D', 'A', 'U')) / parameters$human_population
  logit_prev <- log(infectious_prev / (1 - infectious_prev))
  logit_rdt <- parameters$rdt_intercept + parameters$rdt_coeff * logit_prev
  exp(logit_rdt) / (1 + exp(logit_rdt))
}
