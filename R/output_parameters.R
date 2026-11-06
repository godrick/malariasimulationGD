#' Set age groups for epidemiological outputs
#'
#' Supply age breaks for each output you want to update. Omitting an argument
#' or passing `NULL` keeps that output's existing age groups in `parameters`.
#'
#' @details Give breaks in increasing order, in whole days. Each group starts
#'   at one break and ends the day before the next. For example, `c(0, 10, 100)`
#'   defines two groups: ages 0 through 9 days and 10 through 99 days.
#'
#'   With default [get_parameters()] settings, clinical incidence has no age
#'   groups and prevalence has one group covering ages 730 through 3650 days
#'   inclusive. To specify that same prevalence group using breaks, supply
#'   `prevalence = c(730, 3651)`.
#'
#' @param parameters Model parameters, usually from [get_parameters()].
#' @param age_group Age breaks for population counts.
#' @param incidence Age breaks for incidence (D + Tr + A).
#' @param clinical_incidence Age breaks for symptomatic incidence.
#' @param severe_incidence Age breaks for severe incidence (falciparum only).
#' @param prevalence Age breaks for infections detectable by PCR or light microscopy.
#' @param n_with_hypnozoites Age breaks for counts of people with hypnozoites (vivax only).
#' @param ica Age breaks for average acquired clinical immunity.
#' @param icm Age breaks for average maternal clinical immunity.
#' @param iva Age breaks for average acquired immunity to severe disease (falciparum only).
#' @param ivm Age breaks for average maternal immunity to severe disease (falciparum only).
#' @param id Age breaks for average immunity to detectability (falciparum only).
#' @param ib Age breaks for average blood immunity (falciparum only).
#' @param iaa Age breaks for average acquired antiparasite immunity (vivax only).
#' @param iam Age breaks for average maternal antiparasite immunity (vivax only).
#' @param hypnozoites Age breaks for the average number of hypnozoite batches (vivax only).
#' @return The updated parameter list. For each supplied output,
#'   `*_rendering_min_ages` and `*_rendering_max_ages` hold the inclusive bounds.
#' @examples
#' p <- get_parameters()
#' p <- set_epi_outputs(p, clinical_incidence = c(0, 10, 100))
#' p$clinical_incidence_rendering_min_ages  # 0, 10
#' p$clinical_incidence_rendering_max_ages  # 9, 99
#'
#' # Changing population counts later keeps the clinical incidence groups.
#' p <- set_epi_outputs(p, age_group = c(0, 365))
#' p$clinical_incidence_rendering_max_ages  # 9, 99
#' @export
set_epi_outputs <- function(parameters,
                            age_group = NULL,
                            incidence = NULL,
                            clinical_incidence = NULL,
                            severe_incidence = NULL,
                            prevalence = NULL,
                            n_with_hypnozoites = NULL,
                            ica = NULL,
                            icm = NULL,
                            iva = NULL,
                            ivm = NULL,
                            id = NULL,
                            ib = NULL,
                            iaa = NULL,
                            iam = NULL,
                            hypnozoites = NULL
){
  
  input <- list(
    age_group = age_group,
    incidence = incidence,
    clinical_incidence = clinical_incidence,
    severe_incidence = severe_incidence,
    prevalence = prevalence,
    n_with_hypnozoites = n_with_hypnozoites,
    ica = ica,
    icm = icm,
    iva = iva,
    ivm = ivm,
    id = id,
    ib = ib,
    iaa = iaa,
    iam = iam,
    hypnozoites = hypnozoites
  )
  
  if(parameters$parasite == "falciparum"){
    if(!any(vapply(list(n_with_hypnozoites, iaa, iam, hypnozoites), is.null, logical(1)))){
      stop("Some selected epi outputs do not apply to this parasite.")
    }
    
  } else if (parameters$parasite == "vivax"){
    if(!any(vapply(list(severe_incidence, iva, ivm, id, ib), is.null, logical(1)))){
      stop("Some selected epi outputs do not apply to this parasite.")
    }
  }
  
  input <- input[!sapply(input, is.null)]
  
  for(i in seq_along(input)){
    name <- names(input)[i]
    ages <- input[[i]]
    min_ages <- ages[-length(ages)]
    max_ages <- ages[-1] - 1
    parameters[[paste0(name, "_rendering_min_ages")]] <- min_ages
    parameters[[paste0(name, "_rendering_max_ages")]] <- max_ages
  }

  return(parameters)
}
