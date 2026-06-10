# ------------------------------------------------------------------------------
# homing_drive.R
# ------------------------------------------------------------------------------
# Minimal homing-drive cube helpers for the one-node native backend example.
# The released genotype is HH, and any genotype containing H is treated as a
# drive carrier in the quick output summaries.
# ------------------------------------------------------------------------------

set_genotype_trait_multiplicative <- function(cube, allele_values, trait = "s") {
  allele_values <- unlist(allele_values)
  if (is.null(cube[[trait]])) {
    cube[[trait]] <- stats::setNames(rep(1, length(cube$genotypesID)), cube$genotypesID)
  }

  cube[[trait]] <- stats::setNames(
    vapply(names(cube[[trait]]), function(g) {
      alleles <- strsplit(g, "", fixed = TRUE)[[1]]
      if (!all(alleles %in% names(allele_values))) {
        stop(sprintf(
          "Genotype '%s' contains alleles missing from allele_values: %s",
          g,
          paste(setdiff(alleles, names(allele_values)), collapse = ", ")
        ), call. = FALSE)
      }
      prod(allele_values[alleles])
    }, numeric(1)),
    names(cube[[trait]])
  )

  cube
}

build_one_node_drive_cube <- function(
    baseline_b = 0.55,
    drive_b = 0,
    allele_s = c(W = 1, H = 0.8, B = 0.5, R = 1)
) {
  if (!requireNamespace("MGDrivE", quietly = TRUE)) {
    stop("Package 'MGDrivE' is required for the one-node homing-drive example.",
         call. = FALSE)
  }

  cube <- MGDrivE::cubeHomingDrive(
    cM  = 1.00,
    cF  = 1.00,
    chM = 0.95,
    chF = 0.95,
    crM = 0.05,
    crF = 0.05
  )
  cube$releaseType <- "HH"

  drive_carriers <- grepl("H", cube$genotypesID, fixed = TRUE)
  cube$b <- stats::setNames(
    ifelse(drive_carriers, drive_b, baseline_b),
    cube$genotypesID
  )
  cube$c <- stats::setNames(rep(1, length(cube$genotypesID)), cube$genotypesID)

  if (!is.null(cube$omega)) {
    cube$omega <- stats::setNames(rep(1, length(cube$genotypesID)), cube$genotypesID)
  }

  set_genotype_trait_multiplicative(
    cube = cube,
    allele_values = allele_s,
    trait = "s"
  )
}
