# ------------------------------------------------------------------------------
# homing_drive.R
# ------------------------------------------------------------------------------
# Build the MGDrivE homing-drive cube the example releases.
#
# IMPORTANT — MGDrivE's `cubeHomingDrive()` zero-argument default is NOT a
# working homing drive: the package defaults are `chM = chF = 0` (no homing
# conversion), so cleaved wild-type alleles all become broken (B) and the
# released H allele just dilutes via Mendelian inheritance. To get a drive
# that actually sweeps, you must set the homing rates `chM` and `chF` to
# something > 0.
#
# We use chM = chF = 0.95 (high but not absolute), with cM = cF = 1.0
# (germline cleavage at near-saturation) and a small per-cleavage failure
# rate (crM = crF = 0.05) that produces an R (resistance) allele instead of
# B. These are illustrative high-homing protective-drive settings for a compact
# package example. The resulting cube's genotype set is
#   WW, WH, WR, WB, HH, HR, HB, RR, RB, BB
# where W = wildtype, H = drive, R = resistance, B = broken. The H-carriers
# (i.e. drive carriers) are WH, HH, HR, HB.
# ------------------------------------------------------------------------------

build_seven_node_drive_cube <- function(baseline_b = 0.55, drive_b = 0.0) {
  if (!requireNamespace("MGDrivE", quietly = TRUE)) {
    stop("Package 'MGDrivE' is required.", call. = FALSE)
  }
  cube <- MGDrivE::cubeHomingDrive(
    cM  = 1.00,   # cleavage rate, male germline
    cF  = 1.00,   # cleavage rate, female germline
    chM = 0.95,   # homing rate given cleavage, male germline
    chF = 0.95,   # homing rate given cleavage, female germline
    crM = 0.05,   # resistance allele probability given failed homing, male
    crF = 0.05    # resistance allele probability given failed homing, female
  )
  cube$releaseType <- "HH"

  # Make the drive protective: per-genotype mosquito-to-human transmission
  # probability `b`. Drive-positive genotypes (any genotype with an H allele)
  # carry an anti-Plasmodium effector so their `b` is drive_b (default 0,
  # i.e. fully transmission-blocking). All other genotypes keep the
  # wild-type `b = baseline_b`. Without this step, drive-carriers still
  # transmit at the wild-type rate and the drive sweep would not reduce
  # malaria in the diagnostic panels.
  H_carriers <- grepl("H", cube$genotypesID, fixed = TRUE)
  b_vec <- rep(baseline_b, length(cube$genotypesID))
  b_vec[H_carriers] <- drive_b
  names(b_vec) <- cube$genotypesID
  cube$b <- b_vec
  cube$c <- rep(1, length(cube$genotypesID))   # uniform human-to-mosquito
  names(cube$c) <- cube$genotypesID

  cube
}
