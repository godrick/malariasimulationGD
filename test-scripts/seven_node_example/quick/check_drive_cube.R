####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
outFolder <- "mgdrive"

####################
# Simulation Parameters
####################
# days to run the simulation, 2 years
tMax <- 365*10

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- as.matrix(1)

# parameters of the population equilibrium
adultPopEquilibrium <- 1000
sitesNumber <- nrow(moveMat)
patchPops <- rep(adultPopEquilibrium,sitesNumber)

####################
# Basic Inheritance pattern
####################
# Reciprocal translocation cube with standard (default) parameters
# cube <- cubeReciprocalTranslocations()

cube <- cubeHomingDrive(
  cM  = 1.00,   # cleavage rate, male germline
  cF  = 1.00,   # cleavage rate, female germline
  chM = 0.95,   # homing rate given cleavage, male germline
  chF = 0.95,   # homing rate given cleavage, female germline
  crM = 0.05,   # resistance allele probability given failed homing, male
  crF = 0.05    # resistance allele probability given failed homing, female
)


set_genotype_trait_multiplicative <- function(cube, allele_values, trait = "s") {
  allele_values <- unlist(allele_values)
  
  cube[[trait]] <- sapply(names(cube[[trait]]), function(g) {
    alleles <- strsplit(g, "")[[1]]
    prod(allele_values[alleles])
  })
  
  cube
}

allele_s <- c(W = 1, H = 0.8, B = 0.5, R = 1)

cube <- set_genotype_trait_multiplicative(
  cube = cube,
  allele_values = allele_s,
  trait = "s"
)

cube$s

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur once a week, for 6 weeks.
#  There are 100 mosquitoes released every time.
releasesParameters <- list(releasesStart=100,
                           releasesNumber=1,
                           releasesInterval=1,
                           releaseProportion=10)

# generate release vector
releasesVector <- generateReleaseVector(driveCube=cube,
                                        releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- releasesVector
patchReleases[[1]]$femaleReleases <- releasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5),
                                      numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                              beta=bioParameters$betaK, muAd=bioParameters$muAd,
                              popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                              tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                              AdPopEQ=patchPops, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=outFolder,
                          verbose = FALSE)
# run simulation
MGDrivESim$oneRun(verbose = FALSE)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                 remFile = TRUE, verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)