## load functions, set image output directory
source("simulationFunctions.R")
source("mgiFunctions.R")
imgDir <- "../img/"
dataDir <- "../data/"

## READ/PROCESS THE DATA #############################################

##' another source of experimental data is the MGI website, which has
##' a number of "panels" of mice generated in experimental work and
##' then measured, typically the crosses used are backcrosses
##' the functionality to extract all of this and the reference website
##' are in the simulationFunctions file, we start by loading the MGI
##' provided database of markers (which is rather large)
mgiMarkers <- readMGIlists() # descriptions of all markers

##' next we can consider the panels, first read the data
mgiPanels <- lapply(paste0(mgiUrl, mgiPNames), readMGIrpt)
names(mgiPanels) <- names(mgiPNames)
##' match the marker names back to the marker reference
mgiPanels.mrkr <- lapply(mgiPanels,
                         function(panel) {
                             match(names(panel$data),
                                   mgiMarkers$`Marker Symbol`)
                         })
##' get cM positions
mgiPanels.cMs <- lapply(mgiPanels.mrkr, processcMs,
                        tab = mgiMarkers$`cM Position`)
##' Infs indicate markers localized to a chromosome but without a
##' known cM position, while NAs indicate markers missing from the
##' reference file or marked as unknown there

##' filter the panels by those markers with known positions
mgiFiltered <- mapply(filterPanelBycM, mgiPanels, mgiPanels.cMs,
                      SIMPLIFY = FALSE)

##' identify the type of cross of each using the summaries
mgiPanel.cross <- c("backcross", "backcross", "backcross",
                    "backcross", "backcross", "backcross",
                    "backcross", "unclear", "backcross",
                    "backcross", "intercross", "unclear",
                    "backcross", "backcross")
names(mgiPanel.cross) <- names(mgiFiltered)

## remove unnecessary data
rm(list = c("mgiMarkers", "mgiPanels.cMs", "mgiPanels.mrkr",
            "mgiPanels"))


## CORRELATION FOR ALL DATA ##########################################

## get observed correlation matrices
mgiCorrs <- lapply(mgiFiltered, mgiCorrelation)
## get theoretical correlations
mgiCorrs.th <- mapply(mgiTheory, mgiFiltered, mgiPanel.cross)


## FILTERING AND SIMULATING A SUBSET #################################

## simulate those that have large enough sample (removing bad markers)
## preparation
simPanels <- mgiFiltered[c("jax.bsb", "jax.bss", "mit", "ucla.bsb")]
simPanels <- lapply(simPanels, mgiDropZeroPanelist)
simP.filt <- lapply(simPanels, mgiDropBadMarker)
simP.corr <- lapply(simP.filt, mgiCorrelation, use = "all.obs")
simP.th <- mapply(mgiTheory, simP.filt,
                  mgiPanel.cross[names(simPanels)])

## simulating these is a huge memory suck, so use the loop here
set.seed(30211)
pnlnm <- "ucla.bsb" # name of the panel
#pnlSim <- simulateMGI(simP.filt[[pnlnm]]$markers,
#                      nrow(simP.filt[[pnlnm]]$data),
#                      reps = 1000)

nsim <- 5000
pnlQuants <- matrix(0, nrow(simP.corr[[pnlnm]]),
                    nrow(simP.corr[[pnlnm]]))
for (ii in 1:nsim) {
    pnlSim <- simulateMGI(simP.filt[[pnlnm]]$markers,
                      nrow(simP.filt[[pnlnm]]$data),
                      reps = 1)
    pnlQuants <- pnlQuants + (pnlSim[[1]] <= simP.corr[[pnlnm]])
    if ((ii %% 100) == 0) {
        cat("\r Done ", ii)
    }
}

png(paste0(gsub("\\.", "", pnlnm), "_sim.png"),
    width = 540, height = 540, type = "cairo")
par(mar = c(0.1,0.9,0.9,0.1))
chrOrder <- order(as.numeric(simP.filt[[pnlnm]]$markers$chr))
chrTabOrder <- order(as.numeric(unique(simP.filt[[pnlnm]]$markers$chr)))
corrImg(pnlSim[[1]][chrOrder, chrOrder],
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(41),
        breaks = seq(-1, 1, length.out = 42),
        axes = FALSE)
addChromosomes(simP.filt[[pnlnm]]$markers, chrTabOrder)
dev.off()

png(paste0(gsub("\\.", "", pnlnm), "_quant.png"),
    width = 540, height = 540, type = "cairo")
par(mar = c(0.1,0.9,0.9,0.1))
chrOrder <- order(as.numeric(simP.filt[[pnlnm]]$markers$chr))
chrTabOrder <- order(as.numeric(unique(simP.filt[[pnlnm]]$markers$chr)))
corrImg(pnlQuants[chrOrder, chrOrder],
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(3),
        breaks = c(0, 0.025, 0.975, 1)*nsim,
        axes = FALSE)
addChromosomes(simP.filt[[pnlnm]]$markers, chrTabOrder)
dev.off()

png(paste0(gsub("\\.", "", pnlnm), ".png"),
    width = 540, height = 540, type = "cairo")
par(mar = c(0.1,0.9,0.9,0.1))
chrOrder <- order(as.numeric(simP.filt[[pnlnm]]$markers$chr))
chrTabOrder <- order(as.numeric(unique(simP.filt[[pnlnm]]$markers$chr)))
corrImg(simP.corr[[pnlnm]][chrOrder, chrOrder],
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(41),
        breaks = seq(-1, 1, length.out = 42),
        axes = FALSE)
addChromosomes(simP.filt[[pnlnm]]$markers, chrTabOrder)
dev.off()

## what does suppression of zeroes do to the correlations
#mgiCorrs.zero <- lapply(mgiCorrs, zeroEigSuppress)
