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
pnlnm <- "jax.bss" # name of the panel
pnlSim <- simulateMGI(simP.filt[[pnlnm]]$markers,
                      nrow(simP.filt[[pnlnm]]$data),
                      reps = 300)
pnlQuants <- matrix(0, nrow(pnlSim[[1]]), nrow(pnlSim[[1]]))
for (ii in seq_along(pnlSim)) {
    pnlQuants <- pnlQuants + (pnlSim[[ii]] <= simP.corr[[pnlnm]])
}

pdf(paste0(gsub("\\.", "", pnlnm), ".pdf"),
    width = 720, height = 720)
corrImg(pnlQuants/length(pnlSim) - diag(0.5, nrow(pnlQuants)),
        col = c("steelblue", "white", "firebrick"),
        breaks = c(0, 0.025, 0.975, 1))
dev.off()

## what does suppression of zeroes do to the correlations
#mgiCorrs.zero <- lapply(mgiCorrs, zeroEigSuppress)
