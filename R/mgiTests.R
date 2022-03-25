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

## common markers for the bsb panels
bsb.common <- intersect(simP.filt$jax.bsb$markers$symbol,
                        simP.filt$ucla.bsb$markers$symbol)
bsb.markers <- simP.filt$jax.bsb$markers[
                   simP.filt$jax.bsb$markers$symbol %in% bsb.common,]
bsb.Order <- order(as.numeric(bsb.markers$chr))
bsb.TabOrder <- order(as.numeric(unique(bsb.markers$chr)))


## PLOT BSB PANELS ###################################################
pnlnm <- "ucla.bsb" # name of the panel
nsim <- 10000
pal <- colorRampPalette(c("steelblue", "white", "firebrick"))(41)

## plot common markers
png(paste0(gsub("\\.", "", pnlnm), "_common.png"),
    width = 540, height = 540, type = "cairo")
corrImg(simP.corr[[pnlnm]][bsb.common[chrOrder],
                           bsb.common[chrOrder]],
        col = pal,
        breaks = seq(-1, 1, length.out = 42),
        axes = FALSE)
addChromosomes(bsb.markers, chrTabOrder)
dev.off()

## a correlation test plot of chromosomes 2 and 4 of the common markers
bsb.24 <- bsb.markers[bsb.markers$chr %in% c("2","4"),]
## simulate the cross
bsb.24sim <- simulateMGI(bsb.24,
                         ceiling(mean(c(nrow(simP.filt[["jax.bsb"]]$data),
                                        nrow(simP.filt[["ucla.bsb"]]$data)))),
                         reps = nsim)
bsb.24cor <- array(0, dim = c(8, 8, nsim))
for (ii in seq_along(bsb.24sim)) bsb.24cor[,,ii] <- bsb.24sim[[ii]]
rm(bsb.24sim)
## get quantiles
bsb.jaxcorr <- simP.corr[["jax.bsb"]][bsb.24$symbol, bsb.24$symbol]
bsb.uclacorr <- simP.corr[["ucla.bsb"]][bsb.24$symbol, bsb.24$symbol]
## TO FIX: GET QUANTILES
jaxquants <- apply(
uclaquants <- Reduce(function(x, y) { x + (y <= bsb.uclacorr)},
                     bsb.24sim, init = matrix(0, nrow = nrow(bsb.24),
                                              ncol = nrow(bsb.24)))


## Stepping stones:
## correlation matrix -> test plot w/ -1,1 range and coloured cells
## -> current test plot (more diagnostic, almost significance test)


## the correlation test plot for these two data sets combined
png("bsbCorrTest.png", width = 720, height = 720, type = "cairo")
markerNames <- bsb.24$symbol
par(mfrow = c(8,8), mar = c(0.1,0.1,0.1,0.1))
for (ii in 1:8) {
    for (jj in 1:8) {
        tempmean <- mean(c(bsb.jaxcorr[ii,jj],
                           bsb.uclacorr[ii,jj]))
        tempq <- sum(tempmean <= bsb.24cor[ii,jj,])
        if (ii == jj) {
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            text(0.5, 0.5, markerNames[ii], cex = 1.5)
        } else if (ii < jj) {
            par(mar = c(2.1,0.5,0.5,0.5))
            tempdens <- density(bsb.24cor[ii,jj,])
            plot(NA, xlim = range(tempdens$x), ylim = range(tempdens$y),
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            axis(1, at = seq(min(tempdens$x), max(tempdens$x),
                             length.out = 5),
                 labels = round(seq(min(tempdens$x), max(tempdens$x),
                                    length.out = 5),2))
            polygon(tempdens, col = "gray70")
            if (tempq > 0.5*nsim) {
                shadInd <- tempdens$x <= tempmean
            } else {
                shadInd <- tempdens$x >= tempmean
            }
            polygon(c(tempdens$x[shadInd], tempmean),
                    c(tempdens$y[shadInd], 0),
                    col = "gray50")
            abline(v = bsb.jaxcorr[ii,jj], col = "black", lwd = 1)
            abline(v = bsb.uclacorr[ii,jj], col = "black", lwd = 1)
            abline(v = mean(c(bsb.jaxcorr[ii,jj],
                              bsb.uclacorr[ii,jj])),
                   col = "firebrick", lwd = 2)
        } else {
            par(mar = c(0.1,0.1,0.1,0.1))
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            if (tempq > nsim*0.975) {
                tempcol <- "firebrick"
            } else if (tempq < nsim*0.025) {
                tempcol <- "steelblue"
            } else {
                tempcol = NA
            }
            rect(0, 0, 1, 1, col = tempcol)
            text(0.5, 0.5, tempq, cex = 2)
        }
    }
}
dev.off()

## simulating panel crosses
set.seed(30211)
pnlQuants <- matrix(0, nrow(simP.corr[[pnlnm]]),
                    nrow(simP.corr[[pnlnm]]))
for (ii in 1:nsim) { # big memory suck: this is more efficient
    pnlSim <- simulateMGI(simP.filt[[pnlnm]]$markers,
                      nrow(simP.filt[[pnlnm]]$data),
                      reps = 1)
    pnlQuants <- pnlQuants + (pnlSim[[1]] <= simP.corr[[pnlnm]])
    if ((ii %% 100) == 0) {
        cat("\r Done ", ii)
    }
}
diag(pnlQuants) <- nsim/2

png(paste0(gsub("\\.", "", pnlnm), "_sim.png"),
    width = 540, height = 540, type = "cairo")
par(mar = c(0.1,0.9,0.9,0.1))
chrOrder <- order(as.numeric(simP.filt[[pnlnm]]$markers$chr))
chrTabOrder <- order(as.numeric(unique(simP.filt[[pnlnm]]$markers$chr)))
corrImg(pnlSim[[1]][chrOrder, chrOrder],
        col = pal,
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
        col = pal,
        breaks = seq(-1, 1, length.out = 42),
        axes = FALSE)
addChromosomes(simP.filt[[pnlnm]]$markers, chrTabOrder)
dev.off()

## what does suppression of zeroes do to the correlations
#mgiCorrs.zero <- lapply(mgiCorrs, zeroEigSuppress)
