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
bsb.markers <- bsb.markers[bsb.Order,]
 
## the BSB simulations for chromosomes 2, 4, and 18
## set up the simulated marker set
bsb.sub <- bsb.markers[bsb.markers$chr %in% c("2","4","18"),]
ncom <- nrow(bsb.sub) # number of markers
nsim <- 10000 # number of simulated crosses
pal <- colorRampPalette(c("steelblue", "white", "firebrick"))(41)
bsb.subthry <- theoryCor(lapply(split(bsb.sub$cMs,
                                      bsb.sub$chr)[c(2,3,1)],
                                diff))

## simulate the cross with 80 mice
set.seed(30211)
bsb.sim <- simulateMGI(bsb.sub,
                       ceiling(mean(c(nrow(simP.filt[["jax.bsb"]]$data),
                                      nrow(simP.filt[["ucla.bsb"]]$data)))),
                       chrOrd = c(2,3,1),
                       reps = nsim) # simulated crosses
bsb.cor <- array(0, dim = c(ncom, ncom, nsim)) # place in an array
for (ii in seq_along(bsb.sim)) bsb.cor[,,ii] <- bsb.sim[[ii]]
rm(bsb.sim) # remove list (save mem)
## get quantiles
bsb.jaxcorr <- simP.corr[["jax.bsb"]][bsb.sub$symbol, bsb.sub$symbol]
bsb.uclacorr <- simP.corr[["ucla.bsb"]][bsb.sub$symbol, bsb.sub$symbol]
jaxquants <- matrix(rowSums(apply(bsb.cor, 3,
                                  function(mat) mat <= bsb.jaxcorr)),
                    ncol = ncom)
uclaquants <-  matrix(rowSums(apply(bsb.cor, 3,
                                  function(mat) mat <= bsb.uclacorr)),
                      ncol = ncom)

## alternatively: simulate the individual crosses to recombine later
set.seed(10340504)
jaxSim <- simulateMGI(bsb.sub, nrow(simP.filt[["jax.bsb"]]$data),
                      chrOrd = c(2,3,1),
                      reps = nsim)
jaxSimCor <- array(0, dim = c(ncom, ncom, nsim))
for (ii in seq_along(jaxSim)) jaxSimCor[,,ii] <- jaxSim[[ii]]
rm(jaxSim)
uclaSim <- simulateMGI(bsb.sub, nrow(simP.filt[["ucla.bsb"]]$data),
                       chrOrd = c(2,3,1),
                       reps = nsim)
uclaSimCor <- array(0, dim = c(ncom, ncom, nsim))
for (ii in seq_along(uclaSim)) uclaSimCor[,,ii] <- uclaSim[[ii]]
rm(uclaSim) 

## PLOTTING ##########################################################

## density helper to zero for polygons
zeroDens <- function(dens) {
    dens$y <- c(0, dens$y, 0)
    dens$x <- c(dens$x[1], dens$x, dens$x[length(dens$x)])
    dens
}

## start with the distributions for strongly/weakly associated markers
weakEx <- c(1,4)
strngEx <- c(7,8)
## simple hypotheticals: one, two, three, crossovers etc.
corCross <- function(x, ntwo, npop = 80) {
    nonneg <- max(ntwo - x, 0)
    refMrk <- c(rep(1, npop - ntwo), rep(2, ntwo))
    xMrk <- c(rep(1, npop - nonneg), rep(2, nonneg))
    suppressWarnings(cor(refMrk, xMrk))
}
ntwoSeq <- seq(25, 55, by = 5)
xThry <- matrix(0, 5, length(ntwoSeq))
for (ii in 1:5) { for (jj in seq_along(ntwoSeq)) {
                      xThry[ii, jj] <- corCross(ii - 1, ntwoSeq[jj])
                  }
}

## weak density and barplot
#res <- 480
#png("weakDen.png", width = res, height = res, type = "cairo")
side <- 5
ppi <- 150
png("weakDen.png", width = side, height = side, units = "in", res = ppi,
    type = "cairo")
tempdens <- density(bsb.cor[weakEx[1], weakEx[2], ])
plot(NA, xlim = range(tempdens$x), ylim = range(tempdens$y),
     xlab = "Correlation", ylab = "Density")
polygon(tempdens, col = "gray70")
dev.off()
#png("weakBar.png", width = res, height = res, type = "cairo")
png("weakBar.png", width = side, height = side, units = "in", res = ppi,
    type = "cairo")
temptab <- table(bsb.cor[weakEx[1], weakEx[2], ])/nsim
plot(NA, xlim = range(tempdens$x), ylim = range(temptab),
     xlab = "Correlation", ylab = "Proportion")
for (ii in seq_along(temptab)) {
    lines(rep(as.numeric(names(temptab)[ii]), 2),
          c(0, temptab[ii]))
}
dev.off()

## strong density and barplot
#res <- 480
#png("strngDen.png", width = res, height = res, type = "cairo")
png("strngDen.png", width = side, height = side, units = "in", res = ppi,
    type = "cairo")
tempdens <- zeroDens(density(bsb.cor[strngEx[1], strngEx[2], ]))
plot(NA, xlim = range(tempdens$x), ylim = range(tempdens$y),
     xlab = "Correlation", ylab = "Density")
polygon(tempdens, col = "gray70")
dev.off()
#png("strngBar.png", width = res, height = res, type = "cairo")
png("strngBarLndscp.png", width = side, height = 0.8*side, units = "in", res = ppi,
    type = "cairo")
temptab <- table(bsb.cor[strngEx[1], strngEx[2], ])/nsim
plot(NA, xlim = range(as.numeric(names(temptab)[-length(temptab)])),
     ylim = range(temptab[-length(temptab)]),
     xlab = "Correlation", ylab = "Proportion")
pal <- hcl.colors(nrow(xThry), "Set 2")
for (ii in seq_along(pal)) {
    abline(v = xThry[ii,], col = adjustcolor(pal[ii], 0.7))
}
for (ii in seq_along(temptab)) {
    lines(rep(as.numeric(names(temptab)[ii]), 2),
          c(0, temptab[ii]))
}
legend(x = "topleft", legend = 1:(length(pal)-1), fill = pal[-length(pal)],
       title = "Recombinant count")
dev.off()
#png("strngBarClose.png", width = res, height = res, type = "cairo")
png("strngBarClose.png", width = side, height = side, units = "in", res = ppi,
    type = "cairo")
pal <- colorRampPalette(c("steelblue", "gray50", "firebrick"))(length(ntwoSeq))
plot(NA, xlim = c(0.968, 0.976), ylim = c(0,0.07),
     xlab = "Correlation", ylab = "Proportion")
for (ii in seq_along(ntwoSeq)) {
    abline(v = xThry[,ii], col = pal[ii])
}
for (ii in seq_along(temptab)) {
    lines(rep(as.numeric(names(temptab)[ii]), 2),
          c(0, temptab[ii]))
}
legend(x = "topleft", legend = round(2*ntwoSeq/80 +
                                     (80-ntwoSeq)/80, 2),
       fill = pal, title = "Mean genetic score")
dev.off(
)

## plot these cross over curves directly
nCross <- 0:50
ntwoSeq <- seq(1, 80, by = 1)
xThry <- matrix(0, length(nCross), length(ntwoSeq))
for (ii in seq_along(nCross)) {
    for (jj in seq_along(ntwoSeq)) {
                      xThry[ii, jj] <- corCross(nCross[ii], ntwoSeq[jj])
                  }
}
pal <- colorRampPalette(c("black", "firebrick"))(length(nCross))
#png("crossCurves.png", width = 540, height = 540, type = "cairo")
png("crossCurves.png", width = side, height = side, units = "in", res = ppi,
    type = "cairo")
plot(NA, xlim = range(ntwoSeq),
     ylim = range(c(xThry, 1.03), na.rm = TRUE),
     xlab = "Mean score", ylab = "Correlation", xaxt = "n")
axis(1, at = seq(0, 80, by = 10),
     labels = round(1+seq(0, 80, by = 10)/80,2))
text(x = 0, y = 1.02, "Recombinant offspring", pos = 4)
text(x = nCross, y = apply(xThry, 1, min, na.rm = TRUE),
     labels = nCross)
for (ii in seq_along(nCross)) lines(ntwoSeq, xThry[ii,], col = pal[ii])
dev.off()

## add some observed markers
jaxCommon <- mgiSubset(simP.filt$jax.bsb, bsb.sub$symbol)
jaxPts <- mgiCorrMeans(jaxCommon)
xThry <- matrix(0, length(nCross), length(ntwoSeq))
for (ii in seq_along(nCross)) {
    for (jj in seq_along(ntwoSeq)) {
        xThry[ii, jj] <- corCross(nCross[ii], ntwoSeq[jj],
                                  npop = nrow(jaxCommon$data))
                  }
}
#png("jaxcrossCurves.png", width = 540, height = 540, type = "cairo")
png("jaxcrossCurves.png", width = side, height = side, units = "in", res = ppi,
    type = "cairo")
plot(NA, xlim = range(ntwoSeq),
     ylim = range(c(xThry, 1.03), na.rm = TRUE),
     xlab = "Mean score", ylab = "Correlation", xaxt = "n")
axis(1, at = seq(0, 80, by = 10),
     labels = round(1+seq(0, 80, by = 10)/80,2))
text(x = 0, y = 1.02, "Recombinant offspring", pos = 4)
inds <- seq(1, length(nCross), by = 5)
text(x = nCross[inds], y = apply(xThry, 1, min, na.rm = TRUE)[inds],
    labels = nCross[inds])
for (ii in seq_along(nCross)) lines(ntwoSeq, xThry[ii,], col = pal[ii])
for (ii in 1:3) points(80*(jaxPts$means[[ii]]-1), jaxPts$corr[[ii]])
dev.off()



## corr dist (change ncom = 2 for 2x2, ncom = 8 chromosomes 2 & 4)
ncom <- 8 # nrow(bsb.sub) # 2 # CHANGE TO CHANGE DISPLAY
#res <- 720 # better for big # 540 # better for 2x2
pal <- colorRampPalette(c("steelblue", "white", "firebrick"))(41)
#png("bsbCorrDist.png", width = res, height = res, type = "cairo")
## width/height = 10 for big, 2.5 for 2x2
png("bsbCorrDist.png", width = 8, height = 8, units = "in", res = ppi,
    type = "cairo")
markerNames <- bsb.sub$symbol
cuts <- seq(-1, 1, length.out = 42)
par(mfrow = c(ncom,ncom), mar = c(0.1,0.1,0.1,0.1))
for (ii in 1:ncom) { # 7:8 for 2x2
    for (jj in 1:ncom) { # 7:8 for 2x2
        tempmean <- mean(bsb.cor[ii,jj,])
        if (ii == jj) {
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            text(0.5, 0.5, markerNames[ii], cex = 1.5)
        } else if (ii < jj) {
            par(mar = c(2.1,0.5,0.5,0.5))
            tempdens <- zeroDens(density(bsb.cor[ii,jj,]))
            plot(NA, xlim = c(-1,1), ylim = c(0,7),
                 yaxt = "n", xlab = "", ylab = "")
            polygon(tempdens, col = "gray70")
            abline(v = 0)
            abline(v = bsb.subthry[ii,jj], col = "firebrick",
                   lty = 1)
            abline(v = tempmean, lty = 2)
        } else {
            par(mar = c(0.1,0.1,0.1,0.1))
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            rect(0, 0, 1, 1,
                 col = pal[sum(cuts < tempmean)])
            text(0.5, 0.5, round(tempmean, 2), cex = 2)
        }
    }
}
dev.off()

## the correlation test plot for these two data sets combined
#png("bsbCorrTest.png", width = res, height = res, type = "cairo")
png("bsbCorrTest.png", width = 10, height = 10, units = "in", res = ppi,
    type = "cairo")
par(mfrow = c(ncom,ncom), mar = c(0.1,0.1,0.1,0.1))
for (ii in 1:ncom) { 
    for (jj in 1:ncom) {
        tempmean <- mean(c(bsb.jaxcorr[ii,jj],
                           bsb.uclacorr[ii,jj]))
        tempq <- sum(tempmean >= bsb.cor[ii,jj,])
        if (ii == jj) {
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            text(0.5, 0.5, markerNames[ii], cex = 1.5)
        } else if (ii < jj) {
            par(mar = c(2.1,0.5,0.5,0.5))
            tempdens <- zeroDens(density(bsb.cor[ii,jj,]))
            plot(NA, xlim = range(tempdens$x), ylim = range(tempdens$y),
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            axis(1, at = seq(min(tempdens$x), max(tempdens$x),
                             length.out = 5),
                 labels = round(seq(min(tempdens$x), max(tempdens$x),
                                    length.out = 5),2))
            polygon(tempdens, col = "gray70")
            #if (tempq <= 0.5*nsim) {
                shadInd <- tempdens$x <= tempmean
            #} else {
            #    shadInd <- tempdens$x >= tempmean
            #}
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

## another one: this one giving the mean distribution across sims
#png("bsbCorrTest.png", width = 720, height = 720, type = "cairo")
png("bsbCorrTest.png", width = 8, height = 8, units = "in", res = ppi,
    type = "cairo")
markerNames <- bsb.sub$symbol
cuts <- seq(-1, 1, length.out = 42)
par(mfrow = c(ncom,ncom), mar = c(0.1,0.1,0.1,0.1))
for (ii in 1:ncom) {
    for (jj in 1:ncom) {
        tempjax <- mean(jaxSimCor[ii,jj,])
        tempucla <- mean(uclaSimCor[ii,jj,])
        tempmean <- mean(c(bsb.jaxcorr[ii,jj], bsb.uclacorr[ii,jj]))
        tempcomb <- 0.5*jaxSimCor[ii,jj,] + 0.5*uclaSimCor[ii,jj,]
        if (ii == jj) {
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            text(0.5, 0.5, markerNames[ii], cex = 1.5)
        } else if (ii < jj) {
            par(mar = c(2.1,0.5,0.5,0.5))
            tempdens <- density(tempcomb)
            plot(NA, xlim = range(tempdens$x), ylim = range(tempdens$y),
                 yaxt = "n", xlab = "", ylab = "")
            polygon(tempdens, col = "gray70")
            shadInd <- tempdens$x <= tempmean
            polygon(c(tempdens$x[shadInd], tempmean),
                    c(tempdens$y[shadInd], 0),
                    col = "gray50")
            abline(v = bsb.jaxcorr[ii,jj], lty = 2)
            abline(v = bsb.uclacorr[ii,jj], lty = 4)
            abline(v = tempmean, col = "firebrick", lwd = 2)
        } else { 
            tempq <- sum(tempcomb <= mean(c(bsb.jaxcorr[ii,jj],
                                           bsb.uclacorr[ii,jj])))
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

## PLOT BSB PANELS ###################################################

pnlnm <- "ucla.bsb" # name of the panel
## plot common markers
png(paste0(gsub("\\.", "", pnlnm), "_common.png"),
    width = 540, height = 540, type = "cairo")
corrImg(simP.corr[[pnlnm]][bsb.common[bsb.Order],
                           bsb.common[bsb.Order]],
        col = pal,
        breaks = seq(-1, 1, length.out = 42),
        axes = FALSE)
addChromosomes(bsb.markers, bsb.TabOrder)
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

## the correlation distribution plot of these results
png("pairedSimBig.png", width = 1440, height = 1440, type = "cairo")
markerNames <- bsb.sub$symbol
cuts <- seq(-1, 1, length.out = 42)
par(mfrow = c(ncom,ncom), mar = c(0.1,0.1,0.1,0.1))
for (ii in 1:ncom) {
    for (jj in 1:ncom) {
        tempjax <- mean(jaxSimCor[ii,jj,])
        tempucla <- mean(uclaSimCor[ii,jj,])
        if (ii == jj) {
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            text(0.5, 0.5, markerNames[ii], cex = 1.5)
        } else if (ii < jj) {
            par(mar = c(2.1,0.5,0.5,0.5))
            tempjaxdens <- density(jaxSimCor[ii,jj,])
            tempucladens <- density(uclaSimCor[ii,jj,])
            plot(NA, xlim = range(tempjaxdens$x, tempucladens$x),
                 ylim = range(tempjaxdens$y, tempucladens$y),
                                        #xlim = c(-1,1), ylim = c(0,7),
                 yaxt = "n", xlab = "", ylab = "")
            #abline(v = bsb.subthry[ii,jj], col = "firebrick",
            #       lty = 1)
            polygon(tempjaxdens, lty = 2,
                    col = adjustcolor("#d95f02", alpha.f = 0.5))
            polygon(tempucladens, lty = 3,
                    col = adjustcolor("#1b9e77", alpha.f = 0.5))
            abline(v = bsb.jaxcorr[ii,jj], col = "#d95f02", lwd = 1)
            abline(v = bsb.uclacorr[ii,jj], col = "#1b9e77", lwd = 1)
            abline(v = mean(c(bsb.jaxcorr[ii,jj],
                              bsb.uclacorr[ii,jj])),
                   col = "firebrick", lwd = 2)
            #abline(v = 0)
            #abline(v = tempjax, lty = 2)
            #abline(v = tempucla, lty = 3)
        } else {
            tempqjax <- sum(jaxSimCor[ii,jj,] < bsb.jaxcorr[ii,jj])
            tempqucla <- sum(uclaSimCor[ii,jj,] < bsb.uclacorr[ii,jj])
            par(mar = c(0.1,0.1,0.1,0.1))
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            if (tempqjax > nsim*0.975) {
                tempcoljax <- "firebrick"
            } else if (tempqjax < nsim*0.025) {
                tempcoljax <- "steelblue"
            } else {
                tempcoljax = NA
            }
            if (tempqucla > nsim*0.975) {
                tempcolucla <- "firebrick"
            } else if (tempqucla < nsim*0.025) {
                tempcolucla <- "steelblue"
            } else {
                tempcolucla = NA
            }
            rect(0, 0, 1, 0.5, col = tempcoljax)
            text(0.5, 0.25, tempqjax, cex = 2, col = "#d95f02")
            rect(0, 0.5, 1, 1, col = tempcolucla)
            text(0.5, 0.75, tempqucla, cex = 2, col = "#1b9e77")
        }
    }
}
dev.off()
