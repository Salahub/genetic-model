## load package with data, genome functions
library(simpleGenome)

## FUNCTIONS TO ANALYZE, SIMULATE MGI PANELS #########################
## drop bad panelists
mgiDropZeroPanelist <- function(panel) {
    filterPopulation(panel,
                     rule = function(enc) !all(is.na(enc[,1])) &
                                          !all(is.na(enc[,2])))
}

## filter bad markers
mgiDropBadMarker <- function(panel, prop = 0) {
    bad <- apply(do.call(cbind, panel$encodings), 1,
                 function(row) {
                     sum(is.na(row))/(2*length(panel$encodings)) > prop
                 })
    subsetPopulation(panel, inds = !bad)
}

## simulate panel correlation based on its setting and cM distances
simulateMGICor <- function(panel, reps = 1000, meioseArgs = list(),
                           setting = c("backcross", "intercross"),
                           asArray = FALSE) {
    setting <- match.arg(setting) # identify case
    M <- with(panel, makeGenome(location, alleles, chromosome,
                                markerFun = markerPureDom))
    F <- with(panel, makeGenome(location, alleles, chromosome,
                                markerFun = markerPureRec))
    F1 <- do.call(sex, args = c(genome1 = M, genome2 = F, meioseArgs))
    allcors <- vector("list", reps)
    if (setting == "intercross") {
        for (ii in 1:reps) {
            pop <- asPopulation(replicate(npop, sex(F1, F1),
                                          simplify = FALSE))
            allcors[[ii]] <- popCorrelation(pop)
            if (ii %% 100 == 0) cat("\r -- Simulated population ", ii)
        }
    } else if (setting == "backcross") {
        for (ii in 1:reps) {
            pop <- asPopulation(replicate(npop, sex(F1, M),
                                          simplify = FALSE))
            allcors[[ii]] <- popCorrelation(pop)
            if (ii %% 100 == 0) cat("\r -- Simulated population ", ii)
        }
    } else {
        stop("Setting unknown")
    }
    if (asArray) { # place in an appropriate array
        tempcor <- array(0, dim = c(ncom, ncom, nsim))
        for (ii in seq_along(allcors)) tempcor[,,ii] <- allcors[[ii]]
        allcors <- tempcor # for output
    }
    allcors
}

## suppress zeros and reconstruct a correlation matrix
zeroEigSuppress <- function(sigma) {
    decom <- eigen(sigma)
    eigs <- decom$values
    eigs[eigs < 0] <- 0 # zero negatives
    recon <- decom$vectors %*% diag(eigs) %*% t(decom$vectors)
    recon
}

## compute mgi correlation and report population mean for marker pairs
mgiCorrMeans <- function(panel, use = "pairwise.complete.obs",
                         method = "pearson") {
    dropLowTri <- function(mat) { # keep upper entries
        mat[!upper.tri(mat)] <- NA
        mat
    }
    srtd <- panel$data[, order(panel$markers$chr)] # group chromosomes
    dotDrop <- apply(srtd, 2, # turn dots to NA
                     function(mrk) {
                         temp <- mrk
                         temp[temp == "."] <- NA
                         temp
                     })
    numer <- apply(dotDrop, 2, function(mrk) unclass(factor(mrk)))
    splt <- lapply(split(panel$markers$symbol, panel$markers$chr),
                   function(mrks) numer[, mrks]) # intra chromosome
    cors <- lapply(splt, cor, use = use, method = method)
    cors <- lapply(cors, dropLowTri)
    means <- lapply(splt, colMeans) # means by marker
    pwmns <- lapply(means, function(mn) outer(mn, mn, pmax)) # ref mean
    pwmns <- lapply(pwmns, dropLowTri)
    list(corr = cors, means = pwmns)
}


## FILTERING AND DISPLAYING A SUBSET #################################

## load in the pre-processed panels, focusing on a few with large
## samples when controlling for complete data alone
data(jax_bsb); data(jax_bss); data(ucla_bsb)
panels <- list("jax_bsb" = jax_bsb, "jax_bss" = jax_bss,
               "ucla_bsb" = ucla_bsb)
## remove bad markers, panelists
panelsFilt <- lapply(panels, mgiDropZeroPanelist)
panelsComplete <- lapply(panelsFilt, mgiDropBadMarker)
## compute correlation
panelsCorr <- lapply(panelsComplete, popCorrelation)
## visualize a chosen panel
corrImg(panelsCorr[["ucla_bsb"]], xaxt = "n", yaxt = "n")
addChromosomeLines(panelsComplete[["ucla_bsb"]], lncol = "black")

## get theoretical correlations
panelsTheory <- mapply(theoryCorrelation, genome = panelsComplete,
                       setting = "backcross")
## visualize again
corrImg(panelsTheory[["ucla_bsb"]], xaxt = "n", yaxt = "n")
addChromosomeLines(panelsComplete[["ucla_bsb"]], lncol = "black")


## SIMULATIONS #######################################################

## subset BSB panels by common markers
bsbCommon <- intersect(panelsComplete$jax_bsb$marker,
                       panelsComplete$ucla_bsb$marker)
bsbPanels <- with(panelsComplete,
                  list(jax = subsetPopulation(jax_bsb,
                                inds = jax_bsb$marker %in% bsbCommon),
                       ucla = subsetPopulation(ucla_bsb,
                                inds = ucla_bsb$marker %in% bsbCommon)))

## perform simulations for chromosomes 2, 4, and 18
## set up the simulated marker set
bsbSubset <- lapply(bsbPanels, subsetPopulation,
                    inds = bsbPanels$jax$chromosome %in% c("2","4","18"))
nsim <- 10000 # number of simulated crosses
pal <- colorRampPalette(c("steelblue", "white", "firebrick"))(41)
## theoretical correlation
bsb.subthry <- theoryCorrelation(bsbSubset$jax, setting = "backcross")

## simulate the cross with 80 mice
set.seed(30211)
bsb.cor <- simulateMGI(bsb.sub,
                       ceiling(mean(c(nrow(simP.filt[["jax.bsb"]]$data),
                                      nrow(simP.filt[["ucla.bsb"]]$data)))),
                       chrOrd = c(2,3,1),
                       reps = nsim, asArray = TRUE) # simulated crosses
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
jaxSimCor <- simulateMGI(bsb.sub, nrow(simP.filt[["jax.bsb"]]$data),
                         chrOrd = c(2,3,1),
                         reps = nsim, asArray = TRUE)
uclaSimCor <- simulateMGI(bsb.sub, nrow(simP.filt[["ucla.bsb"]]$data),
                          chrOrd = c(2,3,1),
                          reps = nsim, asArray = TRUE)

## PLOTTING ##########################################################

## density helper to zero polygon bases
zeroDens <- function(dens) {
    dens$y <- c(0, dens$y, 0)
    dens$x <- c(dens$x[1], dens$x, dens$x[length(dens$x)])
    dens
}

## simple recombinant hypotheticals
corCross <- function(x, ntwo, npop = 80) {
    nonneg <- max(ntwo - x, 0) # cannot have negative values
    refMrk <- c(rep(1, npop - ntwo), rep(2, ntwo)) # reference
    xMrk <- c(rep(1, npop - nonneg), rep(2, nonneg)) # swapped
    suppressWarnings(cor(refMrk, xMrk)) # avoid NaN warnings
}

## start with the distributions for strongly/weakly associated markers
weakEx <- c(1,4)
strngEx <- c(7,8)
ntwoSeq <- seq(25, 55, by = 5)
xThry <- matrix(0, 5, length(ntwoSeq))
for (ii in 1:5) { for (jj in seq_along(ntwoSeq)) {
                      xThry[ii, jj] <- corCross(ii - 1, ntwoSeq[jj])
                  }
}

## weakly related density
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

## barplot for weakly related markers
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

## strongly related density
#res <- 480
#png("strngDen.png", width = res, height = res, type = "cairo")
png("strngDen.png", width = side, height = side, units = "in", res = ppi,
    type = "cairo")
tempdens <- zeroDens(density(bsb.cor[strngEx[1], strngEx[2], ]))
plot(NA, xlim = range(tempdens$x), ylim = range(tempdens$y),
     xlab = "Correlation", ylab = "Density")
polygon(tempdens, col = "gray70")
dev.off()

## strongly related barplot
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
    lines(rep(as.numeric(names(temptab)[ii]), 2), c(0, temptab[ii]))
}
legend(x = "topleft", legend = 1:(length(pal)-1), fill = pal[-1],
       title = "Recombinant count")
dev.off()

## focused view on one recombinant individual
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
    lines(rep(as.numeric(names(temptab)[ii]), 2), c(0, temptab[ii]))
}
legend(x = "topleft", legend = round(2*ntwoSeq/80 + (80-ntwoSeq)/80, 2),
       fill = pal, title = "Mean genetic score")
dev.off(
)

## plotting the recombinant curves directly
nCross <- 0:50 # denser numbers of recombinants
ntwoSeq <- seq(1, 80, by = 1) # denser means
xThry <- matrix(0, length(nCross), length(ntwoSeq))
for (ii in seq_along(nCross)) {
    for (jj in seq_along(ntwoSeq)) {
                      xThry[ii, jj] <- corCross(nCross[ii], ntwoSeq[jj])
                  }
}
pal <- colorRampPalette(c("black", "firebrick"))(length(nCross))

## plot the above
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

## all together
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
