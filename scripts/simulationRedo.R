## load functions, set image output directory
library(toyGenomeGenR)

## SIMULATIONS #######################################################
## start incredibly simple: a single chromosome, markers equidistant
nmark <- 20
locs <- list(cumsum(rep(15, nmark)))
alls <- rep(list(c("A","a")), nmark)
chr <- factor(rep("1", nmark))
mother <- makeGenome(locs, alls, chr, markerPureDom(nmark))
father <- makeGenome(locs, alls, chr, markerPureRec(nmark))
F1 <- sex(mother, father) # this always gives the same genome
npop <- 1e5
interCross1 <- asPopulation(replicate(npop, sex(F1, F1),
                                      simplify = FALSE))
backCross1 <- asPopulation(replicate(npop, sex(F1, mother),
                                     simplify = FALSE))
interCor1 <- popCorrelation(interCross1)
backCor1 <- popCorrelation(backCross1)
corPal <- colorRampPalette(c("steelblue", "white", "firebrick"))(41)
corBrk <- seq(-1, 1, length.out = 42)
png( "inter1.png") # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(interCor1, xaxt = "n", yaxt = "n", main = "",
        col = corPal, breaks = corBrk)
dev.off() # for output
png("back1.png") # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(backCor1, xaxt = "n", yaxt = "n", main = "",
        col = corPal, breaks = corBrk)
dev.off() # for output

## now consider the case where markers are not equidistant
locs <- list(cumsum(c(rep(2,5), rep(5,5), rep(10,5),
                      rep(20,5)))) # and varied dists
mother <- makeGenome(locs, alls, chr,
                     markerPureDom(sum(sapply(locs, length))))
father <- makeGenome(locs, alls, chr,
                     markerPureRec(sum(sapply(locs, length))))
F1 <- sex(mother, father)
interCross2 <- asPopulation(replicate(npop, sex(F1, F1),
                                      simplify = FALSE))
backCross2 <- asPopulation(replicate(npop, sex(F1, mother),
                                     simplify = FALSE))
interCor2 <- popCorrelation(interCross2)
backCor2 <- popCorrelation(backCross2)
png("inter2.png") # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(interCor2, xaxt = "n", yaxt = "n", main = "",
        col = corPal, breaks = corBrk)
dev.off() # for output
png("back2.png") # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(backCor2, xaxt = "n", yaxt = "n", main = "",
        col = corPal, breaks = corBrk)
dev.off() # for output

## now multiple chromosomes, but keep distances constant
loc <- list(cumsum(rep(5, 8)), cumsum(rep(15, 12)))
chr <- factor(c(rep("1", 8), rep("2", 12)))
alls <- rep(list(c("A", "a")), 20)
mother <- makeGenome(loc, alls, chr, markerPureDom(length(chr)))
father <- makeGenome(loc, alls, chr, markerPureRec(length(chr)))
F1 <- sex(mother, father)
interCross3 <- asPopulation(replicate(npop, sex(F1, F1),
                                      simplify = FALSE))
backCross3 <- asPopulation(replicate(npop, sex(F1, mother),
                                     simplify = FALSE))
interCor3 <- popCorrelation(interCross3)
backCor3 <- popCorrelation(backCross3)
png("inter3.png") # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(interCor3, xaxt = "n", yaxt = "n", main = "",
        col = corPal, breaks = corBrk)
dev.off() # for output
png("back3.png") # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(backCor3, xaxt = "n", yaxt = "n", main = "",
        col = corPal, breaks = corBrk)
dev.off() # for output

## SIMULATIONS EMULATING OTHER WORK ##################################
## Cheverud 2001: factorial design with one chromosome and different
## distances, this is no different than the first case above, but
## insert his particulars
set.seed(201311) # reproducibility
npop <- 500
loc <- list(cumsum(rep(6.25, 17)))
alls <- rep(list(c("A", "a")), 17)
chr <- factor(rep("1", 17))
mother <- makeGenome(loc, alls, chr, markerPureDom(length(chr)))
father <- makeGenome(loc, alls, chr, markerPureRec(length(chr)))
F1 <- sex(mother, father)
cheverud <- asPopulation(replicate(npop, sex(F1, F1),
                                   simplify = FALSE))
cheverudCor <- popCorrelation(cheverud)
png("chevSim.png") # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(cheverudCor, axes = FALSE, main = "",
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(41),
        breaks = seq(-1, 1, length.out = 42))
dev.off()
## theory
cheverudCorTh <- theoryCorrelation(F1)
png("chevSimTheory.png") # for output
par(mar = c(0.5,0.5,0.5,0.5)) # if no title is desired
corrImg(cheverudCorTh, axes = FALSE, main = "",
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(41),
        breaks = seq(-1, 1, length.out = 42))
dev.off()

## Landner and Botstein 1989: 12 chromosomes, each 100 cM long with
## markers placed along each chromosome at 20 cM intervals, backcross
set.seed(23031) # reproducibility
npop <- 250
chr <- factor(rep(1:20, each = 6))
alls <- rep(list(c("A","a")), length(chr))
loc <- lapply(table(chr), function(x) 20*(0:(x-1)))
mother <- makeGenome(loc, alls, chr, markerPureDom(length(chr)))
father <- makeGenome(loc, alls, chr, markerPureRec(length(chr)))
F1 <- sex(mother, father)
landbot <- asPopulation(replicate(npop, sex(F1, mother),
                                  simplify = FALSE))
landbotCor <- popCorrelation(landbot)
png("LBSim.png") # for output
par(mar = c(0.1,0.9,0.9,0.1)) # if no title is desired
corrImg(landbotCor, axes = FALSE, main = "",
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(41),
        breaks = seq(-1, 1, length.out = 42))
addChromosomeLines(F1)
dev.off()
## theory
landbotCorTh <- theoryCorrelation(F1)
png("LBSimTheory.png") # for output
par(mar = c(0.1,0.9,0.9,0.1)) # if no title is desired
corrImg(landbotCorTh, axes = FALSE, yaxt = "n", main = "",
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(41),
        breaks = seq(-1, 1, length.out = 42))
addChromosomeLines(F1)
dev.off()

## Li and Ji 2005: ten independent regions which have equidistant
## elements, "LD r^2 = 0.8" (r^2 = (p_AB - p_A*p_B)/
## (p_A(1-p_A)*p_B(1-p_B))
## interesting note: suggest a different pairwise measure (chi^2)
## - Hardy-Weinberg equilibrium on all alleles (p_A = 1 - p_a)
## - minor alleles all have same frequency (i.e. p_a = p_b)
## making some assumptions:
## - AB was the initial form of these SNPs (occurred together)
## - the population we are modelling is one generation later
##   (for consistency with other one-generation sims)
## we get p_r = 1 - sqrt((r^2))
npop <- 400
pr <- 1 - sqrt(0.8)
chr <- factor(rep(1:10, each = 5))
alls <- rep(list(c("A", "a")), length(chr))
loc <- lapply(table(chr), # invert prob to distance
              function(x) cumsum(rep(invHaldane(pr), times = x)))
mother <- makeGenome(loc, alls, chr, markerPureDom(length(chr)))
father <- makeGenome(loc, alls, chr, markerPureRec(length(chr)))
F1 <- sex(mother, father)
liji <- asPopulation(replicate(npop, sex(F1, F1),
                               simplify = FALSE))
lijiCor <- popCorrelation(liji)
png(paste(imgDir, "liji.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(lijiCor)
dev.off()

## VERIFYING OUR DISTANCE CALCULATION ################################
##' cheverud cites earlier work by him and his colleagues in which
##' they perform exactly the F2 intercross of mice strains that is
##' the basis of these models, this reproduces the results
chevPos <- c(0, 9.4, 29.5, 41.4, 52.1, 73.3, 115.2, 131.7)
xpos <- seq(0, 1, length.out = length(chevPos))
chr <- factor(rep("1", length(chevPos)))
chevDists <- outer(chevPos, chevPos,
                   FUN = function(x,y) abs(x-y))
chevGenome <- makeGenome(list(chevPos),
                         alleles = rep(list(c("A","a")), length(chr)),
                         chromosome = chr,
                         markerPureDom(length(chr)))
chevLT <- matrix(c(0.5, 0.82, 0.64, 0.48, 0.32, 0.12, 0.08, 0.02,
                   0, 0.5, 0.65, 0.51, 0.35, 0.12, 0.1, 0.02,
                   0, 0, 0.5, 0.75, 0.57, 0.27, 0.15, 0.06,
                   0, 0, 0, 0.5, 0.75, 0.45, 0.15, 0.03,
                   0, 0, 0, 0, 0.5, 0.6, 0.19, 0.09,
                   0, 0, 0, 0, 0, 0.5, 0.33, 0.35,
                   0, 0, 0, 0, 0, 0, 0.5, 0.61,
                   0, 0, 0, 0, 0, 0, 0, 0.5),
                 nrow = length(chr))
chevCorr <- chevLT + t(chevLT)
chevCorrTheory <- theoryCorrelation(chevGenome)
chevCorrDiff <- chevCorr - chevCorrTheory
png("chevCorr.png")
par(mar = c(0.1, 0.1, 0.1, 0.1))
## experimental output
corrImg(chevCorr,
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(41),
        breaks = seq(-1, 1, length.out = 42), axes = FALSE,
        main = "", xlab = "",ylab = "")
text(x = rep(xpos, times = length(chevPos)),
     y = rep(xpos, each = length(chevPos)),
     labels = round(t(apply(chevCorr, 1, rev)), 2))
dev.off()
## theoretical output
png("chevCorrTheory.png")
par(mar = c(0.1, 0.1, 0.1, 0.1))
corrImg(chevCorrTheory,
        col = colorRampPalette(c("steelblue", "white", "firebrick"))(41),
        breaks = seq(-1, 1, length.out = 42),
        axes = FALSE, main = "", xlab = "",
        ylab = "")
text(x = rep(xpos, times = length(chevPos)),
     y = rep(xpos, each = length(chevPos)),
     labels = round(t(apply(chevCorrTheory, 1, rev)), 2))
dev.off()
## difference output
png(paste(imgDir, "chevCorrDiff.png", sep = "/"))
par(mar = c(0.1, 0.1, 0.1, 0.1))
corrImg(chevCorrDiff, xaxt = "n", yaxt = "n", main = "", xlab = "",
        ylab = "", breaks = seq(-0.5, 0.5, length.out = 12),
        col = colorRampPalette(c("steelblue", "white",
                                 "firebrick"))(11))
text(x = rep(xpos, times = length(chevPos)),
     y = rep(xpos, each = length(chevPos)),
     labels = round(t(apply(chevCorrDiff, 1, rev)), 2))
dev.off()

## but is this random? perform a lineup test
truPos <- sample(1:25, 1)
npop <- 510
mother <- with(chevGenome,
               makeGenome(location, alleles, chromosome,
                          markerPureDom(length(chromosome))))
father <- with(chevGenome,
               makeGenome(location, alleles, chromosome,
                          markerPureRec(length(chromosome))))
F1 <- sex(mother, father)
par(mfrow = c(5,5), mar = c(1, 0.1, 0.1, 0.1))
for (ii in 1:25) {
    if (ii == truPos) {
        corrImg(chevCorr - chevCorrTheory, xaxt = "n", yaxt = "n",
                main = "", xlab = "", ylab = "",
                breaks = seq(-0.5, 0.5, length.out = 12),
                col = colorRampPalette(c("steelblue", "white",
                                         "firebrick"))(11))
        mtext(ii, side = 1, cex = 0.8)
    } else {
        tempPop <- asPopulation(replicate(npop, sex(F1, F1),
                                          simplify = FALSE))
        corrImg(popCorrelation(tempPop) - chevCorrTheory,
                xaxt = "n", yaxt = "n", main = "", xlab = "",
                ylab = "", breaks = seq(-0.5, 0.5, length.out = 12),
                col = colorRampPalette(c("steelblue", "white",
                                         "firebrick"))(11))
        mtext(ii, side = 1, cex = 0.8)
    }
}

## alternatively, sample more and colour by quantile
nsim <- 10000
simCorrs <- replicate(nsim,
                      popCorrelation(
                          asPopulation(
                              replicate(npop,
                                        sex(F1, F1),
                                        simplify = FALSE))))
quants <- apply(simCorrs, c(1,2), quantile, probs = seq(0,1,0.025))
lessthan <- sapply(1:nsim, function(ii) simCorrs[,,ii] <= chevCorr)
quantiles <- matrix(apply(lessthan, 1, sum), nrow = 8)

corrImg(quantiles - diag(5000, nrow = 8, ncol = 8),
        col = colorRampPalette(c("steelblue", "white",
                                 "firebrick"))(3),
                                        #breaks = seq(0,10000,length.out = 20))
        breaks = c(0, 250, 9750, 10000))
text(x = rep(xpos, times = length(chevPos)),
     y = rep(xpos, each = length(chevPos)),
     labels = round(t(apply(quantiles/10000, 1, rev)), 2))

quantilePal <- colorRampPalette(c("steelblue", "white",
                                  "firebrick"))(3)
quantileCuts <- matrix(as.numeric(cut(quantiles, c(-1, 250, 9750, 10000))),
                       nrow = 8)
png("chevCorrTest.png", width = 720, height = 720, type = "cairo")
markerNames <- c("D1Mit3", "D1Mit20", "D1Mit74", "D1Mit7", "D1Mit11",
                 "D1Mit14", "D1Mit17", "D1Mit155")
par(mfrow = c(8,8), mar = c(0.1,0.1,0.1,0.1))
for (ii in 1:8) {
    for (jj in 1:8) {
        if (ii == jj) {
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            text(0.5, 0.5, markerNames[ii], cex = 1.5)
        } else if (ii < jj) {
            par(mar = c(2.1,0.5,0.5,0.5))
            tempdens <- density(simCorrs[ii,jj,])
            plot(NA, xlim = range(tempdens$x), ylim = range(tempdens$y),
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            axis(1, at = seq(min(tempdens$x), max(tempdens$x),
                             length.out = 5),
                 labels = round(seq(min(tempdens$x), max(tempdens$x),
                                    length.out = 5),2))
            polygon(tempdens, col = "gray70")
            abline(v = chevCorr[ii,jj], col = "firebrick", lwd = 2)
            abline(v = chevCorrTheory[ii,jj], col = "black", lwd = 2)
            #if (chevCorr[ii,jj] <= chevCorrTheory[ii,jj]) {
                shadInd <- tempdens$x <= chevCorr[ii,jj]
            #} else {
            #    shadInd <- tempdens$x >= chevCorr[ii,jj]
            #}
            polygon(c(tempdens$x[shadInd], chevCorr[ii,jj]),
                    c(tempdens$y[shadInd], 0),
                    col = "gray50")
        } else {
            par(mar = c(0.1,0.1,0.1,0.1))
            plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
            rect(0, 0, 1, 1, col = quantilePal[quantileCuts[ii,jj]])
            text(0.5, 0.5, quantiles[ii,jj], cex = 2)
        }
    }
}
dev.off()
