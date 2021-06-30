## load functions, set image output directory
source("simulationFunctions.R")
imgDir <- "../img/"

## SIMULATIONS #######################################################
## start incredibly simple: a single chromosome, markers equidistant
nchrom <- 1
nmark <- 20
dists <- list(c(rep(15,19)))
mother <- abiogenesis(nchrom, nmark, dists = dists, allele = 1)
father <- abiogenesis(nchrom, nmark, dists = dists, allele = 0)
F1 <- sex(mother, father) # this always gives the same genome
npop <- 1e5
interCross1 <- replicate(npop, sex(F1, F1), simplify = FALSE)
backCross1 <- replicate(npop, sex(F1, mother), simplify = FALSE)
interCor1 <- popCorrelation(interCross1)
backCor1 <- popCorrelation(backCross1)
png(paste(imgDir, "inter1.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(interCor1, xaxt = "n", yaxt = "n", main = "")
dev.off() # for output
png(paste(imgDir, "back1.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(backCor1, xaxt = "n", yaxt = "n", main = "")
dev.off() # for output

## now consider the case where markers are not equidistant
dists <- list(c(rep(2,5), rep(5,5), rep(10,5),
                rep(20,4))) # and varied dists
mother <- abiogenesis(nchrom, nmark, dists = dists, allele = 1)
father <- abiogenesis(nchrom, nmark, dists = dists, allele = 0)
F1 <- sex(mother, father)
interCross2 <- replicate(npop, sex(F1, F1), simplify = FALSE)
backCross2 <- replicate(npop, sex(F1, mother), simplify = FALSE)
interCor2 <- popCorrelation(interCross2)
backCor2 <- popCorrelation(backCross2)
png(paste(imgDir, "inter2.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(interCor2, xaxt = "n", yaxt = "n", main = "")
dev.off() # for output
png(paste(imgDir, "back2.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(backCor2, xaxt = "n", yaxt = "n", main = "")
dev.off() # for output

## now multiple chromosomes, but keep distances constant
nchrom <- 2
nmark <- c(8, 12)
dists <- list(rep(5, 7), rep(15, 11))
mother <- abiogenesis(nchrom, nmark, dists = dists, allele = 1)
father <- abiogenesis(nchrom, nmark, dists = dists, allele = 0)
F1 <- sex(mother, father)
interCross3 <- replicate(npop, sex(F1, F1), simplify = FALSE)
backCross3 <- replicate(npop, sex(F1, mother), simplify = FALSE)
interCor3 <- popCorrelation(interCross3)
backCor3 <- popCorrelation(backCross3)
png(paste(imgDir, "inter3.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(interCor3, xaxt = "n", yaxt = "n", main = "")
dev.off() # for output
png(paste(imgDir, "back3.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(backCor3, xaxt = "n", yaxt = "n", main = "")
dev.off() # for output

## SIMULATIONS EMULATING OTHER WORK ##################################
## Cheverud 2001: factorial design with one chromosome and different
## distances, this is no different than the first case above, but
## insert his particulars
npop <- 500
nchrom <- 1
nmark <- 16
dists <- list(rep(6.25, 16))
mother <- abiogenesis(nchrom, nmark, dists = dists, allele = 1)
father <- abiogenesis(nchrom, nmark, dists = dists, allele = 0)
F1 <- sex(mother, father)
cheverud <- replicate(npop, sex(F1, F1), simplify = FALSE)
cheverudCor <- popCorrelation(cheverud)
png(paste(imgDir, "cheverud.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(cheverudCor, xaxt = "n", yaxt = "n", main = "")
dev.off()

## Landner and Botstein 1989: 12 chromosomes, each 100 cM long with
## markers placed along each chromosome at 20 cM intervals, backcross
npop <- 250
nchrom <- 12
nmark <- rep(6, nchrom)
dists <- lapply(nmark, function(x) rep(20, times = x-1))
mother <- abiogenesis(nchrom, nmark, dists = dists, allele = 1)
father <- abiogenesis(nchrom, nmark, dists = dists, allele = 0)
F1 <- sex(mother, father)
landbot <- replicate(npop, sex(F1, mother), simplify = FALSE)
landbotCor <- popCorrelation(landbot)
png(paste(imgDir, "landbot.png", sep = "/")) # for output
par(mar = c(0.1,0.1,0.1,0.1)) # if no title is desired
corrImg(landbotCor)
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
nchrom <- 10
nmark <- rep(5, nchrom)
dists <- lapply(nmark, function(x) rep(invHaldane(pr),
                                       times = x - 1))
mother <- abiogenesis(nchrom, nmark, dists = dists, allele = 1)
father <- abiogenesis(nchrom, nmark, dists = dists, allele = 0)
F1 <- sex(mother, father)
liji <- replicate(npop, sex(F1, F1), simplify = FALSE)
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
nmark <- length(chevPos)
xpos <- seq(0, 1, length.out = nmark)
chevDists <- outer(chevPos, chevPos,
                   FUN = function(x,y) abs(x-y))
chevLT <- matrix(c(0.5, 0.82, 0.64, 0.48, 0.32, 0.12, 0.08, 0.02,
                   0, 0.5, 0.65, 0.51, 0.35, 0.12, 0.1, 0.02,
                   0, 0, 0.5, 0.75, 0.57, 0.27, 0.15, 0.06,
                   0, 0, 0, 0.5, 0.75, 0.45, 0.15, 0.03,
                   0, 0, 0, 0, 0.5, 0.6, 0.19, 0.09,
                   0, 0, 0, 0, 0, 0.5, 0.33, 0.35,
                   0, 0, 0, 0, 0, 0, 0.5, 0.61,
                   0, 0, 0, 0, 0, 0, 0, 0.5),
                 nrow = nmark)
chevCorr <- chevLT + t(chevLT)
chevCorrTheory <- theoryCorrelation(chevDists)
chevCorrDiff <- chevCorr - chevCorrTheory
png(paste(imgDir, "chevCorr.png", sep = "/"))
par(mar = c(0.1, 0.1, 0.1, 0.1))
## experimental output
corrImg(chevCorr, xaxt = "n", yaxt = "n",
        main = "", xlab = "",ylab = "")
text(x = rep(xpos, times = nmark), y = rep(xpos, each = nmark),
     labels = round(t(apply(chevCorr, 1, rev)), 2))
dev.off()
## theoretical output
png(paste(imgDir, "chevCorrTheory.png", sep = "/"))
par(mar = c(0.1, 0.1, 0.1, 0.1))
corrImg(chevCorrTheory, xaxt = "n", yaxt = "n", main = "", xlab = "",
        ylab = "")
text(x = rep(xpos, times = nmark), y = rep(xpos, each = nmark),
     labels = round(t(apply(chevCorrTheory, 1, rev)), 2))
dev.off()
## difference ouput
png(paste(imgDir, "chevCorrDiff.png", sep = "/"))
par(mar = c(0.1, 0.1, 0.1, 0.1))
corrImg(chevCorrDiff, xaxt = "n", yaxt = "n", main = "", xlab = "",
        ylab = "", breaks = seq(-0.5, 0.5, length.out = 12),
        col = colorRampPalette(c("steelblue", "white",
                                 "firebrick"))(11))
text(x = rep(xpos, times = nmark), y = rep(xpos, each = nmark),
     labels = round(t(apply(chevCorrDiff, 1, rev)), 2))
dev.off()

## but is this random? perform a lineup test
truPos <- sample(1:25, 1)
npop <- 510
mother <- abiogenesis(1, nmark, dists = list(diff(chevPos)),
                      allele = 1)
father <- abiogenesis(1, nmark, dists = list(diff(chevPos)),
                      allele = 0)
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
        tempPop <- replicate(npop, sex(F1, F1), simplify = FALSE)
        corrImg(popCorrelation(tempPop) - chevCorrTheory,
                xaxt = "n", yaxt = "n", main = "", xlab = "",
                ylab = "", breaks = seq(-0.5, 0.5, length.out = 12),
                col = colorRampPalette(c("steelblue", "white",
                                         "firebrick"))(11))
        mtext(ii, side = 1, cex = 0.8)
    }
}

## alternatively, sample more and colour by quantile
nsim <- 1000
simCorrs <- replicate(nsim,
                      popCorrelation(replicate(npop,
                                               sex(F1, F1),
                                               simplify = FALSE)))
lessthan <- sapply(1:nsim, function(ii) simCorrs[,,ii] <= chevCorr)
quantiles <- matrix(apply(test, 1, sum), nrow = 8)

corrImg(quantiles - diag(500, nrow = 8, ncol = 8),
        col = colorRampPalette(c("steelblue", "white",
                                 "firebrick"))(9),
        breaks = seq(0,1000,length.out = 10))
