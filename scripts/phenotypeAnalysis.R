library(toyGenomeGenR)

## FUNCTIONS #########################################################
narrowPlot <- function(xgrid, ygrid, main = "", xlab = "", ylab = "",
                       xticks = xgrid, yticks = ygrid,
                       mars = c(2.1, 2.1, 1.1, 1.1),
                       xlim = range(xgrid), ylim = range(ygrid),
                       addGrid = TRUE, ...) {
    par(mar = mars) # set narrow margins
    plot(NA, ylim = ylim, xlim = xlim, xaxt = 'n', xlab = "",
         yaxt = 'n', ylab = "", main = "", ...)
    ## add labels
    mtext(main, side = 3, line = 0, cex = 0.8) # main
    mtext(ylab, side = 2, line = 1, cex = 0.8) # ylab
    mtext(xlab, side = 1, line = 1, padj = 0, cex = 0.8) # xlab
    ## add grid lines
    if (addGrid) {
        abline(h = ygrid, v = xgrid, lty = 1,
               col = adjustcolor("gray", alpha.f = 0.4))
    }
    ## and ticks
    mtext(side = 1, at = xgrid, text = "|", line = 0, cex = 0.5,
          padj = -2)
    mtext(text = xticks, at = xgrid, side = 1, cex = 0.8)
    mtext(side = 2, at = ygrid, text = "|", line = 0, cex = 0.5,
          padj = 1)
    mtext(text = yticks, at = ygrid, side = 2, cex = 0.8)
}

## a function to join genotypes and phenotypes
joinPhenoGeno <- function(pheno, geno, phenoCols = 1:4) {
    strains <- factor(pheno$strain) # for easy reference
    geno2pheno <- sapply(levels(strains), # index match back
                         function(str) which(grepl(paste0("^", str,
                                                          "(?=\\:)"),
                                                   colnames(snps),
                                                   perl = TRUE)))
    ## reform genotype with the correct ordering
    geno <- matrix(nrow = nrow(pheno), ncol = nrow(snps))
    colnames(geno) <- snps[, "rs"] # rs SNP names
    for (ii in seq_len(nrow(pheno))) { # fill matrix
        if (length(geno2pheno[[pheno$strain[ii]]]) > 0) {
            geno[ii, ] <- snps[, geno2pheno[[pheno$strain[ii]]]]
        }
    }
    ## drop empty genos
    empty <- sapply(geno2pheno, length) == 0
    ## bind together based on phenoCols
    phenoGeno <- cbind(pheno[!empty, phenoCols],
                       as.data.frame(geno[!empty,],
                                     stringsAsFactors = TRUE))
    phenoGeno
}

## independent pooled chi value
poolChi <- function(p, kap) {
    pchisq(sum(qchisq(p, df = kap, lower.tail = FALSE)),
           df = length(p)*kap, lower.tail = FALSE)
}

## try the method of moments gamma match from ferrari's arxiv note on
## chi-squared sums
gammaApprox <- function(qs, sigma, kap) {
    diag(sigma) <- 0 # for later sum
    M <- length(qs)
    u <- 2*(1 + 2*sum(sigma*kap)/(M*kap))
    pgamma(sum(qs), shape = (M*kap)/u, scale = u,
           lower.tail = FALSE)
}
## can also try a satterthwaite approximation
satterApprox <- function(qs, sigma, kap) {
    diag(sigma) <- 0
    M <- length(qs)
    ex <- M*kap
    varx <- 2*M*kap + 2*kap*sum(sigma)
    f <- 2*ex^2/varx
    c <- varx/(2*ex)
    pchisq(sum(qs)/c, df = f, lower.tail = FALSE)
}
## write a wrapper that incorporates both
poolChiDep <- function(ps, kap, sigma,
                       method = c("gamma", "satterthwaite")) {
    chis <- qchisq(ps, df = kap, lower.tail = FALSE)
    method <- match.arg(method)
    if (method == "gamma") {
        gammaApprox(chis, sigma, kap)
    } else if (method == "satterthwaite") {
        satterApprox(chis, sigma, kap)
    } else {
        stop("Method must be one of 'gamma' or 'satterthwaite'")
    }
}

## and a bin splitting p-value
## this is a faster implementation of the recursive binning of marbR
## written custom for this application
unirandsplit <- function(x, g, n, lim = 10) {
    brk1 <- sample(lim:(n-lim), 1)
    if (brk1-2*lim <= 0) { # maintain minimum bin size
        brk2 <- sample(seq(brk1+lim, n-lim, by = 1), 1)
        brks <- c(0, brk1, brk2, n)
    } else if (n-brk1-2*lim <= 0) {
        brk2 <- sample(seq(lim, brk1-lim, by = 1), 1)
        brks <- c(0, brk2, brk1, n)
    } else {
        brk2 <- sample(c(seq(lim, brk1-lim, by = 1),
                         seq(brk1+lim, n-lim, by = 1)), 1)
        brks <- c(0, min(brk1, brk2), max(brk1, brk2), n)
    }
    xc <- cut(x, breaks = c(0, brk1, brk2, n))
    obs1 <- table(xc, g) # counts
    ex1 <- outer(diff(brks), table(g)/n) # expected densities
    pchisq(sum(obs1^2/ex1) - n,
           df = prod(dim(obs1) - 1), lower.tail = FALSE)
}


##' load the data necessary to fit models that convert the observed
##' correlations between markers into the correlation between
##' chi-transformed random variables so the Satterthwaite
##' approximation can be used
kseq <- exp(seq(-8, 8, by = 0.1))
chiCordf <- readRDS("chiCorrs.Rds")
## fit conversion models
chiCorMods <- lapply(log(kseq), function(k) {
    lm(chicor ~ I(zcor^2) + I(zcor^4) + I(zcor^6) + I(zcor^8) +
           I(zcor^10), data = chiCordf,
       subset = abs(chiCordf$logkap - k)<0.0001)#$coefficients
})
## a wrapper function which uses these models to convert values
convertRho <- function(sigma, kapInd, models = chiCorMods) {
    matrix(predict(chiCorMods[[kapInd]],
                   newdata = data.frame(zcor = c(sigma))),
           ncol = ncol(sigma))
}


## load strain SNPs
snps <- readRDS("./data/strainSNPs.Rds")
snps[snps == ""] <- NA
snps <- snps[!is.na(snps[, "rs"]),] # remove missing names
rownames(snps) <- snps[, "rs"] # index for quick search later
snps <- snps[order(as.numeric(snps[, "cMs"])),] # order by cm dists
snps <- snps[order(snps[, "chr"]),] # return chromosome order


## EXAMPLE: PAIGEN DATA ##############################################
## read in phenotype data
pheno <- read.csv("./data/Paigen4_animaldata.csv")

## join to snps
phenoInds <- 1:7
phenoGeno <- joinPhenoGeno(pheno, snps, phenoCols = 1:7)
## filter out SNPs by distributions
nlev <- sapply(phenoGeno, function(col) sum(levels(col) != ""))
phenoGeno <- phenoGeno[, nlev != 1]
## indexes to aid in later selection
colChr <- snps[names(phenoGeno)[-phenoInds], "chr"]
colRd <- snps[names(phenoGeno)[-phenoInds], "requested"]
tabs <- sapply(phenoGeno, table)

##' aim for data which is not identical for all individuals in the
##' data by looking at which markers are least uniform by API
##' request ("read")
set.seed(15372023)
byRead <- sapply(split(sapply(tabs, sum)[-phenoInds], colRd),
                 function(k) names(k)[which.max(k)])
pGbyR <- phenoGeno[, c(names(phenoGeno)[phenoInds], byRead)]
## select a random subset of these
pGbyRsub <- cbind(pGbyR[, phenoInds],
                  pGbyR[ , sample((max(phenoInds)+1):ncol(pGbyR),
                                  40)])

##' perform the same analysis to get the most varied and complete by
##' chromosome, select the top k
k <- 1 # for the paper, one was chosen to get independent markers
byChr <- sapply(split(sapply(tabs, sum)[-phenoInds], colChr),
                function(k) names(k)[order(k, decreasing = TRUE)[1:1]])
pGbyChr <- phenoGeno[, c(names(phenoGeno)[phenoInds], byChr)]

## use these to get snp column indices (female only)
testSnps <- pGbyR[pGbyR$sex == "f",]
snpCols <- (ncol(pheno)+1):(ncol(testSnps))
## order to restore chromosome ordering and locations
snpOrder <- rownames(snps)[rownames(snps) %in%
                           names(testSnps)[snpCols]]
testSnps <- testSnps[, c(names(testSnps)[1:ncol(pheno)],
                         snpOrder)]
## apply random splitting to obtain a p-value
tests <- lapply(snpCols,
                function(ind) unirandsplit(rank(testSnps[, "nonHDL"]),
                                           g = testSnps[, ind],
                                           n = nrow(testSnps)))
## extract p-values
pvals <- unlist(tests)

##' next: kappa values are swept in the chi-squared pooled p-value
##' to obtain a sequence of kappas which are most powerful for the
##' obtained p-values
nullQuants <- readRDS("curveMinQuantiles.Rds") # comparison values
quantLevs <- c("5%") # quantiles
kseq <- exp(seq(-8, 8, by = 0.1)) # sweep values
pooled <- sapply(kseq, poolChi, p = pvals) # pooled p-values
## not included in the thesis: unadjusted curves of the pooled
## p-values by kappa
png("FemPaigenCurve.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(-3, 3, by = 1.5),
           ygrid = seq(-80, 0, by = 20),
           xlim = c(-3.5, 3.5),
           xlab = expression(log[10]~{"("~kappa~")"}),
           ylab = expression(log[10]~{"("~p~")"}))
lines(log(kseq, 10), log(pooled, 10), type = 'l') # clear min
abline(h = log(nullQuants[quantLevs, as.character(1000)], 10),
       lty = 2, col = "firebrick")
text(x = rep(3.95, 3), labels = quantLevs, cex = 0.6, xpd = NA,
     y = log(nullQuants[quantLevs, as.character(100)], 10),
     adj = c(0.2, 0.5))
dev.off()

##' now: perform the Satterthwaite adjustment
##' compute the correlation matrix
obscorrs <- cor(sapply(testSnps[, snpCols], as.numeric),
                use = "pairwise.complete.obs")
chrs <- snps[names(testSnps)[snpCols], "chr"]
locs <- split(as.numeric(snps[names(testSnps)[snpCols], "cMs"]),
              chrs)
thcorrs <- theoryCorrelation(list(location = locs, chromosome = chrs))
## apply the adjusted p-value function
pooledadj <- sapply(1:length(kseq),
                    function(ii) {
                        poolChiDep(ps = pvals, kap = kseq[ii],
                                   sigma = convertRho(thcorrs,
                                                      kapInd = ii))
                    })

## the adjusted curve corresponds with Figure 6.9
png("FemPaigenCurveAdj2.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(-3, 3, by = 1.5),
           ygrid = seq(-4, 0, by = 1),
           xlim = c(-3.5, 3.5),
           xlab = expression(log[10]~{"("~kappa~")"}),
           ylab = expression(log[10]~{"("~p~")"}))
lines(log(kseq, 10), log(pooledadj, 10), type = 'l') # clear min
abline(h = log(nullQuants[quantLevs, as.character(1000)], 10),
       lty = 2, col = "firebrick")
text(x = rep(3.95, 3), labels = quantLevs, cex = 0.6, xpd = NA,
     y = log(nullQuants[quantLevs, as.character(100)], 10),
     adj = c(0.2, 0.5))
dev.off()

## plot chi 0.18 quantile transformations of p-values (Figure 6.10(b))
png("FemPaigenPvalCutoff.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(0, 1, by = 0.2),
           ygrid = seq(0, 50, by = 10), ylim = c(0, 55),
           ylab = expression(F[chi]^{-1}~{"("~1-p~"; 0.18)"}),
           xlab = "Sample quantile of p")
points(x = ppoints(length(pvals)), cex = 0.6,
       sort(qchisq(pvals, 0.18, lower.tail = FALSE),
            decreasing = TRUE))
abline(v = 0.05, lty = 2, col = "firebrick")
dev.off()

## select and inspect the top 5% of the quantile transformed p-values
pvalOrd <- order(qchisq(pvals, 0.18, lower.tail = FALSE),
                 decreasing = TRUE)
topSnps <- names(testSnps)[snpCols][pvalOrd[1:53]]
## where on the genome are these
png("FemPaigenTable.png", width = 6, height = 3, units = "in",
    res = 480)
dpar(mar = c(4.1, 4.1, 1.1, 1.1))
plot(table(factor(snps[topSnps, "chr"],
                  levels = c(as.character(1:19), "X"))),
     xlab = "Chromosome", ylab = "Frequency")
dev.off()
## compare with paigen data
png("FemPaigenTablePaper.png", width = 6, height = 3, units = "in",
    res = 480)
par(mar = c(4.1, 4.1, 1.1, 1.1))
plot(table(factor(c(rep("1", 4), rep("2", 2), rep("3", 1),
                    rep("4", 3), rep("5", 3), rep("6", 2),
                    rep("7", 3), rep("8", 3), rep("9", 3),
                    rep("10", 1), rep("11", 3), rep("12", 2),
                    rep("13", 1), rep("14", 1), rep("15", 3),
                    rep("16", 0), rep("17", 2), rep("18", 1),
                    rep("19", 4), rep("X", 1)),
                  levels = c(as.character(1:19), "X"))),
     xlab = "Chromosome", ylab = "Frequency")
dev.off()


## the next examples are not included in the thesis: the example above
## was not cherry picked, but it was realized a more complete
## exploration of this first example would be more meaningful than
## the scattershot exploration of several
## example: coat colour data #########################################
## coat data
pheno <- read.csv("./data/strainCoats.csv")

## join to snps
phenoInds <- 1:2
phenoGeno <- joinPhenoGeno(pheno, snps, phenoCols = phenoInds)
## filter out SNPs by distributions
nlev <- sapply(phenoGeno, function(col) sum(levels(col) != ""))
phenoGeno <- phenoGeno[,nlev != 1]
## features to help with genome selection
tabs <- sapply(phenoGeno, table)
colChr <- snps[names(phenoGeno)[-phenoInds], "chr"]
colRd <- snps[names(phenoGeno)[-phenoInds], "requested"]

## a few others: most complete per request
set.seed(151008085)
byRead <- sapply(split(sapply(tabs, sum)[-phenoInds], colRd),
                 function(k) names(k)[which.max(k)])
pGbyR <- phenoGeno[, c(names(phenoGeno)[phenoInds], byRead)]
## random subset of these
pGbyRsub <- cbind(pGbyR[, phenoInds],
                  pGbyR[ , sample((max(phenoInds)+1):ncol(pGbyR),
                                  40)])
## three most complete by chromosome
byChr <- sapply(split(sapply(tabs, sum)[-phenoInds], colChr),
                function(k) names(k)[order(k, decreasing = TRUE)[1:3]])
pGbyChr <- phenoGeno[, c(names(phenoGeno)[phenoInds],
                         byChr)]

## placeholder to easily change target
testSnps <- pGbyR
testSnps$coat <- grepl("agouti", testSnps$coat)
snpCols <- (ncol(pheno)+1):(ncol(testSnps))
## categorical target: chi square test of independence
tests <- lapply(snpCols,
                function(ind) chisq.test(testSnps[, "coat"],
                                         testSnps[, ind]))
                                         #simulate.p.value = TRUE,
                                         #B = 10000))
## extract p-values
pvals <- sapply(tests, function(tst) tst$p.value)
## get a squence of kappas for p-values
kseq <- exp(seq(-8, 8, by = 0.01))
pooled <- sapply(kseq, poolChi, p = pvals)
png("coatCurve.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(-3, 3, by = 1.5),
           ygrid = seq(-250, 0, by = 50),
           xlim = c(-3.5, 3.5),
           xlab = expression(log[10]~{"("~kappa~")"}),
           ylab = expression(log[10]~{"("~p~")"}))
lines(log(kseq, 10), log(pooled, 10), type = 'l') # clear min
abline(h = log(nullQuants[quantLevs, as.character(1000)], 10),
       lty = 2, col = "firebrick")
text(x = rep(3.95, 3), labels = quantLevs, cex = 0.6, xpd = NA,
     y = log(nullQuants[quantLevs, as.character(100)], 10),
     adj = c(0.2, 0.5))
dev.off()
## seemingly robust to random subsamples!

## example 3 hmdp bone density measures on male mice #################
pheno <- read.csv("./data/HMDPpheno1_animaldata.csv")

## join to snps
phenoInds <- 1:6
phenoGeno <- joinPhenoGeno(pheno, snps, phenoCols = 1:6)
## filter out SNPs by distributions
nlev <- sapply(phenoGeno, function(col) sum(levels(col) != ""))
phenoGeno <- phenoGeno[,nlev != 1]
## inds to aid in selection
colChr <- snps[names(phenoGeno)[-phenoInds], "chr"]
colRd <- snps[names(phenoGeno)[-phenoInds], "requested"]
tabs <- sapply(phenoGeno, table)

## a few others: most complete per request
set.seed(834090523)
byRead <- sapply(split(sapply(tabs, sum)[-phenoInds], colRd),
                 function(k) names(k)[which.max(k)])
pGbyR <- phenoGeno[, c(names(phenoGeno)[phenoInds], byRead)]
## random subset of these
pGbyRsub <- cbind(pGbyR[, phenoInds],
                  pGbyR[ , sample((max(phenoInds)+1):ncol(pGbyR),
                                  40)])
## k most complete by chromosome
k <- 1
byChr <- sapply(split(sapply(tabs, sum)[-phenoInds], colChr),
                function(k) names(k)[order(k, decreasing = TRUE)[1:1]])
pGbyChr <- phenoGeno[, c(names(phenoGeno)[phenoInds],
                         byChr)]

## get snp column indices
testSnps <- pGbyRsub
snpCols <- (ncol(pheno)+1):(ncol(testSnps))
## continuous target: kruskal-wallis one-way anova test
tests <- lapply(snpCols,
                function(ind) kruskal.test(testSnps[, "BMD_femur"],
                                           g = testSnps[, ind]))
## extract p-values
pvals <- sapply(tests, function(tst) tst$p.value)
## get a sequence of kappas for p-values
kseq <- exp(seq(-8, 8, by = 0.01))
pooled <- sapply(kseq, poolChi, p = pvals)
plot(log(kseq), log(pooled), type = 'l') # kind of cool.. clear min
hist(pvals)
