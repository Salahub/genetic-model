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

## simple helper
klDivU <- function(tab) {
    tot <- sum(tab)
    unif <- tot/length(tab)
    sum((unif/tot)*log(unif/tab))
}

## chi-square pool function
poolChi <- function(p, k) {
    M <- length(p) # dimension
    pchisq(sum(qchisq(p, df = k, lower.tail = FALSE)), df = M*k,
           lower.tail = FALSE)
}

## load strain SNPs
snps <- readRDS("./data/strainSNPs.Rds")
snps[snps == ""] <- NA
snps <- snps[!is.na(snps[, "rs"]),] # remove missing names
rownames(snps) <- snps[, "rs"] # index for quick search later

## 1st example: coat colour data #####################################
## coat data
pheno <- read.csv("./data/strainCoats.csv")

## join to snps
phenoInds <- 1:2
phenoGeno <- joinPhenoGeno(pheno, snps, phenoCols = phenoInds)
## filter out SNPs by distributions
nlev <- sapply(phenoGeno, function(col) sum(levels(col) != ""))
phenoGeno <- phenoGeno[,nlev != 1]
## features to help with genome selection
colChr <- snps[names(phenoGeno)[-phenoInds], "chr"]
colRd <- snps[names(phenoGeno)[-phenoInds], "requested"]
tabs <- sapply(phenoGeno, table)
## kl divergence
klds <- sapply(tabs, klDivU)
## small kldivs
pGklcutoff <- phenoGeno[, c(phenoInds,
                            which(klds[-phenoInds] < 0.1) + 2)]
## a few others: min per request
minRead <- sapply(split(klds[-phenoInds], colRd),
                  function(k) names(k)[which.min(k)])
pGminR <- phenoGeno[, c(names(phenoGeno)[phenoInds], minRead)]
## min by chromsome
minChr <- sapply(split(klds[-phenoInds], colChr),
                 function(k) names(k)[which.min(k)])
pGminChr <- phenoGeno[, c(names(phenoGeno)[phenoInds],
                          minChr)]

## get snp column indices
testSnps <- pGminR
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
plot(log(kseq), log(pooled), type = 'l') # kind of cool.. clear min

## does this min have something to do with beta densities?
## first assay to get a sense
step <- 0.5
a <- seq(0.5, 5, by = step)
b <- seq(0.5, 5, by = step)
params <- expand.grid(b = b, a = a)
betas <- lapply(seq_len(nrow(params)),
                function(ii) function(n) rbeta(n, params$a[ii],
                                               params$b[ii]))
## generate random p-values
nsim <- 1e4
samps <- lapply(betas, function(f) f(nsim))
## sweep with pooled p kappas
kseq <- exp(seq(-8, 8, by = 0.5))
pools <- lapply(samps, function(smp) sapply(kseq, poolChi, p = smp))

## plot one example
ii <- 13
par(mfrow = c(1,2))
plot(seq(0, 1, 0.01), xlab = "x", ylab = "Density", type = "l",
     main = paste("a =", params$a[ii], ",", "b =", params$b[ii]),
     dbeta(seq(0, 1, 0.01), params$a[ii], params$b[ii]),
     ylim = c(0,5))
plot(log(kseq), log(pools[[ii]]), main = "Chi p-value path",
     xlab = expression(paste0("log(", kappa, ")")),
     ylab = "log(p)", type = "l", ylim = c(-10, 1))
abline(h = log(0.05), lty = 2)

## make this more robust with repeated samples...
simBetaPath <- function(a = 1, b = 1, n = 1e3, nsim = 100,
                        kseq = exp(seq(-8, 8, by = 0.25))) {
    betaSamps <- matrix(rbeta(n*nsim, a, b), ncol = n) # sim matrix
    betaTrans <- lapply(kseq, qchisq, p = betaSamps,
                        lower.tail = FALSE) # quantile values
    paths <- mapply(function(mt, k) pchisq(rowSums(mt), df = k*n,
                                           lower.tail = FALSE),
                    betaTrans, kseq)
    t(paths)
}
## and construct a similar plot
plotBetaPath <- function(a, b, paths = simBetaPath(a, b),
                         kseq = exp(seq(-8, 8, by = 0.25)),
                         ylim1 = c(0,5), ylim2 = c(-10,1)) {
    par(mfrow = c(1,2))
    plot(seq(0, 1, 0.01), xlab = "x", ylab = "Density", type = "l",
         main = paste("Beta(a =", a, ",", "b =", b, ")"),
         dbeta(seq(0, 1, 0.01), a, b), ylim = ylim1)
    plot(NA, xlim = range(log(kseq, base = 10)),
         main = "Chi p-value path",
         xlab = expression(paste(log[10], "(", kappa, ")")),
         ylab = expression(paste(log[10], "(p)")), type = "n",
         ylim = ylim2)
    for (jj in seq_len(ncol(paths))) {
        lines(log(kseq, base = 10), log(paths[,jj], base = 10),
              col = adjustcolor("black", 0.6))
    }
    abline(h = log(0.05, base = 10), lty = 2,
           col = adjustcolor("firebrick", 0.5))
}

## what about a quantile plot?
plotBetaQuants <- function(a, b, qnts = c(0.005, 0.025, 0.25),
                           paths = simBetaPath(a, b),
                           kseq = exp(seq(-8, 8, by = 0.25)),
                           ylim1 = c(0,5), ylim2 = c(-10,1)) {
    par(mfrow = c(1,2))
    plot(seq(0, 1, 0.005), xlab = "x", ylab = "Density", type = "l",
         main = paste("Beta(a =", a, ",", "b =", b, ")"),
         dbeta(seq(0, 1, 0.005), a, b), ylim = ylim1)
    abline(v = seq(0, 1, by = 0.2), h = seq(0, 5, by = 1),
           col = "gray90")
    plot(NA, xlim = range(log(kseq, base = 10)),
         main = "Chi p-value path quantiles",
         xlab = expression(paste(log[10], "(", kappa, ")")),
         ylab = expression(paste(log[10], "(p)")), type = "n",
         ylim = ylim2)
    abline(h = seq(-10, 0, by = 2), v = seq(-3, 3, by = 1),
           col = "gray90")
    for (qn in qnts) {
        polygon(log(c(kseq, rev(kseq)), base = 10),
                c(apply(log(paths, base = 10), 1, quantile, probs = qn),
                  rev(apply(log(paths, base = 10), 1, quantile,
                            probs = 1-qn))),
                col = adjustcolor("gray", 0.25), border = NA)
    }
    lines(x = log(kseq, 10),
          y = apply(log(paths, base = 10), 1, median))
    abline(h = log(0.05, base = 10), lty = 2,
           col = adjustcolor("firebrick", 0.5))
}

## call based on a, b, simulation settings
nsim <- 1e3
n <- 1e3
a <- 1.4
b <- 1.5
kseq <- exp(seq(-8, 8, by = 0.1))
sims <- simBetaPath(a = a, b = b, n = n, nsim = nsim)
plotBetaQuants(a = a, b = b, paths = sims)
##' some interesting settings:
##' a = 0.4, b = 0.2
##' a = 0.82, b = 0.85
##' a = 3, b = 4
##' a = 4, b = 3
##' a = 0.85, b = 0.82
##' a = 0.2, b = 0.4
##' a = 1, b = 2
##' a = 0.9, b = 1.1

## 2nd example: paigen data ##########################################
pheno <- read.csv("./data/Paigen4_animaldata.csv")
## add genotypes

## do the Kruskal-Wallis association test between a targets and snps
tests <- lapply(snpCols,
                function(ind) kruskal.test(phenoGeno[, target],
                                           phenoGeno[, ind]))
