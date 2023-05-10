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
testSnps <- pGbyRsub
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
abline(h = log(0.05), lty = 2, col = "firebrick")
hist(pvals)
## seemingly robust to random subsamples!


## 2nd example: paigen data ##########################################
pheno <- read.csv("./data/Paigen4_animaldata.csv")

## join to snps
phenoInds <- 1:7
phenoGeno <- joinPhenoGeno(pheno, snps, phenoCols = 1:7)
## filter out SNPs by distributions
nlev <- sapply(phenoGeno, function(col) sum(levels(col) != ""))
phenoGeno <- phenoGeno[,nlev != 1]
## inds to aid in selection
colChr <- snps[names(phenoGeno)[-phenoInds], "chr"]
colRd <- snps[names(phenoGeno)[-phenoInds], "requested"]
tabs <- sapply(phenoGeno, table)

## a few others: most complete per request
set.seed(15372023)
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
testSnps <- pGbyChr
snpCols <- (ncol(pheno)+1):(ncol(testSnps))
## continuous target: kruskal-wallis one-way anova test
tests <- lapply(snpCols,
                function(ind) kruskal.test(testSnps[, "TG"],
                                           g = testSnps[, ind]))
## extract p-values
pvals <- sapply(tests, function(tst) tst$p.value)
## get a sequence of kappas for p-values
kseq <- exp(seq(-8, 8, by = 0.01))
pooled <- sapply(kseq, poolChi, p = pvals)
plot(log(kseq), log(pooled), type = 'l') # kind of cool.. clear min
hist(pvals)

## smallest p-values are for chromosomes 1, 4, 10, 19
chr1 <- pGbyR[, c(phenoInds, which(names(pGbyR) %in%
                                   snps[snps[, "chr"] == "19", "rs"]))]
snpCols <- (ncol(pheno)+1):(ncol(chr1))
tests <- lapply(snpCols,
                function(ind) kruskal.test(chr1[, "TG"],
                                           g = chr1[, ind]))
pvals <- sapply(tests, function(tst) tst$p.value)
kseq <- exp(seq(-8, 8, by = 0.01))
pooled <- sapply(kseq, poolChi, p = pvals)
plot(log(kseq), log(pooled), type = 'l')
hist(pvals)

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
