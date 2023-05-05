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
                                         testSnps[, ind])
                                         #simulate.p.value = TRUE,
                                         #B = 10000))
## extract p-values
pvals <- sapply(tests, function(tst) tst$p.value)
## get a squence of kappas for p-values
kseq <- exp(seq(-8, 8, by = 0.01))
pooled <- sapply(kseq, poolChi, p = pvals)
plot(log(kseq), log(pooled), type = 'l') # kind of cool.. clear min

## does this min have something to do with beta densities
step <- 0.1
a <- seq(0.3, 5, by = step)
b <- seq(0.3, 5, by = step)
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

## 2nd example: paigen data ##########################################
pheno <- read.csv("./data/Paigen4_animaldata.csv")
## add genotypes

## do the Kruskal-Wallis association test between a targets and snps
tests <- lapply(snpCols,
                function(ind) kruskal.test(phenoGeno[, target],
                                           phenoGeno[, ind]))
