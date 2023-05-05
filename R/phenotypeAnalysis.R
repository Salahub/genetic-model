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

## get snp column indices
snpCols <- (ncol(pheno)+1):(ncol(phenoGeno))
## categorical target: chi square test of independence
tests <- lapply(snpCols,
               function(ind) chisq.test(phenoGeno[, "coat"],
                                        phenoGeno[, ind]))
## extract p-values
pvals <- sapply(tests, function(tst) tst$p.value)

## select one SNP from each request at random
indsByReq <- split(snps[nlev[-(1:ncol(pheno))] != 1, "rs"],
                   snps[nlev[-(1:ncol(pheno))] != 1, "requested"])


## paigen data
pheno <- read.csv("./data/Paigen4_animaldata.csv")
## add genotypes

## do the Kruskal-Wallis association test between a targets and snps
tests <- lapply(snpCols,
                function(ind) kruskal.test(phenoGeno[, target],
                                           phenoGeno[, ind]))
