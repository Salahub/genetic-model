### Helper Functions for Genome Objects ##############################

##' @title A helper to generate a genetic sequence
##' @param nchrom the integer number of chromosomes in the sequence
##' @param nmark a vector of integers of length nchrom giving the
##' number of markers for each chromosome
##' @param alleles (optional) vector of possible encodings for the
##' markers
##' @param markerFun a function which accepts an integer and vector
##' of alleles and returns a 2-column matrix of marker measurements
##' for a single chromosome to be mapped on all chromosomes
##' @return a list of 2-column matrices giving marker encodings on
##' each chromosome
##' @author Chris Salahub
generateEncoding <- function(nchrom, nmark, alleles = c(1, 1),
                             markerFun = markerPure) {
    ## map the single chromosome function across all chromosomes
    lapply(1:nchrom, function(ind) markerFun(nmark[ind], alleles))
}

##' @title Marker measurement generation functions
##' @param nmark an integer giving the number of markers
##' @param alleles (optional) vector of possible marker measurements
##' @return a 2-column matrix with nmark rows with entries from
##' alleles
##' @author Chris Salahub
markerPure <- function(nmark, alleles) { # pure single chromosome
    matrix(rep(allele, nmark*2), ncol = 2)
}
markerHybrid <- function(nmark, alleles) { # mixed (i.e. Aa on all)
    cbind(rep(1, nmark), rep(0, nmark))
}
markerIntercross <- function(nmark, alleles) { # from an intercross
    allele <- as.numeric(runif(2) < 0.5)
    cbind(rep(allele[1], nmark), rep(allele[2], nmark))
}
abio.back <- function(nmark, alleles) { # no drift
    cbind(rep(allele, nmark), rep(runif(1) < 0.5, nmark))
}
abio.random <- function(nmark, alleles) { # totally random genome
    matrix(as.numeric(runif(2*nmark) < 0.5), ncol = 2)
}


## generates distances using the specified generation function
distGenesis <- function(nchrom, nmark, distFun) {
    lapply(1:nchrom, function(ind) distFun(nmark[ind]-1))
}

## haldane's map function
mapHaldane <- function(d) {
    p <- (1 - exp(-d/50))/2
    p
}
## its inverse
invHaldane <- function(pr) {
    d <- -50*log(1 - 2*pr)
    d
}

## kosambi's map function
mapKosambi <- function(d) {
    p <- (exp(d/25) - 1)/(2*(exp(d/25) + 1))
    p
}
## its inverse
invKosambi <- function(pr) {
    d <- 25*log((1 + 2*pr)/(1 - 2*pr))
    d
}
