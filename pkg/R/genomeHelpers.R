### Helper Functions for Genome Objects ##############################

##' @title Generate a genetic sequence
##' @param nchrom the integer number of chromosomes in the sequence
##' @param nmark a vector of integers of length nchrom giving the
##' number of markers for each chromosome
##' @param alleles (optional) vector of possible encodings for the
##' markers, assumed to be of length 2
##' @param markerFun a function which accepts an integer and vector
##' of alleles and returns a 2-column matrix of marker measurements
##' for a single chromosome to be mapped on all chromosomes
##' @return a list of 2-column matrices giving marker encodings on
##' each chromosome
##' @author Chris Salahub
generateEncoding <- function(nchrom, nmark, alleles = c(1,0),
                             markerFun = markerPure) {
    ## map the single chromosome function across all chromosomes
    lapply(1:nchrom, function(ind) markerFun(nmark[ind], alleles))
}

##' @title Helpers to simulate marker measurements
##' @param nmark an integer giving the number of markers
##' @param alleles vector of possible marker measurements
##' @return a 2-column matrix with nmark rows with entries from
##' alleles
##' @author Chris Salahub
markerPure <- function(nmark, alleles) {
    matrix(rep(alleles[1], nmark*2), ncol = 2)
}
markerHybrid <- function(nmark, alleles) {
    cbind(rep(alleles[1], nmark), rep(alleles[2], nmark))
}
markerIntercross <- function(nmark, alleles) {
    cbind(rep(allele[1], nmark),
          rep(allele[2], nmark))[,sample(1:2, 2)]
}
markerBackcross <- function(nmark, alleles) {
    cbind(rep(allele, nmark),
          rep(alleles[sample(1:2, 1)], nmark))
}
markerRandom <- function(nmark, alleles){
    matrix(sample(alleles, size = 2*nmark, replace = TRUE), ncol = 2)
}

##' @title Generate distances for a genetic sequence
##' @param nchrom integer number of chromosomes in the sequence
##' @param nmark vector of integers of length nchrom giving the number
##' of markers on each chromosome
##' @param distFun function which accepts an integer and returns a
##' vector of floats of length on less than the integer
##' @author Chris Salahub
distGenesis <- function(nchrom, nmark, distFun, ...) {
    lapply(1:nchrom, function(ind) distFun(nmark[ind]-1, ...))
}

##' @title Helpers for simulating genetic distances


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
