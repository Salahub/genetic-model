### Generator Functions for Genome Objects ###########################

##' @title Generate allele settings
##' @param nmark named vector of the number of measured markers by
##' chromosome
##' @return a list the same length as the sum of nmark with each entry
##' c("A", "a")
##' @author Chris Salahub
generateBiAllelic <- function(nmark) {
    rep(list(c("A","a")), sum(nmark))
}

##' @title Generate a genetic sequence
##' @param nmark a named vector of number of measured markers by
##' chromosome
##' @param markerFuns a list of functions which each generate measured
##' markers for a chromosome based on the alleles, assuming the
##' ordering of alleles is meaningful, recycled if length does not
##' match nmark
##' @return a list containing a 2-column matrix of marker encodings,
##' a list of possible alleles, and a vector of chromosome membership
##' @author Chris Salahub
generateEncoding <- function(nmark,
                             markerFuns = markerPureDom) {
    if (!is.list(markerFuns)) {
        enc <- mapply(markerFuns, n = nmark, SIMPLIFY = FALSE)
    } else {
        enc <- mapply(do.call, markerFuns, lapply(nmark, as.list),
                      SIMPLIFY = FALSE)
    }
    do.call(rbind, enc)
}

##' @title Helpers to simulate marker encodings
##' @param n integer number of markers
##' @return a 2-column numeric matrix with entries given by the
##' corresponding setting
##' @author Chris Salahub
markerPureDom <- function(n) {
    cbind(mv = rep(1, n), pv = rep(1, n))
}
markerPureRec <- function(n) {
    cbind(mv = rep(0, n), pv = rep(0, n))
}
markerHybrid <- function(n) {
    cbind(mv = rep(1, n), pv = rep(0, n))[,sample(1:2, 2)]
}
markerRandom <- function(n){
    cbind(mv = sample(c(1,0), size = n, replace = TRUE),
          pv = sample(c(1,0), size = n, replace = TRUE))
}

##' @title Generate locations for a genetic sequence
##' @param nmark vector of integers of giving the number of markers on
##' each chromosome for which distances are generated
##' @param locFuns a list of functions which each accept an integer
##' and return a numeric vector with length equal to that integer,
##' recycled if not of the same length as nmark
##' @param ... (optional) a list of additional arguments to each
##' function in locFuns, recycled if not the same length
##' @return a list of numeric vectors giving the distances between
##' markers on each chromosome
##' @author Chris Salahub
generateLocations <- function(nmark, locFuns, ...) {
    if (!is.list(locFuns)) { # simple handling of single function
        lapply(nmark, locFuns, ...)
    } else { # more complex
        nmark <- lapply(nmark, as.list) # for argument creation
        args <- mapply(c, nmark, ..., SIMPLIFY = FALSE)
        lFs <- rep(locFuns, length.out = length(nmark))
        mapply(do.call, lFs, args, SIMPLIFY = FALSE)
    }
}

##' @title Helpers for simulating genetic distances
##' @param nmark integer number of markers between which distances
##' are generated
##' @param delta constant distance to repeat
##' @param min smallest numeric location, must be zero or greater
##' @param max largest numeric location
##' @return a vector of length nmark - 1 drawn from dists
##' @author Chris Salahub
locationUniform <- function(nmark, min = 0, max = 100) {
    runif(nmark, min, max)
}
locationRegular <- function(nmark, delta = 50) {
    delta*(1:nmark)
}
