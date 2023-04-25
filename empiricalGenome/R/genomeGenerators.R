### Generator Functions for Genome Objects ###########################

##' @title Generate allele settings
##' @details This function, called when alleles are not specified in
##' `makeGenome`, creates a list by repeating `c("A", "a")`.
##' @param nmark named integer vector of the count of measured markers
##' by chromosome
##' @return A list the same length as the sum of nmark with each entry
##' `c("A", "a")`
##' @author Chris Salahub
generateBiAllelic <- function(nmark) {
    rep(list(c("A","a")), sum(nmark))
}

##' @title Generate a genetic encoding
##' @details A wrapper function for marker-generating functions,
##' this takes a list of marker-generating functions and an integer
##' vector of counts of the same length and generates encodings by
##' calling each marker function on the corresponding count. The
##' marker functions are recycled by default.
##' @param nmark named integer vector of measured markers by
##' chromosome
##' @param markerFuns list of functions which each generate encodings
##' for a chromosome given `n`, the number of markers on the
##' chromosome
##' @param ... list of additional arguments to each function in
##' `markerFuns`, recycled if not the same length
##' @return A two-column matrix of encodings for the entire genome.
##' @author Chris Salahub
generateEncoding <- function(nmark, markerFuns = markerPureDom,
                             ...) {
    if (!is.list(markerFuns)) { # single function case
        enc <- lapply(nmark, markerFuns, ...)
    } else { # more complicated
        nmakr <- lapply(nmark, as.list) # argument list
        args <- mapply(c, nmark, ..., SIMPLIFY = FALSE) # add extra
        mFs <- rep(markerFuns, length.out = length(nmark))
        enc <- mapply(do.call, mFs, args, SIMPLIFY = FALSE)
    }
    do.call(rbind, enc)
}

##' Encodings
##' @title Simulate genetic encodings
##' @description These functions simulate the given number of marker
##' encodings under different settings.
##' @details These functions all accept `n`, the number of markers,
##' and generate a two column matrix of encodings. The included
##' functions all correspond to known crosses or settings useful to
##' experiment with recombination.
##' @param n integer number of markers
##' @return A 2-column numeric matrix of encodings.
##' @author Chris Salahub
##' @describeIn encodings Simulating a pure dominant encoding
markerPureDom <- function(n) {
    cbind(mv = rep(1, n), pv = rep(1, n))
}
##' @describeIn encodings Simulating a pure recessive encoding
markerPureRec <- function(n) {
    cbind(mv = rep(0, n), pv = rep(0, n))
}
##' @describeIn encodings Simulating a hybrid encoding
markerHybrid <- function(n) {
    cbind(mv = rep(1, n), pv = rep(0, n))[,sample(1:2, 2)]
}
##' @describeIn encodings Simulating a random encoding
markerRandom <- function(n){
    cbind(mv = sample(c(1,0), size = n, replace = TRUE),
          pv = sample(c(1,0), size = n, replace = TRUE))
}

##' @title Generate locations for markers
##' @details A wrapper function for location-generating functions that
##' accepts a list of functions and an integer vector of counts of the
##' same length and returns a list of locations resulting from the
##' application of the functions to the corresponding counts.
##' @param nmark vector of integers giving the number of locations by
##' chromosome
##' @param locFuns list of functions which each accept an integer
##' and return a numeric vector with length equal to that integer,
##' recycled if not of the same length as `nmark`
##' @param ... list of additional arguments to each function in
##' `locFuns`, recycled if not the same length
##' @return A list of numeric vectors giving the locations of markers
##' on each chromosome.
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

##' Locations
##' @title Simulate marker locations
##' @description These functions generate locations given a number of
##' markers based on different settings.
##' @details These functions accept an integer `n` and potentially
##' other arguments and return a vector of length `n` giving locations
##' on a chromosome. The default included functions are
##' `locationUniform`, which randomly samples uniform locations
##' between `min` and `max` cM, and `locationRegular`, which assumes
##' a constant distance of `delta` cM between each marker.
##' @param nmark integer number of markers
##' @param delta constant distance to repeat
##' @param min smallest numeric location, must be zero or greater
##' @param max largest numeric location
##' @return A vector of length `nmark` giving marker locations.
##' @author Chris Salahub
##' @describeIn locations Simulating uniform locations
locationUniform <- function(nmark, min = 0, max = 100) {
    runif(nmark, min, max)
}
##' @describeIn locations Simulating evenly-spaced locations
locationRegular <- function(nmark, delta = 50) {
    delta*(1:nmark)
}
