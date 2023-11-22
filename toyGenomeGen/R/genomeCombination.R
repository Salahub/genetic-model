### Functions to Combine and Cross Genome Objects ####################

##' Maps
##' @title Genetic mapping functions
##' @description These functions support the conversion of vectors
##' between map distance and probability of recombination.
##' @details These functions convert map distances (assumed to be in
##' centiMorgans) to probabilities of recombination or the inverse.
##' Only the Haldane and Kosambi distances are defined by default.
##' @param d numeric vector of additive map distances
##' @param pr numeric vector of probabilities
##' @return Returns a numeric vector of recombination probabilities
##' based on the corresponding map distance units.
##' @author Chris Salahub
##' @describeIn maps Haldane's map distance
mapHaldane <- function(d) {
    (1 - exp(-d/50))/2
}
##' @describeIn maps The inverse of Haldane's map distance
invHaldane <- function(pr) {
    -50*log(1 - 2*pr)
}
##' @describeIn maps Kosambi's map distance
mapKosambi <- function(d) {
    tanh(d/50)/2 # (exp(d/25) - 1)/(2*(exp(d/25) + 1))
}
##' @describeIn maps The inverse of Kosambi's map distance
invKosambi <- function(pr) {
    atanh(2*pr)*50 # 25*log((1 + 2*pr)/(1 - 2*pr))
}

##' @title Functions to simulate crossing over
##' @details These functions are core to the action of `meiose` by
##' accepting a numeric vector of probabilities of recombination and
##' a numeric vector of locations and returning the indices in the
##' location vector where cross over events occur. By default, only
##' a simple independent cross over function is included, but a
##' chi-squared based function is shown in the `customFuns` demo.
##' @param probs vector of probabilities of cross overs between
##' markers
##' @param locs vector of locations of markers on a chromosome
##' @return A vecotr giving the indices of cross over events.
##' @author Chris Salahub
crossIndep <- function(probs, locs) {
    breaks <- runif(length(probs))
    which(breaks < probs)
}

##' @title Simulating meiosis
##' @details A helper function for `sex`, `meiose` dictates the
##' recombination that occurs within a `genome` and is meant to
##' mimic the action of meiosis in the creation of gametes. By
##' changing `crossFun` to reflect a given model of crossing over,
##' effectively any internal genetic recombination model can be
##' mimicked by this function.
##' @param genome a `genome` to be recombined
##' @param probs numeric vector of probabilities of recombination
##' which override distance calculations
##' @param crossFun function which accepts a vector of probabilities
##' and a vector of marker locations and returns the indices where
##' crossovers occur
##' @param map function which accepts distances between locations in
##' `genome` and returns the probabilities of recombination for each
##' distance
##' @return A list of encoding sub-matrices split by chromosome with
##' certain rows of the encoding recombined as dictated by `probs`,
##' `crossFun`, and `map`.
##' @author Chris Salahub
meiose <- function(genome, probs = NULL, crossFun = crossIndep,
                   map = mapHaldane) {
    encodings <- split(genome$encoding, genome$chromosome)
    encodings <- lapply(encodings, # fix structure
                        function(vec) matrix(vec, ncol = 2))
    dists <- lapply(genome$location, diff)
    ## write a helper to drift a single chromosome
    chromDrift <- function(copies, probs, locs) {
        copy1 <- copies[,1]
        copy2 <- copies[,2]
        crossovers <- crossFun(probs, locs)
        for (ii in seq_along(crossovers)) {
            crossPos <- crossovers[ii]
            inter <- copy2[1:crossPos]
            copy2[1:crossPos] <- copy1[1:crossPos]
            copy1[1:crossPos] <- inter
        }
        cbind(copy1, copy2)
    }
    ## get probabilities
    if (is.null(probs)) probs <- lapply(dists, map)
    ## take the above and apply it across the genome
    drifted <- Map(chromDrift, encodings, probs, genome$location)
    drifted
}

##' @title Crossing two genomes
##' @details As its name suggests, `sex` simulates sexual reproduction
##' between two `genome`s. It first calls `meiose` on the objects
##' independently, passing in the optional arguments to `meiose`.
##' Next, independent assortment randomly assigns one copy from each
##' of the recombined gametes from the meiosis calls into a new
##' `genome` and returns it.
##' @param genome1 a `genome`
##' @param genome2 a `genome` with markers observed at the same
##' locations as `genome1`
##' @param probs1 vector of probabilities to be passed to `meiose`
##' for `genome1`
##' @param probs2 vector of probabilities to be passed to `meiose`
##' for `genome2`
##' @param crossFun function dictating crossing over in `meiosis` for
##' both genomes
##' @param map function which converts map distances to probabilities
##' @return A `genome` with the same marker locations as `genome1`
##' and `genome2` with encodings recombined according to `crossFun`.
##' @author Chris Salahub
sex <- function(genome1, genome2, probs1 = NULL, probs2 = NULL,
                map = mapHaldane, crossFun = crossIndep) {
    ## perform a distance check
    if (!identical(genome1$location, genome2$location) |
        !identical(genome1$chromosome, genome2$chromosome)) {
        stop("Markers don't match")
    }
    ## meiose alleles
    gamete1 <- meiose(genome1, probs1, crossFun = crossFun, map = map)
    gamete2 <- meiose(genome2, probs2, crossFun = crossFun, map = map)
    ## pick from the copies for each genome
    chosenCopies <- replicate(length(gamete1),
                              sample(c(1,2), size = 2, replace = TRUE),
                              simplify = FALSE)
    ## select and reassort
    offspringEnc <- Map(function(g1, g2, cps) cbind(g1[,cps[1]],
                                                    g2[,cps[2]]),
                        gamete1, gamete2, chosenCopies)
    ## return the offspring
    offspring <- list(encoding = do.call(rbind, offspringEnc),
                      alleles = genome1$alleles,
                      chromosome = genome1$chromosome,
                      location = genome1$location)
    class(offspring) <- "genome"
    offspring
}
