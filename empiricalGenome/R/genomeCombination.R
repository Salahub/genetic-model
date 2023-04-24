### Functions to Combine and Cross Genome Objects ####################

##' @title Computing recombination probability from map distances
##' @param d numeric vector of additive map distances, for Kosambi and
##' Haldane's map functions in centiMorgans or their equivalent
##' @return numeric vector of recombination probabilities based on the
##' corresponding map distance units
##' @author Chris Salahub
mapHaldane <- function(d) {
    (1 - exp(-d/50))/2
}
invHaldane <- function(pr) {
    -50*log(1 - 2*pr)
}
mapKosambi <- function(d) {
    tanh(d/50)/2 # (exp(d/25) - 1)/(2*(exp(d/25) + 1))
}
invKosambi <- function(pr) {
    atanh(2*pr)*50 # 25*log((1 + 2*pr)/(1 - 2*pr))
}

##' @title Functions to simulate crossing over
##' @param probs vector of probabilities of cross overs between
##' markers
##' @param locs vector of locations of markers on a chromosome
##' @return indices of cross over events
##' @author Chris Salahub
crossIndep <- function(probs, locs) {
    breaks <- runif(length(probs))
    which(breaks < probs)
}

##' @title Simulating meiosis under non-interference
##' @param genome a genome object which will be meiosed
##' @param probs (optional) custom probabilities of recombination
##' which override distance calculations
##' @param crossFun (optional) function which accepts a vector of
##' probabilities and marker positions and returns the indices where
##' crossovers occurred
##' @param map function which accepts distances stored in the genome
##' object and returns probabilities of recombination
##' @return a list of encodings split by chromosome with certain rows
##' of the encoding recombined
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
##' @param genome1 a genome object
##' @param genome2 a genome object with markers observed at the same
##' locations as on genome 1
##' @param probs1 (optional) probabilities to be passed into meiose
##' for genome1
##' @param probs2 (optional) probabilities to be passed into meiose
##' for genome2
##' @param crossFun (optional) function to dictate crossing over for
##' both genomes
##' @param map function to convert map distances to probabilities
##' @return a genome object with the same marker locations as genome1
##' and genome2 with its encodings recombined
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
