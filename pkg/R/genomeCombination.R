### Functions to Combine and Cross Genome Objects ####################

##' @title Computing recombination probability from map distances
##' @param d numeric vector of additive map distances, for Kosambi and
##' Haldane's map functions in centiMorgans or their equivalent
##' @return numeric vector of recombination probabilities based on the
##' corresponding map distance units
##' @author Chris Salahub
mapHaldane <- function(d) {
    1 - exp(-d/50))/2
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
#mapRenewal <- function(d, density = function(x) dchisq(x, df = 4)/4,
#                       ...) {
#    c0 <- integrate(function(y) integrate(density, y, Inf), d, Inf)
#    (1 - c0)/2
#}

##' @title Simulating meiosis
##' @param genome a genome object which will be meiosed
##' @param probs (optional) custom probabilities of recombination
##' which override distance calculations
##' @param map function which accepts distances stored in the genome
##' object and returns probabilities of recombination
##' @return a new genome object with certain rows of the encoding
##' recombined
##' @author Chris Salahub
meiose <- function(genome, probs = NULL, map = mapHaldane) {
    encodings <- split(genome$alleles, genome$chromosome)
    dists <- genome$distances
    ## write a helper to drift a single chromosome
    chromDrift <- function(copies, probs) {
        copy1 <- copies[,1]
        copy2 <- copies[,2]
        breaks <- runif(length(probs))
        crossovers <- which(breaks < probs)
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
    drifted <- Map(chromDrift, alleles, probs)
    newGenome <- list(encoding = do.call(rbind, drifted),
                      distances = genome$distances,
                      alleles = alleles,
                      chromosome = genome$chromosome)
    class(newGenome) <- "genome"
    newGenome
}

## now a function that crosses two given genomes (sex)
sex <- function(genome1, genome2, map = mapHaldane) {
    ## perform a distance check
    if (!identical(genome1$dists, genome2$dists)) {
        stop("Markers don't match")
    }
    ## meiose alleles
    genome1 <- meiose(genome1, map = map)
    genome2 <- meiose(genome2, map = map)
    ## get alleles
    alleles1 <- genome1$alleles
    alleles2 <- genome2$alleles
    ## pick from the copies for each genome
    chosenCopies <- replicate(length(alleles1),
                              sample(c(1,2), size = 2, replace = TRUE),
                              simplify = FALSE)
    ## select and reassort
    offspring <- Map(function(g1, g2, cps) cbind(g1[,cps[1]],
                                                 g2[,cps[2]]),
                     alleles1, alleles2, chosenCopies)
    ## return the offspring
    offspring <- list(alleles = offspring, dists = genome1$dists)
    class(offspring) <- "genome"
    offspring
}
