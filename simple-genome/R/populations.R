### Functions to Compute Correlation for a Population of Genomes #####

##' @title Create a population object
##' @param genomes list of genomes measured at the same locations
##' @return a population object with the same slots as a genome object
##' but where encoding is a list of encodings for all entries
##' @author Chris Salahub
asPopulation <- function(genomes) {
    ## perform structural check that all measure the same spots
    consistent <- lapply(genomes[2:length(genomes)],
                         function(g) {
                             identical(genomes[[1]]$location,
                                       g$location) &
                                 identical(genomes[[1]]$alleles,
                                           g$alleles)
                             }
    if (!all(consistent)) {
        stop("Not all genomes measure at the same locations")
    }
    ## extract encodings
    encodings <- lapply(genomes, function(g) g$encoding)
    ## make object
    pop <- list(encodings = encodings, chromosome = genomes[[1]]$chr,
                alleles = genomes[[1]]$alleles,
                location = genomes[[1]]$location)
    class(pop) <- "population"
    pop
}

##' @title Functions for scoring genomes
##' @param genome genome object to have its encoding summarized into
##' a vector of scores
##' @return numeric vector of scores
##' @author Chris Salahub
scoreAdditive <- function(genome) {
    genome$encoding[,1] + genome$encoding[,2]
}
scoreDominance <- function(genome) {
    apply(genome$encoding, 1, max)
}
scoreHomozygous <- function(genome) {
    genome$encoding[,1] == genome$encoding[,2]
}

##' @title Observed correlation between markers for a population
##' @param population a list of genome objects measured at the same
##' nmark genetic locations
##' @param scoring function to summarize each genome into a vector
##' @return a correlation matrix of nmark x nmark
##' @author Chris Salahub
popCorrelation <- function(population, scoring = scoreAdditive) {
    popScores <- sapply(population, scoring)
    cor(t(popScores))
}

##' @title Theoretical correlation for a genome
##' @param genome genome object with specified distances
##' @param map function to convert distances to probabilities of
##' recombination
##' @param setting population setting (one of backcross, intercross,
##' or halfback) that determines the constant multiplying the
##' probabilities to be correlations
##' @return the correlation derived from the distances of genome
##' @author Chris Salahub
theoryCorrelation <- function(genome, map = mapHaldane,
                              setting = "intercross") {
    gamma <- switch(setting, # setting coefficient
                    intercross = 1,
                    backcross = 1,
                    halfback = 1/sqrt(2),
                    0)
    dists <- diff(genome$location) # only distances are used
    chrinds <- c(0, cumsum(table(genome$chromosome)))
    pos <- lapply(dists, function(el) c(0, cumsum(el)))
    M <- chrinds[length(dists)+1] # dimension of matrix
    mat <- matrix(0, nrow = M, ncol = M)
    for (ii in 1:(length(chrinds)-1)) { # fill each block
        mat[(chrinds[ii]+1):chrinds[ii+1],
        (chrinds[ii]+1):chrinds[ii+1]] <-
            (1 - 2*map(abs(outer(pos[[ii]], pos[[ii]], FUN = `-`))))
    }
    gamma*mat # multiply by constant and return
}

##' @title Visualizing genomic correlation
##' @param corrs correlation matrix
##' @param ... optional arguments to pass to image
##' @return nothing, but plot the correlations in a heat map
##' @author Chris Salahub
## write a small wrapper for the image function to place the diagonal
corrImg <- function(corrs, ...) {
    newcorr <- t(apply(corrs, 1, rev))
    image(newcorr, ...)
}

##' @title Adding chromosome boundaries and labels to the plot
##' @param genome genome object representing the correlation plot
##' @param borders boolean, should borders around each chromosome
##' block be added?
##' @return nothing, but add lines separating chromosome blocks to a
##' plotted correlation heat map
##' @author Chris Salahub
addChromosomeLines <- function(genome, borders = FALSE, ...) {
    runs <- rle(genome$chromosome) # run lengths by chromosome
    nmarks <- runs$lengths # counts of each chromosome
    wid <- 1/(sum(nmarks)-1) # relative width unit
    chrTab <- nmarks[ord]*wid # chromosome widths
    postns <- cumsum(c(-wid/2, chrTab))
    for (ii in 1:length(chrTab)) { # add lines, labels
        abline(v = c(postns[ii], postns[ii+1]),
               h = c(1 - postns[ii], 1 - postns[ii+1]),
               col = "gray70")
        mtext(runs$values[ii], side = 3,
              at = postns[ii]/2 + postns[ii+1]/2)
        mtext(runs$values[ii], side = 2,
              at = 1 - postns[ii]/2 - postns[ii+1]/2,
              las = 1)
    }
    if (borders) { # add borders around chromosomes
        for (ii in 1:length(chrTab)) {
            rect(postns[ii], 1 - postns[ii], postns[ii+1],
                 1 - postns[ii+1], ...)
        }
    }
}

