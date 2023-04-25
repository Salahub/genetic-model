### Functions to Compute Correlation for a Population of Genomes #####

##' @title Create a population object
##' @details Provided a list of `genome`s, this function strips
##' redundant information from the genomes to convert them to a
##' `population`, which has the same slots as `genome` except that
##' `encoding` is now `encodings` to reflect that it is a list of
##' encoding matrices and `marker` is added to remove redundant row
##' names.
##' @param genomes list of `genome`s with the same locations
##' @return A population object.
##' @author Chris Salahub
asPopulation <- function(genomes) {
    ## perform structural check that all measure the same spots
    consistent <- lapply(genomes[2:length(genomes)],
                         function(g) {
                             identical(genomes[[1]]$location,
                                       g$location) &
                                 identical(genomes[[1]]$alleles,
                                           g$alleles)
                             })
    if (!all(unlist(consistent))) {
        stop("Not all genomes measure at the same locations")
    }
    ## get marker names
    marker <- rownames(genomes[[1]]$encoding)
    ## extract encodings without names
    encodings <- lapply(genomes,
                        function(g) {
                            enc <- g$encoding
                            rownames(enc) <- NULL
                            enc
                        })
    ## make object
    pop <- list(encodings = encodings, chromosome = genomes[[1]]$chr,
                alleles = genomes[[1]]$alleles, marker = marker,
                location = genomes[[1]]$location)
    class(pop) <- "population"
    pop
}

##' @title Subsetting populations by markers
##' @details Like the `genome` analog, this subsets markers, drops
##' unused chromosomes, and returns a new population object. Unlike
##' the single `genome` case, this includes a second argument to
##' drop individual genomes from the population.
##' @param population `population` to be subset by markers
##' @param markInd indices of markers to keep
##' @param genInd indices of individuals to keep
##' @return A population object including only the encoding matrices
##' of `genInd` subset to the markers of `markInd`.
##' @author Chris Salahub
subsetPopulation <- function(population,
                             markInd = (1:length(population$chromosome)),
                             genInd = (1:length(population$encoding))) {
    newchr <- droplevels(population$chromosome[markInd])# subset chroms
    newloc <- split(unlist(population$location,
                           use.names = FALSE)[markInd], newchr)
    structure(list(encodings = lapply(population$encodings[genInd],
                                      function(enc) enc[markInd,]),
                   chromosome = newchr,
                   marker = population$marker[markInd],
                   alleles = population$alleles[markInd],
                   location = newloc),
              class = "population")
}

##' @title Selecting a genome from a population
##' @details This function takes an individual encoding from a
##' population and converts it back to a `genome` object for
##' plotting and manipulation.
##' @param population `population` to be subset by markers
##' @param ind index of individual encoding to extract
##' @return A `genome` with the encoding of the individual at
##' `ind` in `population`.
##' @author Chris Salahub
selectGenome <- function(population, ind = 1) {
    enc <- population$encodings[[ind]]
    rownames(enc) <- population$marker
    structure(list(encoding = enc,
                   chromosome = population$chromosome,
                   alleles = population$alleles,
                   location = population$location),
              class = "genome")
}

##' @title Remove encodings from a `population` based on a rule
##' @details A simple wrapper function which applies a function that
##' accepts an encoding matrix and returns a logical to every encoding
##' in a `population` and keeps only those which return `TRUE`
##' @param population `population` with encodings to be filtered
##' @param rule function to apply to each encoding that returns a
##' boolean, defaults to checking if any values are NA
##' @return A population with only those encodings for which rule
##' returns `TRUE`
##' @author Chris Salahub
filterPopulation <- function(population,
                             rule = function(mat) any(is.na(mat))) {
    trues <- sapply(population$encodings, rule)
    structure(list(encodings = population$encodings[trues],
                   chromosome = population$chromosome,
                   alleles = population$alleles,
                   marker = population$marker,
                   location = population$location),
              class = "population")
}

##' Scoring
##' @title Summarizing an encoding ("scoring" the genome)
##' @description These functions convert two-column matrices of
##' encodings into summary vectors based on an operation applied to
##' each row.
##' @details These functions accept encoding matrices and return
##' numeric vectors that summarize the marker measurements at each
##' row.
##' @param encoding matrix with two columns giving marker encodings
##' @return A numeric vector of scores with one numeric for each row
##' of `encoding`.
##' @author Chris Salahub
##' @describeIn scoring The additive score
scoreAdditive <- function(encoding) {
    encoding[,1] + encoding[,2]
}
##' @describeIn scoring The dominance score
scoreDominance <- function(encoding) {
    apply(encoding, 1, max)
}
##' @describeIn scoring The homozygous score
scoreHomozygous <- function(encoding) {
    as.numeric(encoding[,1] == encoding[,2])
}

##' @title Compute observed correlation for a `population`
##' @details A wrapper which computes the correlation for a
##' `population` given a scoring function to apply to each
##' encoding.
##' @param population `population` object
##' @param scoring function which accepts an encoding matrix and
##' returns a numeric vector with one numeric for each row of the
##' encoding matrix
##' @param ... additional arguments to pass to `cor`
##' @return A correlation matrix between markers in `population`.
##' @author Chris Salahub
popCorrelation <- function(population, scoring = scoreAdditive, ...) {
    popScores <- sapply(population$encodings, scoring)
    cor(t(popScores), ...)
}

##' @title Compute the theoretical correlation in a genome
##' @details Using a given map distance function and population cross
##' setting, returns the theoretical correlation, proportional to
##' $$1/2 - p_r(d)/2,$$
##' where $p_r(d)$ is the probability of recombination for markers
##' with distance $d$ between them.
##' @param genome `genome` or `population` with slots giving
##' chromosomes and locations of markers
##' @param map function which accepts a vector of numerics giving
##' distances between locations and returns a vector of probabilities
##' of recombination for the corresponding distances
##' @param setting population setting (one of backcross, intercross,
##' or halfback) that determines the constant for computing
##' correlation
##' @return The correlation derived from the distances between markers
##' in `genome` of `population`.
##' @author Chris Salahub
theoryCorrelation <- function(genome, map = mapHaldane,
                              setting = "intercross") {
    gamma <- switch(setting, # setting coefficient
                    intercross = 1,
                    backcross = 1,
                    halfback = 1/sqrt(2),
                    0)
    chrinds <- c(0, cumsum(table(genome$chromosome)))
    pos <- genome$location # shorthand
    pr <- lapply(pos, # matrix of recombination probabilities
                 function(ls) 1-2*map(abs(outer(ls, ls, FUN = `-`))))
    M <- chrinds[length(pos)+1] # dimension of matrix
    mat <- matrix(0, nrow = M, ncol = M)
    for (ii in seq_along(pr)) { # fill each block
        if (length(pr[[ii]]) != 0) { # safe replacement for zero len
            mat[(chrinds[ii]+1):chrinds[ii+1],
            (chrinds[ii]+1):chrinds[ii+1]] <- pr[[ii]]
        }
    }
    gamma*mat # multiply by constant and return
}

##' @title Visualizing genomic correlation
##' @details A simple wrapper for `image` that places the correlation
##' in a more intuitive form where the 1,1 entry is in the top left.
##' @param corrs a numeric correlation matrix
##' @param ... optional arguments to pass to `image`
##' @return Returns nothing, but plots the correlations in a heat map
##' in the active device.
##' @author Chris Salahub
corrImg <- function(corrs, ...) {
    newcorr <- t(apply(corrs, 1, rev))
    image(newcorr, ...)
}

##' @title Add chromosome boundaries and labels to a plot
##' @details Meant to be used after calling `corrImg`, this function
##' adds vertical and horizontal lines to a plot to visually separate
##' the chromosomes of a correlation matrix and labels the bordered
##' regions accordingly.
##' @param genome `genome` or `population` corresponding to the
##' correlation matrix plotted on the active device
##' @param borders logical, should extra borders around each
##' chromosome block on the diagonal be added?
##' @param lncol character giving the colour of lines
##' @param ... additional arguments to pass to rect if borders = TRUE
##' @return Returns nothing, but adds lines separating chromosomes
##' to a plotted correlation heat map
##' @author Chris Salahub
addChromosomeLines <- function(genome, borders = FALSE,
                               lncol = "gray70", ...) {
    tab <- table(genome$chromosome)
    wid <- 1/(sum(tab)-1) # relative width unit
    chrTab <- tab*wid # chromosome widths
    postns <- cumsum(c(-wid/2, chrTab))
    for (ii in 1:length(chrTab)) { # add lines, labels
        abline(v = c(postns[ii], postns[ii+1]),
               h = c(1 - postns[ii], 1 - postns[ii+1]),
               col = lncol)
        mtext(names(tab)[ii], side = 3,
              at = postns[ii]/2 + postns[ii+1]/2)
        mtext(names(tab)[ii], side = 2,
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

##' @title S3 methods for `population`
##' @description Print method
##' @param population `population` object
##' @return Returns nothing, but prints a quick summary of the
##' populations to console
##' @author Chris Salahub
print.population <- function(population) {
    ne <- length(population$encoding)
    nm <- length(population$chromosome)
    nc <- length(population$location)
    dtab <- table(population$chromosome)
    mis <- 100*round(sum(is.na(unlist(population$encoding)))/
                     length(unlist(population$encoding)), 2)
    if (mis > 0) {
        cat("A population of", ne, "genomes encoded at", sum(dtab),
            "markers across", nc, "chromosomes, distributed:\n", dtab,
            "\nRoughly", mis, "% of the data is missing.\n")
    } else {
        cat("A population of", ne, "genomes encoded at", sum(dtab),
            "markers across", nc, "chromosomes, distributed:\n", dtab,
            "\n")

    }
}
