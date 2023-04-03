library(grid)

## MAIN FUNCTIONS ####################################################
## generate a genome based on specifications
abiogenesis <- function(nchrom, nmark, setting = "pure", dists = NULL,
                        allele = 1,
                        distFun = function(n) runif(n, 50, 100)) {
    ## check marker specs
    if (length(nmark) != nchrom) nmark <- rep(nmark, length.out = nchrom)
    ## generate genome object
    genome <- list(alleles = genomeGenesis(nchrom, nmark,
                                         setting = setting,
                                         allele = allele),
                   dists = if (is.null(dists)) {
                               distGenesis(nchrom, nmark,
                                           distFun = distFun)
                           } else dists)
    class(genome) <- "genome"
    genome
}

## write a function to drift a given genome (meiosis event)
meiose <- function(genome, probs = NULL, map = mapHaldane) {
    alleles <- genome$alleles
    dists <- genome$dists
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
    newGenome <- list(alleles = drifted, dists = dists)
    class(newGenome) <- "genome"
    newGenome
}

## now a function that crosses two given genomes (sex)
sex <- function(genome1, genome2, map = mapHaldane) {
    ## perform a distance check
    if (!identical(genome1$dists, genome2$dists)) {
        stop("Allele distances don't match")
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

##' a "scoring" function, consider the possible variants
##' - additive: the count of A's at each site
##' - dominant: whether the gene is heterozygous or homozygous
##' TODO: read this better... still unclear
scoreGenome <- function(genome, score = "additive") {
    scoreFun <- switch(score,
                       additive = `+`,
                       dominant = `==`,
                       stop(paste(score, "not implemented")))
    lapply(genome$alleles,
           function(chrom) scoreFun(chrom[,1], chrom[,2]))
}

## write a simple plot function
plot.genome <- function(genome) {
    xpos <- ppoints(length(genome$alleles))
    delta <- if (length(xpos) == 1) 0.1 else diff(xpos)[1]/5
    yscales <- lapply(genome$dists, function(el) c(0, cumsum(el)))
    yabs <- max(unlist(yscales))
    ypos <- lapply(yscales,
                   function(el) 0.95 - (el/yabs)*0.9)
    grid.newpage()
    for (ii in seq_along(genome$alleles)) {
        chrom <- genome$alleles[[ii]]
        grid.text(chrom,
                  x = xpos[ii] + rep(c(-1,1)*delta, each = nrow(chrom)),
                  y = rep(ypos[[ii]], times = 2),
                  gp = gpar(col = c("steelblue", "firebrick")[chrom+1]))
    }
}

## calculate correlation for a series of genome scores
popCorrelation <- function(population, scoring = scoreGenome) {
    ## check if genomes or scores have been provided
    if (all(sapply(population, class) == "genome")) { # score genomes
        population <- lapply(population, scoring)
    }
    ## unlist the scores, generate correlation matrix
    fullScore <- t(sapply(population, unlist))
    cor(fullScore)
}

## write a small wrapper for the image function to place the diagonal
corrImg <- function(corrs, ...) {
    newcorr <- t(apply(corrs, 1, rev))
    image(newcorr, ...)
}

## add chromosome boundaries to the plot
addChromosomes <- function(marks, ord, ...) {
    wid <- 1/(nrow(marks)-1)
    chrTab <- table(marks$chr)[ord]*wid
    postns <- cumsum(c(-wid/2, chrTab))
    for (ii in 1:length(chrTab)) {
        abline(v = c(postns[ii], postns[ii+1]),
               h = c(1 - postns[ii], 1 - postns[ii+1]),
               col = "gray70")
        mtext(names(chrTab)[ii], side = 3,
              at = postns[ii]/2 + postns[ii+1]/2)
        mtext(names(chrTab)[ii], side = 2,
              at = 1 - postns[ii]/2 - postns[ii+1]/2,
              las = 1)
    }
    #for (ii in 1:length(chrTab)) {
    #    rect(postns[ii], 1 - postns[ii], postns[ii+1],
    #         1 - postns[ii+1], ...)
    #}
}

## add a theoretical correlation calculation function
theoryCorrelation <- function(d, map = mapHaldane,
                              setting = "intercross") {
    gamma <- switch(setting, # setting coefficient
                    intercross = 1,
                    backcross = 1,
                    halfback = 1/sqrt(2),
                    0)
    gamma*(1 - 2*map(d))
}

## apply it to generate theoretical correlation matrices
theoryCor <- function(dists, map = mapHaldane,
                      setting = "intercross") {
    inds <- c(0, cumsum(sapply(dists,
                               function(el) length(el) + 1)))
    pos <- lapply(dists, function(el) c(0, cumsum(el)))
    M <- inds[length(dists)+1]
    mat <- matrix(0, nrow = M, ncol = M)
    for (ii in 1:(length(inds)-1)) {
        mat[(inds[ii]+1):inds[ii+1],
            (inds[ii]+1):inds[ii+1]] <-
               theoryCorrelation(abs(outer(pos[[ii]],
                                           pos[[ii]],
                                           FUN = `-`)),
                                 map,
                                 setting)
    }
    mat
}
