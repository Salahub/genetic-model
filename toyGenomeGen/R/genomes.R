## MAIN FUNCTIONS ####################################################

##' @title Simulate a genome object
##' @description `simGenome` simulates a genome object based on the
##' distribution of markers in `nmark` and the marker and location
##' generating functions in `markerFuns` and `locFuns`
##' @details `alleles` must be a list of the same length as
##' `sum(nmark)`, with each element providing the possible alleles for
##' the corresponding marker. If it is not provided, the Mendelian
##' A/a biallelic annotation is assumed for every marker. If either of
##' `markerFuns` or `locFuns` provides a function rather than a list
##' of functions, this function is recycled across every chromosome.
##' Otherwise, the length of the list must match the number of
##' chromosomes to provide the different functions used to generate
##' encodings and locations on each chromosome.
##' @param nmark vector of integers giving the count of measured
##' markers by chromosome, if it is named it is assumed the names
##' correspond to the chromosomes, otherwise chromosomes are numbered
##' in order
##' @param markerFun list of functions that simulate marker
##' encodings by chromosome, recycled if it is not as long as nmark
##' @param locFuns list of functions that generate locations for
##' markers by chromosome, recycled if it is not as long as nmark
##' @param alleles vector of possible allele annotations, if missing
##' this is taken to be `c("A","a")` for every marker site
##' @return A genome object with encodings and locations generated
##' by the corresponding functions for each chromosome.
##' @examples
##' simGenome(c("1" = 10, "2" = 5, "X" = 19))
##' @author Chris Salahub
simGenome <- function(nmark, alleles, markerFuns = markerHybrid,
                      locFuns = locationRegular) {
    ## check marker names and generate a chromosome indicator
    if (is.null(names(nmark))) {
        chr <- factor(rep(1:length(nmark), times = nmark))
    } else {
        chr <- factor(rep(names(nmark), times = nmark),
                      levels = names(nmark))
    }
    ## generate alleles if necessary
    if (missing(alleles)) alleles <- generateBiAllelic(nmark)
    ## create encodings
    encoding <- generateEncoding(nmark, markerFuns)
    ## create distances
    locations <- generateLocations(nmark, locFuns)
    ## collect into genome object
    genome <- list(encoding = encoding, chromosome = chr,
                   alleles = alleles, location = locations)
    class(genome) <- "genome"
    genome
}

##' @title Making a genome from provided slots
##' @description `makeGenome` checks arguments for conformity and
##' then stores them in a `genome` object. `genome` is an S3 class
##' representing the encoding matrix of a diploid genome.
##' @details A simple constructor that checks arguments and provides
##' descriptive errors to ensure a valid genome is constructed.
##' @param location list of marker locations by chromosome
##' @param alleles list of alleles by marker location
##' @param chromosome factor giving chromosome values with levels
##' ordered as they should appear in plots
##' @param encoding numeric two-column matrix with the same number
##' of rows as the length of `chromosome`, `alleles`
##' @return A genome object with slots given by the provided
##' arguments.
##' @author Chris Salahub
makeGenome <- function(location, alleles, chromosome,
                       encoding = markerHybrid(length(chromosome))) {
    ## check arguments match
    if (!all(sapply(location, length) == table(chromosome))) {
        stop("Locations and chromosomes imply differing structure")
    }
    if (length(alleles) != length(chromosome)) {
        stop("Alleles and chromosomes imply differing structure")
    }
    if (nrow(encoding) != length(chromosome)) {
        stop("Encoding does not match structure of chromosome, location")
    }
    if (ncol(encoding) != 2) {
        stop("Encoding must have two columns")
    }
    if (!is.factor(chromosome)) {
        warning("Warning: chromosome is not an ordered factor")
        chromosome <- factor(chromosome, levels = unique(chromosome))
    }
    structure(list(encoding = encoding, chromosome = chromosome,
                   alleles = alleles, location = location),
              class = "genome")
}

##' @title Subsetting the markers of a genome object
##' @details A helper function that subsets every slot of a `genome`
##' and drops unused `chromosome` levels by a set of marker indices.
##' @param genome object to be subset
##' @param inds indices of markers to be kept
##' @return A genome object containing only the markers at `inds`,
##' with any chromosomes fully removed dropped as levels of
##' `genome$chromosome`.
##' @author Chris Salahub
subsetGenome <- function(genome, inds) {
    newchr <- droplevels(genome$chromosome[inds]) # subset chromosomes
    newloc <- split(unlist(genome$location)[inds], newchr)
    structure(list(encoding = genome$encoding[inds,],
                   chromosome = newchr,
                   alleles = genome$alleles[inds],
                   location = newloc),
              class = "genome")
}

##' @title Checking if an object is a valid genome
##' @details Performs checks of slot data types, dimensions, and
##' conformity to ensure `genome` is a valid instance of the S3 class
##' `genome`.
##' @param genome object to be checked
##' @return TRUE if the object is a genome, otherwise text outlining
##' what check failed.
##' @author Chris Salahub
checkGenome <- function(genome) {
    isgen <- TRUE # change if false
    ## check names of elements
    if (sort(names(genome)) != c("alleles", "chromosome", "location",
                                 "encoding")) {
        isgen <- FALSE
        cat("Slots namely incorrectly or missing\n")
    } else if (!is.list(genome$alleles) | # check types of elements
               !is.vector(genome$chromosome) |
               !is.list(genome$location) |
               !is.matrix(genome$encoding)) {
        isgen <- FALSE
        cat("Slots contain data of the incorrect type\n")
    } else if (length(genome$alleles) != nrow(genome$encoding)) {
        isgen <- FALSE
        cat("Number of alleles does not match encoding\n")
    } else if (!identical(table(genome$chromosome),
                          sapply(genome$location, length))) {
        isgen <- FALSE
        cat("Locations do not match number of markers by chromosome\n")
    } else if (length(genome$alleles) != length(genome$chromosome)) {
        isgen <- FALSE
        cat("Distances and alleles imply differing counts of markers\n")
    } else if (!is.factor(genome$chromosome)) {
        cat("Chromosome is not an ordered factor\n")
    }
    isgen
}

##' @title Simple encoding function
##' @details This function takes observed annotations, a vector of
##' possible annotations, and a vector of values for each, and returns
##' an encoding given by replacing annotations with values.
##' @param annotation character vector giving observed annotations
##' @param alleles character vector of possible annotations
##' @param values numeric vector of the same length as alleles giving
##' the value of each allele when encoded
##' @return A numeric vector of encodings.
##' @author Chris Salahub
encode <- function(annotation, alleles, values) {
    names(values) <- alleles
    values[annotation]
}

##' @title Convert a `data.frame` to a `genome`
##' @details Given a `data.frame` with annotations in a particular
##' structure, potentially with missing annotations denoted by a set
##' list of characters, this function generates the corresponding
##' `genome` by encoding the data and changing its format.
##' @param df `data.frame` with columns `mv`, `pv`, `chr`, and `pos`
##' @param alleles list with `length(alleles) == nrow(df)` that
##' gives the potential alleles for each row
##' @param values list of identical form to `alleles` giving numeric
##' values for the encoding of each possible allele
##' @param missing character vector giving string values to be
##' interpreted as missing
##' @param encoder function which accepts an annotation and its
##' possible values and returns an encoding
##' @return A `genome` with locations given by `pos` across
##' chromosomes given by `chr` with encodings given by the values in
##' `mv` and `pv`
##' @author Chris Salahub
asGenome <- function(df, alleles, values = c(1,0),
                     missing = c(".", "-"), encoder = encode) {
    ## safety checks
    stopifnot(all(c("mv", "pv", "chr", "pos") %in% names(df)))
    if (is.list(alleles) & (length(alleles) != nrow(df))) {
        stop("Length of list of alleles must match number of rows in df")
    }
    if (!all(sapply(alleles, length) == sapply(values, length))) {
        stop("Lengths of alleles must match lengths of values")
    }
    ## replace missing values
    df$mv[df$mv %in% missing] <- NA
    df$pv[df$pv %in% missing] <- NA
    ## missing alleles
    if (missing(alleles)) {
        warning("alleles not provided, assuming upper/lower case")
        alleles <- lapply(df$mv,
                          function(an) {
                              c(toupper(an), tolower(an))
                          })
    }
    ## alleles/values not lists: encode all at once
    if (!is.list(alleles)) {
        enc <- matrix(encode(as.matrix(df[, c("mv", "pv")]),
                             alleles, values),
                      nrow = nrow(df))
        ## repeat for final object
        alleles <- rep(list(alleles), length = nrow(df))
        values <- rep(list(values), length = nrow(df))
    } else { # encode line by line
        enc <- sapply(1:nrow(df),
                  function(ii) encoder(df[ii, c("mv", "pv")],
                                       alleles[[ii]], values[[ii]]))
    }
    rownames(enc) <- rownames(df) # keep marker names if present
    chrom <- if (is.factor(df$chr)) {
                 df$chr }
             else factor(df$chr, levels = unique(df$chr))
    ## collect into genome object
    genome <- list(encoding = enc,
                   chromosome = chrom,
                   alleles = alleles,
                   location = split(df$pos, df$chr))
    class(genome) <- "genome"
    genome
}

##' Methods
##' @title S3 methods for `genome`
##' @description The `print` and `plot` methods outlined here support
##' the quick description and exploration of a `genome`.
##' @details These methods for printing and plotting genome objects
##' are helpful to get a quick overview of the encoding contained in
##' a genome and the structure of the locations encoded.
##' @param genome instance of the genome class
##' @param chrLens numeric vector of lengths of each chromosome, if
##' missing taken as the largest location on each chromosome
##' @param scale logical: if chrLens is provided, should the locations
##' be scaled to [0,1] when plotted?
##' @param epch vector of R pch values for the potential encoding
##' combinations
##' @param elev vector of possible encodings generated by pasting
##' the two columns of the encoding matrix together, defaults to
##' "0 0", "0 1", "1 0", "1 1"
##' @param add.legend if TRUE, draw a legend in the top right corner,
##' if FALSE suppress legend
##' @return Nothing is returned, but a quick summary of the genome is
##' printed to console or a visualization of the genome is plotted
##' on the active device.
##' @author Chris Salahub
##' @describeIn methods Print method for `genome`
print.genome <- function(genome) {
    nm <- length(genome$chromosome)
    nc <- length(genome$location)
    dtab <- table(genome$chromosome)
    mis <- 100*round(sum(is.na(genome$encoding))/ # missing data?
                     length(genome$encoding), 2)
    if (mis > 0) {
        cat("A genome object encoding", nm, "markers across", nc,
            "chromosomes, distributed:\n", dtab, "\nRoughly",
            mis, "% of the data is missing.\n")
    } else {
        cat("A genome object encoding", nm, "markers across", nc,
            "chromosomes, distributed:\n", dtab, "\n")
    }
}
##' @describeIn methods Plot method for `genome`
plot.genome <- function(genome, chrLens, epch = c(1,2,6,0),
                        elevs = c("0 0", "0 1", "1 0", "1 1"),
                        scale = FALSE, add.legend = FALSE, ...) {
    if (missing(chrLens)) { # take largest location as max length
        chrLens <- sapply(genome$location, max)
    }
    if (scale) {
        pos <- mapply(function(locs, mx) locs/max,
                      genome$location, chrLens, SIMPLIFY = FALSE)
        xlab <- "Relative location"
    } else {
        pos <- genome$location
        xlab = "Location"
    }
    ht <- seq(0, 1, length.out = length(pos)) # heights
    encInds <- c(0, cumsum(table(genome$chr))) # encodind rows
    encFac <- factor(paste(genome$encoding[,1], genome$encoding[,2]),
                     levels = elevs) # encoding factor
    plot(NA, xlim = c(0, max(chrLens)), ylim = range(ht), yaxt = "n",
         xlab = xlab, ylab = "Chromosome", ...) # plot area
    axis(2, at = ht, labels = levels(genome$chr)) # custom y axis
    for (ii in seq_along(pos)) { # add chromosome lines
        lines(x = c(0, chrLens[ii]), y = c(ht[ii], ht[ii]))
        points(x = pos[[ii]], y = rep(ht[ii], length(pos[[ii]])),
               pch = epch[unclass(encFac[(encInds[ii]+1):(encInds[ii+1])])])
    }
    if (add.legend) {
        legend(x = "topright", legend = elevs, pch = epch,
               title = "Encoding")
    }
}
