## MAIN FUNCTIONS ####################################################

##' @title Simple encoding functions
##' @param annotation character vector giving observed annotations
##' @param alleles character vector of possible annotations
##' @param values numeric vector of the same length as alleles giving
##' the value of each allele when encoded
##' @return numeric vector of encodings
##' @author Chris Salahub
encode <- function(annotation, alleles, values) {
    values[match(annotation, alleles)]
}

##' @title Convert a data.frame to a genome object
##' @param df data.frame with columns "mv", "pv", "chr", and "pos"
##' @param alleles list with the same length as df that gives the
##' potential alleles for each row
##' @param values list of identical form to alleles giving numeric
##' values for the encoding of each possible allele
##' @param missing string values to be interpreted as missing values
##' @param encoder function which accepts an annotation and its
##' possible values and returns an encoding
##' @return a genome object with distances given by differences in
##' pos across chromosomes given by chr with encodings given by the
##' values in mv and pv
##' @author Chris Salahub
asGenome <- function(df, alleles, values = c(1,0),
                     missing = c(".", "-"), encoder = encode) {
    ## safety check
    stopifnot(all(c("mv", "pv", "chr", "pos") %in% names(df)))
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
    ## create encoding
    enc <- lapply(1:nrow(df),
                  function(ii) encoder(df[ii, c("mv", "pv")],
                                       alleles[[ii]], values[[ii]]))
    ## collect into genome object
    genome <- list(encoding = enc, chromosome = df$chr,
                   alleles = alleles,
                   distances = lapply(split(df$pos, df$chr), diff))
    class(genome) <- "genome"
    genome
}

##' @title Simulate a genome object
##' @param nmark vector of integers of length nchrom giving the
##' count of measured markers on each chromosome, if named it is
##' assumed the names correspond to the chromosomes, otherwise
##' chromosomes are simply numbered
##' @param markerFun list of functions to simulate marker measurements
##' by chromosome, recycled if it is not as long as nmark
##' @param distFuns list of functions that generate distances between
##' markers by chromosome, recycled if it is not as long as nmark
##' @param alleles vector of possible allele annotations, if missing
##' this is taken to be c("A","a") for every marker site
##' @return a genome object
##' @author Chris Salahub
makeGenome <- function(nmark, alleles, markerFuns = markerHybrid,
                       distFuns = distanceRegular) {
    ## check marker names and generate a chromosome indicator
    if (is.null(names(nmark))) {
        chr <- rep(1:length(nmark), each = nmark)
    } else {
        chr <- rep(names(nmark), each = nmark)
    }
    ## generate alleles if necessary
    if (missing(alleles)) alleles <- generateBiAllelic(nmark)
    ## create encodings
    encoding <- generateEncoding(nmark, markerFuns)
    ## create distances
    distances <- generateDistances(nmark, distFuns)
    ## collect into genome object
    genome <- list(encoding = encoding, chromosome = chr,
                   alleles = alleles, distances = distances)
    class(genome) <- "genome"
    genome
}

##' @title Checking if an object is a valid genome
##' @param genome object to be checked
##' @return TRUE if the object is a genome, otherwise text outlining
##' what failed
##' @author Chris Salahub
checkGenome <- function(genome) {
    isgen <- TRUE # change if false
    ## check names of elements
    if (sort(names(genome)) != c("alleles", "chromosome", "distances",
                                 "encoding")) {
        isgen <- FALSE
        cat("Slots namely incorrectly or missing\n")
    } else if (!is.list(genome$alleles) | # check types of elements
               !is.vector(genome$chromosome) |
               !is.list(genome$distances) |
               !is.matrix(genome$encoding)) {
        isgen <- FALSE
        cat("Slots contain data of the incorrect type\n")
    } else if (length(genome$alleles) != nrow(genome$encoding)) {
        isgen <- FALSE
        cat("Number of alleles does not match encoding\n")
    } else if (!identical(table(genome$chromosome),
                          sapply(genome$distances, length) + 1)) {
        isgen <- FALSE
        cat("Number of distances does not match number of markers\n")
    } else if (length(genome$alleles) != length(genome$chromosome)) {
        isgen <- FALSE
        cat("Distances and alleles imply differing counts of markers")
    }
    isgen        
}
