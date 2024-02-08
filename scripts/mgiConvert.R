## load functions, set image output directory
library(toyGenomeGenR)

## FUNCTIONS TO LOAD, PROCESS MGI DATA ###############################
mgiUrl <- "https://www.informatics.jax.org/downloads/reports/"
mgiPNames <- c(cope_jenk = "MGI_Copeland-Jenkins_Panel.rpt",
               eucib_bsb = "MGI_EUCIB_BSB_Panel.rpt",
               eucib_bss = "MGI_EUCIB_BSS_Panel.rpt",
               jax_bsb = "MGI_JAX_BSB_Panel.rpt",
               jax_bss = "MGI_JAX_BSS_Panel.rpt",
               jax_mutbcb = "MGI_JAX_Mouse_Mutant_Resource_BCB_Panel.rpt",
               jax_mutbss = "MGI_JAX_Mouse_Mutant_Resource_BSS_Panel.rpt",
               koz_fvc58 = "MGI_Kozak_FvC58_Panel.rpt",
               koz_fvspr = "MGI_Kozak_FvSpr_Panel.rpt",
               koz_skive = "MGI_Kozak_Skive_Panel.rpt",
               mit = "MGI_MIT_Panel.rpt",
               reev_c16 = "MGI_Reeves_Chr_16_Panel.rpt",
               seldin = "MGI_Seldin_Panel.rpt",
               ucla_bsb = "MGI_UCLA_BSB_Panel.rpt")
mgiMNames <- c("MRK_List2.rpt")
chrOrder <- c(as.character(1:19), "X") # control factor level ordering

## make a well-formed phenome API request from locations, chromosome,
## dataset, strains
makeAPIreq <- function(beg, end, chr, dset = "UCLA1",
                       strn = "all") {
    paste0("https://phenome.jax.org/api/snpdata?dataset=", dset,
           "&region=", chr, ":", beg, "-", end, "&strains=",
           strn)
}

## process API output
processAPIout <- function(strings) {
    strings <- gsub("\"|,$", "", strings)
    strings <- gsub(",$", ", ", strings)
    splt <- strsplit(strings, ",")
    mat <- do.call(rbind, splt[-1])
    mat[mat == " "] <- ""
    colnames(mat) <- splt[[1]]
    mat
}

## read in MGI data (http://www.informatics.jax.org/)
readMGIrpt <- function(file) {
    raw <- scan(file, what = character(), sep = "\n")
    leg <- which(raw == "Legend:")
    lenHead <- leg + 4
    if (length(leg) == 0) { # separate legend
        leg <- which(grepl("^CHR", raw))[1]
        lenHead <- leg
    }
    desc <- paste(raw[1:lenHead], collapse = "\n") # data description
    dat <- raw[(lenHead+1):length(raw)] # actual data
    refPos <- regexec("\\tJ\\:[0-9]+(?:, J\\:[0-9]+){0,4}", dat)
    refs <- sapply(regmatches(dat, refPos), # extract references
                   function(el) {
                       if (length(el) == 0) {
                           ""
                       } else gsub("\\t", "", el)})
    data <- unlist(regmatches(dat, refPos, invert = TRUE)) # remove refs
    mat <- do.call(rbind, strsplit(data[data != ""], "\\t"))
    rwnms <- mat[1, -(1:3)] # animal numbers/ids
    colnms <- mat[-1, 3] # symbol field
    colDesc <- mat[-1, 1:3] # symbol details
    colnames(colDesc) <- c("chr", "mgiid", "symbol")
    data <- t(mat[-1,-(1:3)])
    rownames(data) <- rwnms
    colnames(data) <- colnms # final data formatting
    list(summary = desc,
         markers = data.frame(chr = factor(colDesc[, "chr"],
                                           levels = chrOrder),
                              colDesc[,c("mgiid", "symbol")],
                              ref = refs[-1]),
         data = as.data.frame(data)) # return everything
}

## process the MGI reference material
readMGIlists <- function(fileList = paste0(mgiUrl, mgiMNames)) {
    lists <- lapply(fileList, scan, what = character(), sep = "\n")
    lists <- lapply(lists, gsub, pattern = "\t$", replacement= "\t.")
    splits <- lapply(lists, strsplit, split = "\t")
    colnms <- splits[[1]][[1]]
    data <- do.call(rbind,
                    lapply(splits, function(splt) do.call(rbind, splt[-1])))
    colnames(data) <- colnms
    as.data.frame(data)
}

## using indices and a reference table, process centiMorgan positions
processcMs <- function(inds, tab) {
    sel <- tab[inds] # take indices from table
    sel[grepl("syntenic", sel)] <- "Inf"
    suppressWarnings(as.numeric(sel)) # warnings by design
}

## filtering panel data by marker distances
filterPanelBycM <- function(panel, locs) {
    locOrd <- order(locs) # order position
    toKeep <- is.finite(locs[locOrd]) # drop NAs and Infs
    outMrk <- data.frame(panel$markers[locOrd[toKeep],],
                         cMs = locs[locOrd[toKeep]])
    outMrk <- outMrk[order(outMrk$chr),] # order chromosome
    list(summary = panel$summary,
         markers = outMrk, data = panel$data[, outMrk$symbol])
}

## processing to convert mgi data to a genome object
mgiBackCrossToGenome <- function(mgiPanel) {
    legendPattern <- "([a-zA-Z0-9]+) indicates allele from ([^\n]+)"
    alleles <- regmatches(mgiPanel$summary, # extract from legend
                          gregexec(legendPattern,
                                   mgiPanel$summary))[[1]]
    allOrd <- order(alleles[3,]) # consistent order
    values <- c(0,1)[allOrd]  # arbitrary values
    back <- regmatches(mgiPanel$summary, # backcross allele
                       regexec("(?<=\\)F1 x )[^\n]+",
                               mgiPanel$summary,
                               perl = TRUE))
    gens <- lapply(1:nrow(mgiPanel$data),
                   function(row) {
                       asGenome(data.frame(mv = as.character(mgiPanel$data[row,]),
                                           pv = alleles[2, match(back,
                                                                 alleles[3,])],
                                           chr = mgiPanel$markers$chr,
                                           pos = mgiPanel$markers$cMs,
                                           check.rows = FALSE,
                                           row.names = names(mgiPanel$data)),
                                alleles = alleles[2, allOrd],
                                values = values)
                   })
    gens
}

##` a helper to create evenly spaced samples between the minimal and
##' maximal marker locations given a matrix of markers with known
##' locations on a single chromosome
regularSample <- function(markers, n = 40) {
    pos <- markers$loc
    rng <- range(markers$loc)
    evenBrks <- seq(rng[1], rng[2], length.out = n)
    nearestInd <- numeric(n)
    minDiff <- rep(Inf, n) # inf values are never minimal
    for (ii in seq_along(pos)) {
        diffs <- abs(evenBrks - pos[ii]) # compare breaks to positions
        closest <- which.min(diffs) # get closest
        if (minDiff[closest] > diffs[closest]) {
            minDiff[closest] <- diffs[closest]
            nearestInd[closest] <- ii
        }
    }
    markers[nearestInd,]
}

## READ/PROCESS MGI DATA #############################################

##' the MGI website, has a number of "panels" of mice generated in
##' experimental work and then measured, typically the crosses used
##' are backcrosses
##' start by loading the MGI provided database of markers (which is
##' rather large)
mgiMarkers <- readMGIlists() # descriptions of all markers

##' filter and sample these by chromosome to pull genotype data
##' by mouse strain
mgiByChrom <- split(mgiMarkers[, c("cM Position", "Marker Symbol",
                                   "genome coordinate start",
                                   "genome coordinate end")],
                    mgiMarkers$Chr)[chrOrder]
mgiByChrom <- lapply(mgiByChrom, # filter by known positions
                     function(df) {
                         names(df) <- c("loc", "sym", "beg", "end")
                         df$loc <- as.numeric(df$loc)
                         df$beg <- as.numeric(df$beg)
                         df$end <- as.numeric(df$end)
                         df[!is.na(df$loc), ]
                     })
##' get regularly sampled markers across each chromosome
mgiRegSamps <- lapply(mgiByChrom[sapply(mgiByChrom, nrow) > 0],
                      regularSample, n = 60)
mgiRegSamps <- cbind(do.call(rbind, mgiRegSamps), # add chromosomes
                     chr = rep(names(mgiRegSamps),
                               times = sapply(mgiRegSamps, nrow)))
mgiRegSamps <- na.omit(mgiRegSamps)

##' standardize interval lengths, decrease length to make reads
##' shorter and increase to make them longer
intLen <- 2e5
mgiRegSamps$end <- mgiRegSamps$beg + intLen
##' make API requests
phenAPIRequests <- with(mgiRegSamps, makeAPIreq(beg, end, chr))
##' pull the data
phenData <- vector(mode = "list", length = length(phenAPIRequests))
names(phenData) <- phenAPIRequests # name with query
hits <- logical(length(phenData)) # did the query at ii succeed?
for (ii in seq_along(phenAPIRequests)) {
    cat(ii, "of", length(phenAPIRequests), ":\n  ")
    curr <- url(phenAPIRequests[ii]) # open connection
    phenData[[ii]] <- scan(curr, sep = "\n", what = "character") #read
    hits[ii] <- length(phenData[[ii]]) > 1
    close(curr) # close connection
}
## process the API pulls
phenDataMats <- lapply(phenData[hits], processAPIout)

##' compute centiMorgans: this requires interpolation, as the API
##' requests work on the scale of base pairs directly, and so the
##' cM values at base pairs surrounding a marker have to be used to
##' approximate the cM position of the marker
mgiCMrate <- diff(mgiRegSamps$loc)/diff(mgiRegSamps$beg)
mgiCMrate[mgiCMrate < 0] <- mgiCMrate[which(mgiCMrate < 0) - 1]
mgiCMrate <- c(mgiCMrate, mgiCMrate[length(mgiCMrate)])
phencMs <- mapply(function(mat, beg, rt, loc) (as.numeric(mat[, "bp38"]) - beg)*rt + loc,
                  phenDataMats, mgiRegSamps$beg[hits],
                  mgiCMrate[hits], mgiRegSamps$loc[hits])
## put it in a single matrix
phenDataFull <- cbind(cMs = unlist(phencMs), do.call(rbind, phenDataMats))
## drop duplicates (some requests overlap)
phenDataFull <- phenDataFull[match(unique(phenDataFull[, "rs"]),
                                   phenDataFull[, "rs"]),]
rownames(phenDataFull) <- NULL
phenDataFull <- phenDataFull[, c("chr", "bp38", "cMs",
                                 colnames(phenDataFull)[4:ncol(phenDataFull)])]

##' next we can consider the panels, first read the data
mgiPanels <- lapply(paste0(mgiUrl, mgiPNames), readMGIrpt)
names(mgiPanels) <- names(mgiPNames)
##' match the marker names back to the marker reference
mgiPanels.mrkr <- lapply(mgiPanels,
                         function(panel) {
                             match(names(panel$data),
                                   mgiMarkers$`Marker Symbol`)
                         })
##' get cM positions
mgiPanels.cMs <- lapply(mgiPanels.mrkr, processcMs,
                        tab = mgiMarkers$`cM Position`)
##' Infs indicate markers localized to a chromosome but without a
##' known cM position, while NAs indicate markers missing from the
##' reference file or marked as unknown there

##' filter the panels by those markers with known positions
mgiFiltered <- mapply(filterPanelBycM, mgiPanels, mgiPanels.cMs,
                      SIMPLIFY = FALSE)

##' identify the type of cross of each using the summaries
mgiPanel.cross <- c("backcross", "backcross", "backcross",
                    "backcross", "backcross", "backcross",
                    "backcross", "unclear", "backcross",
                    "backcross", "intercross", "unclear",
                    "backcross", "backcross")
names(mgiPanel.cross) <- names(mgiFiltered)

##' convert to genome objects
mgiSelect <- names(mgiFiltered)[which(mgiPanel.cross == "backcross")]
mgiGenomes <- vector(mode = "list", length = length(mgiSelect))
names(mgiGenomes) <- mgiSelect # set names
for (nm in mgiSelect) {
    cat(" - Panel : ", nm, "\n")
    mgiGenomes[nm] <- list(mgiBackCrossToGenome(mgiFiltered[[nm]]))
}
##' convert to populations for smaller size
mgiPanelPops <- lapply(mgiGenomes, asPopulation)
##' these can now be saved as RDA files separately
