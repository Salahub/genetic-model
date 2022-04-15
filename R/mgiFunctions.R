source("simulationFunctions.R")

mgiUrl <- "http://www.informatics.jax.org/downloads/reports/"
mgiPNames <- c(cope.jenk = "MGI_Copeland-Jenkins_Panel.rpt",
               eucib.bsb = "MGI_EUCIB_BSB_Panel.rpt",
               eucib.bss = "MGI_EUCIB_BSS_Panel.rpt",
               jax.bsb = "MGI_JAX_BSB_Panel.rpt",
               jax.bss = "MGI_JAX_BSS_Panel.rpt",
               jax.mutbcb = "MGI_JAX_Mouse_Mutant_Resource_BCB_Panel.rpt",
               jax.mutbss = "MGI_JAX_Mouse_Mutant_Resource_BSS_Panel.rpt",
               koz.fvc58 = "MGI_Kozak_FvC58_Panel.rpt",
               koz.fvspr = "MGI_Kozak_FvSpr_Panel.rpt",
               koz.skive = "MGI_Kozak_Skive_Panel.rpt",
               mit = "MGI_MIT_Panel.rpt",
               reev.c16 = "MGI_Reeves_Chr_16_Panel.rpt",
               seldin = "MGI_Seldin_Panel.rpt",
               ucla.bsb = "MGI_UCLA_BSB_Panel.rpt")
mgiMNames <- c("MRK_List1.rpt", "MRK_List2.rpt")

## read in MGI data (http://www.informatics.jax.org/)
readMGIrpt <- function(file) {
    raw <- scan(file, what = character(), sep = "\n")
    leg <- which(raw == "Legend:")
    lenHead <- leg + 4
    if (length(leg) == 0) {
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
         markers = data.frame(colDesc, ref = refs[-1]),
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

## a correlation helper
mgiCorrelation <- function(panel, use = "pairwise.complete.obs",
                           method = "pearson") {
    srtd <- panel$data[, order(panel$markers$chr)] # group chromosomes
    dotDrop <- apply(srtd, 2, # turn dots to NA
                     function(mrk) {
                         temp <- mrk
                         temp[temp == "."] <- NA
                         temp
                     })
    numer <- apply(dotDrop, 2, function(mrk) unclass(factor(mrk))-1)
    cor(numer, use = use, method = method) # NOT pos def
}

## and one for the theoretical correlation
mgiTheory <- function(panel, setting) {
    if (setting == "unclear") setting <- "backcross"
    pos <- split(panel$markers$cMs, panel$markers$chr)
    diffs <- lapply(pos, diff) # adjacent distances
    theoryCor(diffs, setting = setting)
}

## drop bad panelists
mgiDropZeroPanelist <- function(panel) {
    bad <- apply(panel$data, 1, function(row) all(row == "."))
    panel$data <- panel$data[!bad,]
    panel
}

## drop bad markers
mgiDropBadMarker <- function(panel) {
    bad <- apply(panel$data, 2, function(col) any(col == "."))
    goodMarkers <- names(panel$data)[!bad]
    list(summary = panel$summary,
         markers = panel$markers[panel$markers$symbol %in% goodMarkers,],
         data = panel$data[, goodMarkers])
}

## simulate a panel based on a setting and cM distances
simulateMGI <- function(marks, npop, chrOrd, reps = 1000,
                        setting = c("backcross", "intercross")) {
    if (missing(chrOrd)) {
        chroms <- split(marks$cMs, marks$chr)
    } else {
        chroms <- split(marks$cMs, marks$chr)[chrOrd]
    }
    setting <- match.arg(setting)
    M <- abiogenesis(length(chroms), sapply(chroms, length),
                     dists = lapply(chroms, diff), allele = 1)
    F <- abiogenesis(length(chroms), sapply(chroms, length),
                     dists = lapply(chroms, diff), allele = 0)
    F1 <- sex(M, F)
    allcors <- vector("list", reps)
    if (setting == "intercross") {
        for (ii in 1:reps) {
            allcors[[ii]] <- popCorrelation(replicate(npop,
                                                    sex(F1, F1),
                                                    simplify = FALSE))
            if (ii %% 100 == 0) cat("\r -- Simulated population ", ii)
        }
    } else if (setting == "backcross") {
        for (ii in 1:reps) {
            allcors[[ii]] <- popCorrelation(replicate(npop,
                                                    sex(F1, M),
                                                    simplify = FALSE))
            if (ii %% 100 == 0) cat("\r -- Simulated population ", ii)
        }
    } else {
        stop("Setting unknown")
    }
    allcors
}

## suppress zeros and reconstruct a correlation matrix
zeroEigSuppress <- function(sigma) {
    decom <- eigen(sigma)
    eigs <- decom$values
    eigs[eigs < 0] <- 0 # zero negatives
    recon <- decom$vectors %*% diag(eigs) %*% t(decom$vectors)
    recon
}

## compute mgi correlation and report population mean for marker pairs
mgiCorrMeans <- function(panel, use = "pairwise.complete.obs",
                         method = "pearson") {
    dropLowTri <- function(mat) { # keep upper entries
        mat[!upper.tri(mat)] <- NA
        mat
    }
    srtd <- panel$data[, order(panel$markers$chr)] # group chromosomes
    dotDrop <- apply(srtd, 2, # turn dots to NA
                     function(mrk) {
                         temp <- mrk
                         temp[temp == "."] <- NA
                         temp
                     })
    numer <- apply(dotDrop, 2, function(mrk) unclass(factor(mrk)))
    splt <- lapply(split(panel$markers$symbol, panel$markers$chr),
                   function(mrks) numer[, mrks]) # intra chromosome
    cors <- lapply(splt, cor, use = use, method = method)
    cors <- lapply(cors, dropLowTri)
    means <- lapply(splt, colMeans) # means by marker
    pwmns <- lapply(means, function(mn) outer(mn, mn, pmax)) # ref mean
    pwmns <- lapply(pwmns, dropLowTri)
    list(corr = cors, means = pwmns)
}

## subset a panel to specified markers
mgiSubset <- function(panel, symbols) {
    list(summary = panel$summary,
         markers = panel$markers[panel$markers$symbol %in% symbols,],
         data = panel$data[ ,symbols])
}
