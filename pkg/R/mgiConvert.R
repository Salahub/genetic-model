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
    cor(numer, use = use, method = method)
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
mgiDropBadMarker <- function(panel, prop = 0) {
    bad <- apply(panel$data, 2,
                 function(col) sum(col == ".")/length(col) > prop)
    goodMarkers <- names(panel$data)[!bad]
    list(summary = panel$summary,
         markers = panel$markers[panel$markers$symbol %in% goodMarkers,],
         data = panel$data[, goodMarkers])
}

## simulate a panel based on a setting and cM distances
simulateMGI <- function(marks, npop, chrOrd, reps = 1000,
                        setting = c("backcross", "intercross"),
                        asArray = FALSE) {
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
    if (asArray) { # place in an appropriate array
        tempcor <- array(0, dim = c(ncom, ncom, nsim))
        for (ii in seq_along(allcors)) tempcor[,,ii] <- allcors[[ii]]
        allcors <- tempcor # for output
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
