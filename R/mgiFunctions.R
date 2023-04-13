#source("simulationFunctions.R")

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
