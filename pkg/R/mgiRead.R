mgiUrl <- "https://www.informatics.jax.org/downloads/reports/"
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
