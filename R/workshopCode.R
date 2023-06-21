## recreate the code from the workshop here
library(toyGenomeGenR)


## FUNCTIONS #########################################################
## custom plotting function with narrow margins
narrowPlot <- function(xgrid, ygrid, main = "", xlab = "", ylab = "",
                       xticks = xgrid, yticks = ygrid,
                       mars = c(2.1, 2.1, 1.1, 1.1),
                       xlim = range(xgrid), ylim = range(ygrid),
                       addGrid = TRUE, ...) {
    par(mar = mars) # set narrow margins
    plot(NA, ylim = ylim, xlim = xlim, xaxt = 'n', xlab = "",
         yaxt = 'n', ylab = "", main = "", ...)
    ## add labels
    mtext(main, side = 3, line = 0, cex = 0.8) # main
    mtext(ylab, side = 2, line = 1, cex = 0.8) # ylab
    mtext(xlab, side = 1, line = 1, padj = 0, cex = 0.8) # xlab
    ## add grid lines
    if (addGrid) {
        abline(h = ygrid, v = xgrid, lty = 1,
               col = adjustcolor("gray", alpha.f = 0.4))
    }
    ## and ticks
    mtext(side = 1, at = xgrid, text = "|", line = 0, cex = 0.5,
          padj = -2)
    mtext(text = xticks, at = xgrid, side = 1, cex = 0.8)
    mtext(side = 2, at = ygrid, text = "|", line = 0, cex = 0.5,
          padj = 1)
    mtext(text = yticks, at = ygrid, side = 2, cex = 0.8)
}

## quantile plot of central range
addQuantPoly <- function(mat, qnts = c(0.005, 0.025, 0.25),
                         kseq = exp(seq(-8, 8, by = 0.25)),
                         labpos = c("left", "top"), ...) {
    for (qn in qnts) {
        polygon(log(c(kseq, rev(kseq)), base = 10),
                c(apply(mat, 1, quantile,
                        probs = qn),
                  rev(apply(mat, 1, quantile,
                            probs = 1-qn))),
                col = adjustcolor("gray", 0.25), border = NA)
        if (labpos[1] == "left") {
            xind <- 1
            adj <- c(0, 0.5)
        } else if (labpos[1] == "right") {
            xind <- length(kseq)
            adj <- c(1, 0.5)
        }
        if (labpos[2] == "top") {
            yq <- 1-qn
        } else if (labpos[2] == "bottom") {
            yq <- qn
        }
        text(x = log(kseq[xind], 10),
             y = quantile(mat[xind,], yq),
             labels = 1-2*qn, adj = adj, ...)
    }
    lines(x = log(kseq, 10),
          y = apply(mat, 1, median))
}

## independent pooled chi value
poolChi <- function(p, kap) {
    pchisq(sum(qchisq(p, df = kap, lower.tail = FALSE)),
           df = length(p)*kap, lower.tail = FALSE)
}

## try the method of moments gamma match from ferrari's arxiv note on
## chi-squared sums
gammaApprox <- function(qs, sigma, kap) {
    diag(sigma) <- 0 # for later sum
    M <- length(qs)
    u <- 2*(1 + 2*sum(sigma*kap)/(M*kap))
    pgamma(sum(qs), shape = (M*kap)/u, scale = u,
           lower.tail = FALSE)
}
## can also try a satterthwaite approximation
satterApprox <- function(qs, sigma, kap) {
    diag(sigma) <- 0
    M <- length(qs)
    ex <- M*kap
    varx <- 2*M*kap + 2*kap*sum(sigma)
    f <- 2*ex^2/varx
    c <- varx/(2*ex)
    pchisq(sum(qs)/c, df = f, lower.tail = FALSE)
}
## what does this do to pooled p-values?
poolChiDep <- function(ps, kap, sigma,
                       method = c("gamma", "satterthwaite")) {
    chis <- qchisq(ps, df = kap, lower.tail = FALSE)
    method <- match.arg(method)
    if (method == "gamma") {
        gammaApprox(chis, sigma, kap)
    } else if (method == "satterthwaite") {
        satterApprox(chis, sigma, kap)
    } else {
        stop("Method must be one of 'gamma' and 'satterthwaite'")
    }
}

## filter bad markers
mgiDropBadMarker <- function(panel, prop = 0) {
    bad <- apply(do.call(cbind, panel$encodings), 1,
                 function(row) {
                     sum(is.na(row))/(2*length(panel$encodings)) > prop
                 })
    subsetPopulation(panel, markInd = !bad)
}

## write a simple function to generate pesudo-pairs based on prec and
## repeat a test of association with a trait
simTest <- function(prec, nsim = 1e5, npop = 1e2, probs = c(0.5, 0.5),
                    trait = rnorm(npop, sd = 0.3), test = kruskal) {
    r <- rank(trait) # ranking of traits
    ties <- table(trait) # ties
    init <- matrix(rep(sample(c(0,1), npop*nsim, replace = TRUE,
                              probs), 2),
                   nrow = 2, ncol = npop*nsim, byrow = TRUE)
    rec <- runif(npop*nsim) < 2*prec # recom
    init[2,rec] <- sample(c(0,1), sum(rec), replace = TRUE, probs)
    vals <- array(init, dim = c(2, npop, nsim)) # return array
    ps <- apply(vals, c(1,3), test, x = r, n = npop, TIES = ties)
    ps
}
## custom kruskal test (gutted version of kruskal.test)
kruskal <- function(x, g, n, TIES) {
    g <- factor(g)
    k <- nlevels(g)
    ##r <- rank(x)
    r <- x
    ##TIES <- table(x)
    STATISTIC <- sum(tapply(r, g, sum)^2/tapply(r, g, length))
    STATISTIC <- ((12 * STATISTIC/(n * (n + 1)) - 3 * (n + 1))/(1 -
        sum(TIES^3 - TIES)/(n^3 - n)))
    PARAMETER <- k - 1L
    pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
}

## different tester for random splits
simSplit <- function(prec, nsim = 1e5, npop = 1e2, probs = c(0.5, 0.5),
                     trait = rnorm(npop, sd = 0.3)) {
    r <- rank(trait) # ranking of traits
    init <- matrix(rep(sample(c(0,1), npop*nsim, replace = TRUE,
                              probs), 2),
                   nrow = 2, ncol = npop*nsim, byrow = TRUE)
    rec <- runif(npop*nsim) < 2*prec # recom
    init[2,rec] <- sample(c(0,1), sum(rec), replace = TRUE, probs)
    vals <- array(init, dim = c(2, npop, nsim)) # return array
    ps <- apply(vals, 3, randomsplit, x = r, n = npop)
    ps
}
## and a bin splitting p-value
randomsplit <- function(x, g, n, lim = 10) {
    brk1 <- sample(lim:(n-lim), 1)
    if (brk1-2*lim <= 0) {
        brk2 <- sample(seq(brk1+lim, n-lim, by = 1), 1)
        brks <- c(0, brk1, brk2, n)
    } else if (n-brk1-2*lim <= 0) {
        brk2 <- sample(seq(lim, brk1-lim, by = 1), 1)
        brks <- c(0, brk2, brk1, n)
    } else {
        brk2 <- sample(c(seq(lim, brk1-lim, by = 1),
                         seq(brk1+lim, n-lim, by = 1)), 1)
        brks <- c(0, min(brk1, brk2), max(brk1, brk2), n)
    }
    xc <- cut(x, breaks = c(0, brk1, brk2, n))
    obs1 <- table(xc, g[1,]) # counts
    obs2 <- table(xc, g[2,])
    ex1 <- outer(diff(brks), table(g[1,])/n) # expected densities
    ex2 <- outer(diff(brks), table(g[2,])/n) # expected densities
    pchisq(c(sum(obs1^2/ex1), sum(obs2^2/ex2)) - n,
           df = prod(dim(obs1) - 1), lower.tail = FALSE)
}

## split univariate ranks in a single group
splitCol <- function(x, n, lim = 10) {
    nx <- length(x) # size of margin
    xs <- c(0, sort(x), n+1) # sorted order, add limits
    dif <- diff(xs) # size of gaps indicate maximum
    ind1 <- which.max(dif) # index of first split
    splt1 <- max(c(min(c(xs[ind1], n - lim)), lim))
    if (splt1 >= n - 2*lim | ind1 >= nx) { # can only split below
        ind2 <- which.max(dif[1:(ind1-1)])
        splt2 <- max(c(xs[ind2], lim))
        return(list(ex = c(splt2, splt1-splt2, n-splt1)/n,
                    obs = c(ind2-1, ind1-ind2, nx-ind1+1)))
    } else if (splt1 <= 2*lim | ind1 <= 1) { # can only split above
        ind2 <- which.max(dif[(ind1+1):(nx+1)])
        splt2 <- min(c(xs[ind1+ind2], n-lim))
        return(list(ex = c(splt1, splt2-splt1, n-splt2)/n,
                    obs = c(ind1-1, ind2, nx-ind2-ind1+1)))
    } else {
        ind2 <- which.max(dif[(ind1+1):(nx+1)])
        splt2 <- min(c(xs[ind1+ind2], n-lim))
        ind3 <- which.max(dif[1:(ind1-1)])
        splt3 <- max(c(xs[ind3], lim))
        return(list(ex = c(splt3, splt1-splt3, splt2-splt1,
                           n-splt2)/n,
                    obs = c(ind3-1, ind1-ind3, ind2, nx-ind2-ind1+1)))
    }
}

## bin splitting p-value unpaired
unimaxsplit <- function(x, g, n, lim = 10) {
    g1 <- splitCol(x[g == 1], n = n, lim = lim)
    g2 <- splitCol(x[g == 2], n = n, lim = lim)
    obs <- c(g1$obs, g2$obs) # counts
    ex <- c((g1$ex)*sum(g == 1), g2$ex*sum(g == 2)) # expected
    pchisq(sum(obs^2/ex) - n,
           df = prod(dim(obs) - 1), lower.tail = FALSE)
}

## and a bin splitting p-value
unirandsplit <- function(x, g, n, lim = 10) {
    brk1 <- sample(lim:(n-lim), 1)
    if (brk1-2*lim <= 0) {
        brk2 <- sample(seq(brk1+lim, n-lim, by = 1), 1)
        brks <- c(0, brk1, brk2, n)
    } else if (n-brk1-2*lim <= 0) {
        brk2 <- sample(seq(lim, brk1-lim, by = 1), 1)
        brks <- c(0, brk2, brk1, n)
    } else {
        brk2 <- sample(c(seq(lim, brk1-lim, by = 1),
                         seq(brk1+lim, n-lim, by = 1)), 1)
        brks <- c(0, min(brk1, brk2), max(brk1, brk2), n)
    }
    xc <- cut(x, breaks = c(0, brk1, brk2, n))
    obs1 <- table(xc, g) # counts
    ex1 <- outer(diff(brks), table(g)/n) # expected densities
    pchisq(sum(obs1^2/ex1) - n,
           df = prod(dim(obs1) - 1), lower.tail = FALSE)
}

## SCRIPT ############################################################
## null quantiles
nullQuants <- readRDS("curveMinQuantiles.Rds")

## palette and breaks
pal <- colorRampPalette(c("steelblue", "white", "firebrick"))(40)
breaks <- seq(-1, 1, length.out = 41)

## convert marker correlation to chi covariance
## trait parameters and kappa values
kseq <- exp(seq(-8, 8, by = 0.1))
## define range of correlations to compute over
corrs <- seq(0, 1, by = 0.02)
## for each, convert to a probability of recombination
prec <- (1 - corrs)/2

if (!("chiCorrs.Rds" %in% list.files())) {
    set.seed(8081326)
    ## apply this to the range of recombination probabilities
    precPvals <- replicate(5,
                           simplify2array(lapply(prec, simSplit,
                                                 nsim = 1e4,
                                                 npop = 2.5e2)),
                           simplify = FALSE)
    ## get correlation of chi transform by kappa
    chiCorrs <- simplify2array(lapply(precPvals,
                        function(pairs) {
                             sapply(kseq,
                                    function(k) {
                                        apply(qchisq(pairs,
                                                     k,
                                                     lower.tail = FALSE),
                                              3,
                                              function(mat) {
                                                  cor(t(mat))
                                              })[2,]
                                        })
                                      }))
    ## save in single data frame
    chiCordf <- data.frame(expand.grid(zcor = corrs,
                                       logkap = log(kseq),
                                       rep = 1:length(precPvals)),
                           chicor = c(chiCorrs))
    ## save this as an RDS
    saveRDS(chiCordf, "chiCorrs.Rds")
} else {
    chiCordf <- readRDS("chiCorrs.Rds")
}
## fit 5th order polynomials
chiCorMods <- lapply(log(kseq), function(k) {
    lm(chicor ~ I(zcor^2) + I(zcor^4) + I(zcor^6) + I(zcor^8) +
           I(zcor^10), data = chiCordf,
       subset = abs(chiCordf$logkap - k)<0.0001)#$coefficients
})

## use these to convert a rho matrix to r for a kappa
convertRho <- function(sigma, kapInd, models = chiCorMods) {
    matrix(predict(chiCorMods[[kapInd]],
                   newdata = data.frame(zcor = c(sigma))),
           ncol = ncol(sigma))
}

## plot some example curves
targKap <- 0
exDat <- chiCordf[abs(chiCordf$logkap - targKap) < 0.0005,]
png(paste0("rhoVSr", targKap, ".png"), width = 2.5, height = 2.5,
    units = "in", res = 480)
narrowPlot(xgrid = seq(0, 1, by = 0.25),
           ygrid = seq(0, 1, by = 0.25),
           xlab = expression(rho[ij]), ylab = expression(r[ij]))
points(exDat$zcor, exDat$chicor, cex = 0.5)
dev.off()

## how good are these fits?
png("rhoVSrCoefVar.png", width = 3, height = 3, units = "in",
    res = 480)
narrowPlot(xgrid = seq(-3, 3, by = 1.5), xlim = c(-3.5, 3.5),
           ygrid = seq(0.6, 1, by = 0.1),
           xlab = expression(log[10]~{"("~kappa~")"}),
           ylab = expression(R^2))
lines(log(kseq, 10),
      sapply(chiCorMods, function(mod) summary(mod)$r.squared))
dev.off()

## load the jax bsb data
data(jax_bsb)
jax_bsb
## drop NA markers
jax_bsb <- mgiDropBadMarker(jax_bsb)
jax_bsb
## get middle marker in each chromosome
midChr <- function(chroms) ceiling(cumsum(table(chroms)) -
                                   table(chroms)/2)
## extract markers from jax_bsb
jax_bsbChr <- subsetPopulation(jax_bsb,
                               midChr(jax_bsb$chromosome))
png("workshopIndependence.png", width = 3, height = 3, res = 480,
    units = "in")
par(mar = c(1.1, 1.1, 1.1, 1.1))
corrImg(popCorrelation(jax_bsbChr), col = pal, breaks = breaks,
        xaxt = "n", yaxt = "n")
dev.off()

## parameters
beta <- 1
sd <- 0.3
## define vector of counts
ns <- c(1, 5, 10, 15, 20)
## compute scores
jax_scores <- sapply(jax_bsbChr$encodings, scoreAdditive)
## repeat generation and testing many times
nsim <- 1e3
nulls <- curves <- array(dim = c(length(ns), length(kseq), nsim))
nullpmat <- pmat <- array(dim = c(length(ns), nrow(jax_scores),
                                  nsim))
set.seed(11200905)
for (ii in 1:nsim) {
    for (jj in seq_along(ns)) {
        trait <- beta*colMeans(jax_scores[1:ns[jj],,#sample(1:20, ns[jj]),,
                                          drop = FALSE]) +
            rnorm(94, sd = sd)
        ps <- apply(jax_scores, 1,
                    function(g) unirandsplit(x = rank(trait),
                                             g, n = 94))
        null <- apply(jax_scores, 1,
                      function(g) unirandsplit(x = rank(sample(trait)),
                                               g, n = 94))
        pools <- sapply(kseq, poolChi, p = ps)
        nullpl <- sapply(kseq, poolChi, p = null)
        curves[jj, , ii] <- pools
        nulls[jj, , ii] <- nullpl
        pmat[jj, , ii] <- ps
        nullpmat[jj, , ii] <- null
    }
    if (ii %% 50 == 0) cat("\r Done", ii, "of", nsim)
}

## plot observed curves against quantiles
ind <- 5
mat <- curves
png(paste0("wrkshpIndepCurves", ind, ".png"), width = 2.5,
    height = 2.5, units = "in", res = 480)
quantLevs <- c("5%", "1%", "0.1%")
narrowPlot(xgrid = seq(-3, 3, by = 1.5), xlim = c(-3.5, 3.5),
           ygrid = seq(-16, 0, by = 4),
           xlab = expression(log[10]~{"("~kappa~")"}),
           ylab = expression(log[10]~{"("~p~")"}))
abline(h = log(nullQuants[quantLevs, as.character(100)], 10),
       lty = 2, col = "firebrick")
text(x = rep(3.95, 3), labels = quantLevs, cex = 0.6, xpd = NA,
     y = log(nullQuants[quantLevs, as.character(100)], 10),
     adj = c(0.2, 0.5))
addQuantPoly(log(mat[ind,,], 10), kseq = kseq, cex = 0.6,
             labpos = c("right", "bottom"))
dev.off()

## change dependence to a block structure
## get 10 markers over the first two chromosomes
k <- 10
chrInds <- cumsum(table(jax_bsb$chromosome))[1:2]
steps <- diff(c(0,chrInds))/k
inds <- floor(rep(c(0, chrInds[1:(length(chrInds)-1)]), each = k) +
              (1:k)*rep(steps, each = k))
## extract markers from jax_bsb
jax_bsbBlock <- subsetPopulation(jax_bsb, markInd = inds)
jax_blockCor <- popCorrelation(jax_bsbBlock)
png("workshopBlock.png", width = 3, height = 3, res = 480,
    units = "in")
par(mar = c(1.1, 1.1, 1.1, 1.1))
corrImg(jax_blockCor, col = pal, breaks = breaks, xaxt = "n",
        yaxt = "n")
addChromosomeLines(jax_bsbBlock)
dev.off()

## the same analysis
ns <- c(1/length(inds), 0.25, 0.5, 0.75, 1)*length(inds)
beta <- 1
sd <- 0.3
jax_scores <- sapply(jax_bsbBlock$encodings, scoreAdditive)
## repeat generation and testing many times
nsim <- 1e3
kInds <- 1:length(kseq)
#kInds <- which(sapply(chiCorMods,
#                      function(mod) summary(mod)$r.squared) > 0.98)
curves <- array(dim = c(length(ns), length(kseq), nsim))
pmat <- array(dim = c(length(ns), nrow(jax_scores), nsim))
adjcurves <- array(dim = c(length(ns), length(kInds), nsim))
set.seed(14521105)
for (ii in 1:nsim) {
    for (jj in seq_along(ns)) {
        trait <- beta*colMeans(jax_scores[1:ns[jj],,
                                          drop = FALSE]) +
            rnorm(94, sd = sd)
        ps <- apply(jax_scores, 1,
                    function(g) unirandsplit(x = rank(trait), g,
                                             n = 94))
        pools <- sapply(kseq, poolChi, p = ps)
        padj <- sapply(
            kInds,
            function(ii) poolChiDep(p = ps,
                                    kap = kseq[ii],
                                    sigma = convertRho(jax_blockCor,
                                                       ii),
                                    method = "satterthwaite"))
        curves[jj, , ii] <- pools
        adjcurves[jj, , ii] <- padj
        pmat[jj, , ii] <- ps
    }
    if (ii %% 50 == 0) cat("\r Done", ii, "of", nsim)
}

## plot observed curves against quantiles
ind <- 3
version <- "adj"
if (version == "adj") mat <- adjcurves else mat <- curves
png(paste0("wrkshpDepCurves", version, ind, ".png"),
    height = 2.5, width = 2.5, units = "in", res = 480)
quantLevs <- c("5%", "1%", "0.1%")
narrowPlot(xgrid = seq(-3, 3, by = 1.5), xlim = c(-3.5, 3.5),
           ygrid = seq(-30, 0, by = 6),
           xlab = expression(log[10]~{"("~kappa~")"}),
           ylab = expression(log[10]~{"("~p~")"}))
abline(h = log(nullQuants[quantLevs, as.character(100)], 10),
       lty = 2, col = "firebrick")
text(x = rep(3.95, 3), labels = quantLevs, cex = 0.6, xpd = NA,
     y = log(nullQuants[quantLevs, as.character(100)], 10),
     adj = c(0.2, 0.5))
addQuantPoly(log(mat[ind,,], 10), kseq = kseq, cex = 0.6,
             labpos = c("left", "bottom"))
dev.off()

## plot p-value distributions by ind
xpos <- rep(c(jax_bsbBlock$location[[1]],
              100 + jax_bsbBlock$location[[2]]), nsim)
plot(x = xpos, xaxt = "n",
     y = log(c(pmat[ind,,]), 10), type = "n",
     ylab = expression(paste(log[10], "(", p, ")")),
     xlab = "Marker position (cM)", ylim = c(-17, 0))
axis(1, at = seq(0, 200, by = 50),
     labels = c(0, 50, 0, 50, 100))
for (x in unique(xpos)) {
    boxplot(log(c(pmat[ind,,])[xpos == x], 10), add = TRUE,
            at = x, boxwex = 5,
            ann = FALSE, pars = list(xaxt = "n"),
            col = "gray90")
}
abline(v = 100)
abline(h = log(0.05, 10), lty = 2, col = "firebrick")
abline(h = seq(0, -20, by = -5), v = seq(0, 200, by = 50),
       col = adjustcolor("gray", 0.4))


## dependence in a single chromosome
## get 20 markers over the first chromosome
k <- 20
chrInds <- cumsum(table(jax_bsb$chromosome))[1]
steps <- diff(c(0,chrInds))/k
inds <- floor((1:k)*rep(steps, each = k))
## extract markers from jax_bsb
jax_bsbSingle <- subsetPopulation(jax_bsb, markInd = inds)
jax_bsbSingCor <- popCorrelation(jax_bsbSingle)
corrImg(jax_bsbSingCor, col = pal, breaks = breaks,
        xaxt = "n", yaxt = "n")

## the same analysis
ns <- c(1/length(inds), 0.25, 0.5, 0.75, 1)*length(inds)
beta <- 1
sd <- 0.3
jax_scores <- sapply(jax_bsbSingle$encodings, scoreAdditive)
## repeat generation and testing many times
nsim <- 1e3
kseq <- exp(seq(-10, 10, by = 0.5))
nulls <- curves <- array(dim = c(length(ns), length(kseq), nsim))
nullpmat <- pmat <- array(dim = c(length(ns), nrow(jax_scores), nsim))
set.seed(7521105)
for (ii in 1:nsim) {
    for (jj in seq_along(ns)) {
        trait <- beta*colMeans(jax_scores[1:ns[jj],,
                                          drop = FALSE]) +
            rnorm(94, sd = sd)
        ps <- apply(jax_scores, 1,
                    function(g) kruskal.test(x = trait, g)$p.value)
        null <- apply(jax_scores, 1,
                      function(g) kruskal.test(x = sample(trait),
                                               g)$p.value)
        pools <- sapply(kseq, poolChi, p = ps)
        nullpl <- sapply(kseq, poolChi, p = null)
        curves[jj, , ii] <- pools
        nulls[jj, , ii] <- nullpl
        pmat[jj, , ii] <- ps
        nullpmat[jj, , ii] <- null
    }
    if (ii %% 50 == 0) cat("\r Done", ii, "of", nsim)
}

## plot observed curves and quantiles
ind <- 4
mat <- curves
plot(NA, xlim = range(log(kseq, base = 10)),
     ylim = c(-120,0),
            #range(c(log(jax_mgPooled, 10),
            #        log(jax_pgPooled, 10))),
     ylab = expression(paste(log[10], "(p)")),
     xlab = expression(paste(log[10], "(", kappa, ")")))
abline(h = log(0.05, 10), lty = 2, col = "firebrick")
abline(h = seq(0, -120, by = -40), v = seq(-4, 4, by = 2),
       col = adjustcolor("gray", 0.4))
lines(log(kseq, 10), log(apply(mat[ind,,], 1, median), 10))
polygon(c(log(kseq, 10), log(rev(kseq), 10)),
        c(log(apply(mat[ind,,], 1, quantile, probs = 0.025), 10),
          log(rev(apply(mat[ind,,], 1, quantile, probs = 0.975)), 10)),
        col = adjustcolor("gray", 0.2), border = NA)
polygon(c(log(kseq, 10), log(rev(kseq), 10)),
        c(log(apply(mat[ind,,], 1, quantile, probs = 0.25), 10),
          log(rev(apply(mat[ind,,], 1, quantile, probs = 0.75)), 10)),
        col = adjustcolor("gray", 0.2), border = NA)

## can we improve this performance (i.e. bring curves up)?
depcurves <- array(dim = c(length(ns), length(kseq), nsim))
set.seed(11200905)
for (ii in 1:nsim) {
    for (jj in seq_along(ns)) {
        trait <- beta*colMeans(jax_scores[1:ns[jj],,
                                          drop = FALSE]) +
            rnorm(94, sd = sd)
        ps <- apply(jax_scores, 1,
                    function(g) kruskal.test(x = trait, g)$p.value)
        pools <- sapply(kseq, poolChiDep, p = ps,
                        sigma = jax_bsbSingCor,
                        method = "satterthwaite")
        depcurves[jj, , ii] <- pools
    }
    if (ii %% 50 == 0) cat("\r Done", ii, "of", nsim)
}

## plot observed curves against quantiles
ind <- 4
mat <- depcurves
plot(NA, xlim = range(log(kseq, base = 10)),
     ylim = c(-15,0),
            #range(c(log(jax_mgPooled, 10),
            #        log(jax_pgPooled, 10))),
     ylab = expression(paste(log[10], "(p)")),
     xlab = expression(paste(log[10], "(", kappa, ")")))
abline(h = log(0.05, 10), lty = 2, col = "firebrick")
abline(h = seq(0, -15, by = -5), v = seq(-4, 4, by = 2),
       col = adjustcolor("gray", 0.4))
lines(log(kseq, 10), log(apply(mat[ind,,], 1, median), 10))
polygon(c(log(kseq, 10), log(rev(kseq), 10)),
        c(log(apply(mat[ind,,], 1, quantile, probs = 0.025), 10),
          log(rev(apply(mat[ind,,], 1, quantile, probs = 0.975)), 10)),
        col = adjustcolor("gray", 0.2), border = NA)
polygon(c(log(kseq, 10), log(rev(kseq), 10)),
        c(log(apply(mat[ind,,], 1, quantile, probs = 0.25), 10),
          log(rev(apply(mat[ind,,], 1, quantile, probs = 0.75)), 10)),
        col = adjustcolor("gray", 0.2), border = NA)

## could be interpreted with a number of observations:
## 1: dependence makes detection of "anything" easier
## 2: dependence spreads evidence
xpos <- rep(jax_bsbSingle$location$`1`, nsim)
ind <- 1
plot(x = xpos,
     y = log(c(pmat[ind,,]), 10), type = "n",
     ylab = expression(paste(log[10], "(", p, ")")),
     xlab = "Marker position (cM)", ylim = c(-17, 0))
for (x in unique(xpos)) {
    boxplot(log(c(pmat[ind,,])[xpos == x], 10), add = TRUE,
            at = x, boxwex = 3,
            ann = FALSE, pars = list(xaxt = "n"),
            col = "gray90")
}
abline(h = log(0.05, 10), lty = 2, col = "firebrick")
abline(h = seq(0, -20, by = -5), v = seq(0, 100, by = 20),
       col = adjustcolor("gray", 0.4))
