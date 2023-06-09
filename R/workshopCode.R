## recreate the code from the workshop here
library(toyGenomeGenR)

## independent pooled chi value
poolChi <- function(p, kap) {
    pchisq(sum(qchisq(p, df = kap, lower.tail = FALSE)),
           df = length(p)*kap, lower.tail = FALSE)
}

## can we simply adjust using the results from papers to get higher
## curves?
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

if (!("chiCorrs.Rds" %in% list.files())) {
    set.seed(8081326)
    ## apply this to the range of recombination probabilities
    precPvals <- replicate(5,
                           simplify2array(lapply(prec, simSplit, nsim = 1e4,
                                                 npop = 2.5e2)),
                           simplify = FALSE)
    ## get correlation of chi transform by kappa
    chiCorrs <- simplify2array(lapply(precPvals,
                                      function(pairs) {
                                          sapply(kseq,
                                                 function(k) {
                                                     apply(qchisq(pairs, k, lower.tail = FALSE),
                                                           3, function(mat) cor(t(mat)))[2,]
                                                 })
                                      }))
    ## save in single data frame
    chiCordf <- data.frame(expand.grid(zcor = corrs, logkap = log(kseq),
                                       rep = 1:length(precPvals)),
                           chicor = c(chiCorrs))
    ## save this as an RDS
    saveRDS(chiCordf, "chiCorrs.Rds")
} else {
    chiCorrs <- readRDS(chiCordf)
}
## fit a 5th order polynomial
chiCorMods <- lapply(log(kseq), function(k) {
    lm(chicor ~ I(zcor^2) + I(zcor^4) + I(zcor^6) + I(zcor^8) +
           I(zcor^10), data = chiCordf, subset = abs(chiCordf$logkap - k)<0.0001)$coefficients
    })

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
corrImg(popCorrelation(jax_bsbChr), col = pal, breaks = breaks)

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
                                  ncol = nsim))
set.seed(11200905)
for (ii in 1:nsim) {
    for (jj in seq_along(ns)) {
        trait <- beta*colMeans(jax_scores[sample(1:20, ns[jj]),,
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

## plot observed curves against quantiles
ind <- 3
mat <- curves
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
        col = adjustcolor("gray", 0.5), border = NA)

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
ind <- 5
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
corrImg(jax_blockCor, col = pal, breaks = breaks, xaxt = "n",
        yaxt = "n")
addChromosomeLines(jax_bsbBlock)

## the same analysis
ns <- c(1/length(inds), 0.25, 0.5, 0.75, 1)*length(inds)
beta <- 1
sd <- 0.3
jax_scores <- sapply(jax_bsbBlock$encodings, scoreAdditive)
## repeat generation and testing many times
nsim <- 1e3
kseq <- exp(seq(-10, 10, by = 0.5))
curves <- nulls <- array(dim = c(length(ns), length(kseq), nsim))
nullpmat <- pmat <- array(dim = c(length(ns), nrow(jax_scores), nsim))
adjcurves <- array(dim = c(length(ns), length(kseq), nsim))
set.seed(14521105)
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
        padj <- sapply(kseq, poolChiDep, p = ps,
                       sigma = jax_blockCor,
                       method = "satterthwaite")
        nullpl <- sapply(kseq, poolChi, p = null)
        curves[jj, , ii] <- pools
        adjcurves[jj, , ii] <- padj
        nulls[jj, , ii] <- nullpl
        pmat[jj, , ii] <- ps
        nullpmat[jj, , ii] <- null
    }
    if (ii %% 50 == 0) cat("\r Done", ii, "of", nsim)
}

## plot observed curves against quantiles
ind <- 2
mat <- adjcurves
plot(NA, xlim = range(log(kseq, base = 10)),
     ylim = c(-50,0),
            #range(c(log(jax_mgPooled, 10),
            #        log(jax_pgPooled, 10))),
     ylab = expression(paste(log[10], "(p)")),
     xlab = expression(paste(log[10], "(", kappa, ")")))
abline(h = log(0.05, 10), lty = 2, col = "firebrick")
abline(h = seq(0, -50, by = -10), v = seq(-4, 4, by = 2),
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
