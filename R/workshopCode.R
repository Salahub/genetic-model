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
## what does this do to pooled p-values?
poolChiDep <- function(ps, kap, sigma) {
    chis <- qchisq(ps, df = kap, lower.tail = FALSE)
    gammaApprox(chis, sigma, kap)
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

## define vector of counts
ns <- c(1, 5, 10, 15, 20)
beta <- 1
sd <- 0.3
## compute scores
jax_scores <- sapply(jax_bsbChr$encodings, scoreAdditive)
## repeat generation and testing many times
nsim <- 1e3
kseq <- exp(seq(-10, 10, by = 0.5))
curves <- array(dim = c(length(ns), length(kseq), nsim))
nulls <- array(dim = c(length(ns), length(kseq), nsim))
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
curves <- array(dim = c(length(ns), length(kseq), nsim))
nulls <- array(dim = c(length(ns), length(kseq), nsim))
set.seed(11200905)
for (ii in 1:nsim) {
    for (jj in seq_along(ns)) {
        trait <- beta*colMeans(jax_scores[sample(1:nrow(jax_scores),
                                                 ns[jj]),,
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
    }
    if (ii %% 50 == 0) cat("\r Done", ii, "of", nsim)
}

## plot observed curves against quantiles
ind <- 1
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
        trait <- beta*colMeans(jax_scores[sample(1:nrow(jax_scores),
                                                 ns[jj]),,
                                          drop = FALSE]) +
            rnorm(94, sd = sd)
        ps <- apply(jax_scores, 1,
                    function(g) kruskal.test(x = trait, g)$p.value)
        null <- apply(jax_scores, 1,
                      function(g) kruskal.test(x = sample(trait),
                                               g)$p.value)
        pools <- sapply(kseq, poolChiDep, p = ps,
                        sigma = jax_bsbSingCor)
        depcurves[jj, , ii] <- pools
    }
    if (ii %% 50 == 0) cat("\r Done", ii, "of", nsim)
}

## plot observed curves against quantiles
ind <- 5
mat <- depcurves
plot(NA, xlim = range(log(kseq, base = 10)),
     ylim = c(-20,0),
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



## add dependence in a block structure
## get 10 markers over the first two chromosomes
k <- 10
chrInds <- cumsum(table(jax_bsb$chromosome))[1:2]
steps <- diff(c(0,chrInds))/k
inds <- floor(rep(c(0, chrInds[1:(length(chrInds)-1)]), each = k) +
              (1:k)*rep(steps, each = k))
## extract markers from jax_bsb
jax_bsbBlock <- subsetPopulation(jax_bsb, markInd = inds)
corrImg(popCorrelation(jax_bsbBlock), col = pal, breaks = breaks,
        xaxt = "n", yaxt = "n")
addChromosomeLines(jax_bsbBlock)

## the same analysis
ns <- c(1/length(inds), 0.25, 0.5, 0.75, 1)*length(inds)
beta <- 1
sd <- 0.3
jax_scores <- sapply(jax_bsbBlock$encodings, scoreAdditive)
## repeat generation and testing many times
nsim <- 1e3
kseq <- exp(seq(-10, 10, by = 0.5))
curves <- array(dim = c(length(ns), length(kseq), nsim))
nulls <- array(dim = c(length(ns), length(kseq), nsim))
set.seed(11200905)
for (ii in 1:nsim) {
    for (jj in seq_along(ns)) {
        trait <- beta*colMeans(jax_scores[sample(1:nrow(jax_scores),
                                                 ns[jj]),,
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
    }
    if (ii %% 50 == 0) cat("\r Done", ii, "of", nsim)
}

## plot observed curves against quantiles
ind <- 2
mat <- curves
plot(NA, xlim = range(log(kseq, base = 10)),
     ylim = c(-50,0),
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

