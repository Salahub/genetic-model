## simple plotting of effective dimension
genEigenBounds <- function(M, rho) {
    lower <- (1 - rho^2)/(1 - 2*rho*cos((1:M)*pi/(M+1)) + rho^2)
    upper <- (1 - rho^2)/(1 - 2*rho*cos((1:M - 1)*pi/(M+1)) + rho^2)
    if (rho == 1) upper[1] <- lower[1] <- M
    rbind(upper, lower)
}

## and the circulant eigenvalue approximation
genCircApprox <- function(M, rho) {
    expSeq <- 0:(M-1)
    coefs <- rho^expSeq + (expSeq/M)*(rho^(M - expSeq) - rho^expSeq)
    cosines <- outer(expSeq, 0:(M-1), function(x,y) cos(2*x*y*pi/M))
    sort(colSums(cosines*coefs), decreasing = TRUE)
}

## exact case for equidistant markers on a single chromosome from
## Narayan and Shastry
qseq <- function(M, kappa) { # helper to get roots
    qFn <- function(x, nu, N = M) { # root function of N. & S.
        (N - 1)*x + 2*atan2(sin(x), cos(x) - exp(-kappa)) - nu*pi
    }
    sapply(1:M, function(ii) uniroot(qFn, # search intervals
                                     interval = c(ii-1, ii)*pi/(M+1),
                                     nu = ii)$root)
}
exactc1edEig <- function(M, rho = exp(-kappa), kappa = -log(rho)) {
    qs <- qseq(M, kappa = kappa) # the zeros
    sinh(kappa)/(cosh(kappa) - cos(qs)) # the eigenvalues
}

## cheeky version that uses the asymptotic limit instead of
## solving N equations
qseqAsym <- function(M, kappa) {
    (1:M)/(M+1)*pi
}
asymc1edEig <- function(M, rho = exp(-kappa), kappa = -log(rho)) {
    qs <- qseqAsym(M, kappa = kappa) # the asymptotic zeros
    sinh(kappa)/(cosh(kappa) - cos(qs))
}

## M eff functions
MeffChev <- function(lambdas) length(lambdas) + 1 - mean(lambdas^2)

MeffLiJi <- function(lambdas) sum((lambdas >= 1) + (lambdas %% 1))

MeffED <- function(lambdas, p = 1/length(lambdas)) {
    sum((lambdas/max(lambdas))^p)
}

MeffGalwey <- function(lambdas) {
    (MeffED(lambdas, 0.5))^2/MeffED(lambdas, p = 1)
}

MeffFuns <- list(cheverud = MeffChev,
                 liji = MeffLiJi,
                 ed = MeffED,
                 galwey = MeffGalwey)

## first inspect these bounds
plotBounds <- function(bnds, colSeq, alpha = 0.5, ylim, ...) {
    neig <- ncol(bnds[[1]])
    xs <- 1:neig
    if (missing(ylim)) ylim <- range(unlist(bnds))
    plot(NA, xlim = c(1,neig), ylim = ylim, ...)
    for (ii in seq_along(bnds)) {
        polygon(c(xs, rev(xs)), c(bnds[[ii]][1,], rev(bnds[[ii]][2,])),
                col = adjustcolor(colSeq[ii], alpha.f = alpha),
                border = colSeq[ii])
        #lines(xs, abs(bnds[[ii]][1,]-bnds[[ii]][2,]),
        #      col = colSeq)
    }
}
