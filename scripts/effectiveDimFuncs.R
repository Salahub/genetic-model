## Chan's circulant eigenvalue approximation
chanCircEig <- function(M, rho) {
    expSeq <- 0:(M-1)
    coefs <- rho^expSeq + (expSeq/M)*(rho^(M - expSeq) - rho^expSeq)
    cosines <- outer(expSeq, 0:(M-1), function(x,y) cos(2*x*y*pi/M))
    sort(colSums(cosines*coefs), decreasing = TRUE)
}
## Gray's circulant eigenvalues approximation
grayCircEig <- function(M, rho) {
    expSeq <- 0:(M-1)
    coefs <- 1/(1 - rho^M)*(rho^expSeq + rho^(M - expSeq))
    cosines <- outer(expSeq, 0:(M-1), function(x,y) cos(2*x*y*pi/M))
    inter <- coefs*cosines
    inter[1,] <- 1
    sort(colSums(inter), decreasing = TRUE)
}

## Stroeker's approximation of eigenvalues
stroekerEig <- function(M, rho) {
    if (abs(rho) < 0.5) {
        angles <- 1:M*pi/(M+1)
        inv <- 1 - 2*rho*cos(angles) + rho^2 - 4/(M+1)*rho^2*sin(angles)^2
    } else if (rho > 0) {
        angles <- 1:(M-1)*pi/M
        inv <- c((1 - rho)^2 + 2/M*rho*(1 - rho),
                 1 - 2*rho*cos(angles) + rho^2 +
                 2/M*rho*(1-rho)*(1 + cos(angles)))
    } else {
        angles <- 1:(M-1)*pi/M
        inv <- c(1 - 2*rho*cos(angles) + rho^2 -
                 2/M*rho*(1+rho)*(1 - cos(angles)),
                 (1 + rho)^2 - 2/M*rho*(1 + rho))
    }
    (1 - rho^2)/inv
}

## exact case for equidistant markers on a single chromosome from
## Narayan and Shastry
qseq <- function(M, kappa) { # helper to get roots
    qFn <- function(x, nu, N = M) { # root function of N. & S.
        (N - 1)*x + 2*atan2(sin(x), cos(x) - exp(-kappa)) - nu*pi
    }
    sapply(1:M, function(ii) uniroot(qFn, # search intervals
                                     interval = c(ii-1, ii)*pi/(M+1),
                                     nu = ii,
                                     tol = .Machine$double.eps^0.75)$root)
}
exactc1edEig <- function(M, rho = exp(-kappa), kappa = -log(rho)) {
    qs <- qseq(M, kappa = kappa) # the zeros
    sinh(kappa)/(cosh(kappa) - cos(qs)) # the eigenvalues
}

## version that uses the asymptotic limit and one NR iteration instead
qseqAsym <- function(M, kappa) {
    xs <- (1:M)/(M+1)*pi # limit value
    rho <- exp(-kappa)
    adjNum <- 2*atan2(sin(xs), cos(xs) - rho) - 2*xs
    adjDen <- M - 1 + 2*(1 - rho*cos(xs))/(1 - 2*rho*cos(xs) + rho^2)
    xs - adjNum/adjDen
}
asymc1edEig <- function(M, rho = exp(-kappa), kappa = -log(rho)) {
    qs <- qseqAsym(M, kappa = kappa) # the asymptotic zeros
    sinh(kappa)/(cosh(kappa) - cos(qs))
}

## add a function for the asmptotic density of eigenvalues
asymc1edDensity <- function(x, M, rho = exp(-kappa), kappa = -log(rho)) {
    (M + x)/(M*pi*x^2)*
        sinh(kappa)/sqrt(1 - (cosh(kappa) - sinh(kappa)/x)^2)
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

MeffGaoEtAl <- function(lambdas, cutoff = 0.995) {
    sum(cumsum(lambdas/sum(lambdas)) <= cutoff)
}

MeffFuns <- list(cheverud = MeffChev,
                 liji = MeffLiJi,
                 ed = MeffED,
                 galwey = MeffGalwey,
                 gao = MeffGaoEtAl)

## simple plotting of effective dimension
genEigenBounds <- function(M, rho) {
    lower <- (1 - rho^2)/(1 - 2*rho*cos((1:M)*pi/(M+1)) + rho^2)
    upper <- (1 - rho^2)/(1 - 2*rho*cos((1:M - 1)*pi/(M+1)) + rho^2)
    if (rho == 1) upper[1] <- lower[1] <- M
    rbind(upper, lower)
}

## inspect these bounds
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
