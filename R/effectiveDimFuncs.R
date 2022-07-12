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
exactc1edEig <- function(M, rho = exp(-kappa), kappa = -log(rho)) {
    qseq <- (1:M)/(M + 1)*pi
    sinh(kappa)/(cosh(kappa) - cos(qseq))
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
