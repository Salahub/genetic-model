source("effectiveDimFuncs.R")

MeffFuns <- list(cheverud = MeffChev,
                 liji = MeffLiJi,
                 ed = MeffED,
                 galwey = MeffGalwey)

## simple bounds for a sequence of rhos, Ms
rhoSeq <- c(seq(0, 0.95, by = 0.05), 0.99)
Mseq <- c(10, 25, 50, 100, 250, 500, 1000, 2500)
bounds <- lapply(Mseq,
                 function(M) {
                     lapply(rhoSeq, genEigenBounds, M = M)
                 })
## compare these bounds to the circulant approximation
circEigen <- lapply(Mseq,
                    function(M) {
                        lapply(rhoSeq, genCircApprox, M = M)
                    })

## for M = 10, 50, 250 display the results for rho = 0,0.5,0.9,0.99
plotBounds(bounds[[1]][rhoSeq %in% c(0, 0.5, 0.9, 0.99)],
           colSeq = colorRampPalette(c("steelblue3", "black"))(4),
           alpha = 0.3, xlab = "Index", ylab = "Eigenvalue",
           ylim = c(0,20))
legend(x = "topright", legend = c("0", "0.5", "0.9", "0.99"),
       fill = adjustcolor(colorRampPalette(c("steelblue3", "black"))(4),
                          0.3),
       title = expression(rho))

## compare these bounds to the circulant approximation
circDiffs <- lapply(1:length(Mseq),
                    function(iM) {
                        lapply(1:length(rhoSeq),
                               function(ir) {
                                   rbind((bounds[[iM]][[ir]][1,] - circEigen[[iM]][[ir]])/abs(circEigen[[iM]][[ir]]),
                                         (bounds[[iM]][[ir]][2,] - circEigen[[iM]][[ir]])/abs(circEigen[[iM]][[ir]]))
                               })
                    })
## plot these difference bounds
plotBounds(circDiffs[[5]][rhoSeq %in% c(0.5, 0.9, 0.99)],
           colSeq = colorRampPalette(c("steelblue3", "black"))(3),
           alpha = 0.3, xlab = "Index", ylab = "Relative difference",
           ylim = c(-1,1))
abline(h = 0)
legend(x = "topright", legend = c("0.5", "0.9", "0.99"),
       fill = adjustcolor(colorRampPalette(c("steelblue3", "black"))(3), 0.3),
       title = expression(rho))


## increased resolution of bounds for the effective dimension work
rhoSeq <- seq(0, 0.999, by = 0.001)
Mseq <- c(10, 25, 50, 100, 250, 500, 1000, 2500)
bounds <- lapply(Mseq,
                 function(M) {
                     lapply(rhoSeq, genEigenBounds, M = M)
                 })

## get effective dimensions for each M
effectiveDims <- lapply(MeffFuns,
                        function(fun) {
                            lapply(bounds,
                                   function(m) {
                                       sapply(m,
                                              function(r) {
                                                  c(upper = fun(r["upper",]),
                                                    lower = fun(r["lower",]))
                                              })
                                   })
                        })

## look at these
plot(NA, xlim = c(0,1), ylim = c(0, 1000), xlab = expression(rho),
     ylab = "Effective number of tests")
abline(v = seq(0, 1, 0.1), col = "gray70", lty = 2)
abline(h = seq(0, 1000, 100), col = "gray70", lty = 2)
meth <- "galwey" # change for different plots
for (m in Mseq[1:7]) {
    tempData <- effectiveDims[[meth]][[which(Mseq == m)]]
    polygon(c(rhoSeq, rev(rhoSeq)), c(tempData[1,], rev(tempData[2,])),
            col = adjustcolor("steelblue", 0.3), border = "steelblue")
}
