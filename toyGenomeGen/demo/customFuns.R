## the empiricalGenomes package is designed to allow users to easily
## define their own functions once they understand the machinery, this
## demo is meant to guide this exercise
## the core of recombination is the meiose function
meiose
## crossIndep is the default recombination function
crossIndep
## we might instead use a renewal process, say with scaled chi^2
## holding times to account for distances in centiMorgans
crossChi <- function(probs, locs) { # must accept two arguments
    locLen <- length(locs) # largest location
    crosses <- numeric() # cross indices
    currPos <- 100*rchisq(1, df = 2) # current position of process
    while (currPos < locs[1]) { # only care within marker range
        new <- 100*rchisq(1, df = 2) # next holding time
        currPos <- currPos + new # update break spot
    }
    while (currPos <= locs[locLen]) { # breaks within markers
        crosses <- c(crosses, sum(locs < currPos)) # index of split
        new <- 100*rchisq(1, df = 2)
        currPos <- currPos + new
    }
    crosses
}
## let's try it out

## define a simple genome
g <- simGenome(nmark = c("1" = 10, "2" = 10))
plot(g)
## look at the default behaviour
popsize <- 1000
defPop <- asPopulation(replicate(popsize, sex(g, g),
                                 simplify = FALSE))
## try this new custom function
chiPop <- asPopulation(replicate(popsize,
                                 sex(g, g, crossFun = crossChi),
                                 simplify = FALSE))

## plot a comparison
## a suggested divergent palette to use when comparing correlation
pal <- colorRampPalette(c("steelblue", "floralwhite", "firebrick"))(49)
corBrks <- seq(-1, 1, length.out = 50)
oldPar <- par(mfrow = c(1,2), mar = c(2.1, 2.1, 3.1, 0.1))
corrImg(popCorrelation(defPop), xaxt = "n", yaxt = "n",
        main = "Haldane map", breaks = corBrks, col = pal)
addChromosomeLines(defPop)
corrImg(popCorrelation(chiPop), xaxt = "n", yaxt = "n",
        main = "Chi renewal process", breaks = corBrks, col = pal)
addChromosomeLines(chiPop)

## we could also adapt the crossChi function to pick a reference point
## first, assuming there is a centromere where the genes come together
## that is likely to cause crosses
crossChi2 <- function(probs, locs) {
    locLen <- length(locs)
    ref <- rnorm(1, locs[locLen]/2, 5) # centromere random near mid
    locs <- locs - ref # convert to pos/negative
    crosses <- numeric()
    ## handle the negative side first
    currPos <- -100*rchisq(1, 2)
    while (currPos > locs[1]) { # above first marker
        crosses <- c(sum(locs < currPos), crosses) # keep order
        new <- 100*rchisq(1, df = 2)
        currPos <- currPos - new # subtract instead of add
    }
    ## now the positive side
    currPos <- 100*rchisq(1, 2)
    while (currPos < locs[locLen]) { # below last marker
        crosses <- c(crosses, sum(locs < currPos)) # keep order
        new <- 100*rchisq(1, df = 2)
        currPos <- currPos + new
    }
    crosses
}
## do the same experiment
chiPop2 <- asPopulation(replicate(popsize,
                                  sex(g, g, crossFun = crossChi2),
                                  simplify = FALSE))
## plot a comparison
par(mfrow = c(1,2), mar = c(2.1, 2.1, 3.1, 0.1))
corrImg(popCorrelation(chiPop), xaxt = "n", yaxt = "n",
        main = "Chi without centromere", breaks = corBrks, col = pal)
addChromosomeLines(chiPop)
corrImg(popCorrelation(chiPop2), xaxt = "n", yaxt = "n",
        main = "Chi with centromere", breaks = corBrks, col = pal)
addChromosomeLines(chiPop2)
## a wide range of different behaviours can therefore be captured
## with custom functions, changing the combination probabilities
## drastically
par(oldPar)
