source("effectiveDimFuncs.R")
source("simulationFunctions.R")

## global settings
nsim <- 100
npop <- 1e2
bluWhRed <- colorRampPalette(c("steelblue", "white", "firebrick"))
effDimFuncs <- list(cheverud = MeffChev, liji =  MeffLiJi,
                    ed = MeffED, galwey = MeffGalwey)

## start incredibly simple: 3 equidistant markers placed along a single
## chromosome 20 cM apart each
## Theory
c1ed <- list(nchrom = 1, nmark = 4, dists = list(c(rep(15, 3)))) # setting
c1edTheory <- theoryCor(c1ed$dists) # theoretical intercross corr
c1edEig <- eigen(c1edTheory) # eigenvalues
## different effective tests
c1edEffDim <- sapply(effDimFuncs, do.call,
                     args = list(lambdas = c1edEig$values))
## Simulated
c1edPar <- list(genome1 = do.call(abiogenesis,
                                 args = c(c1ed, allele = 1)),
                genome2 = do.call(abiogenesis,
                                 args = c(c1ed, allele = 0)))
c1edF1 <- do.call(sex, args = c1edPar)
c1edCors <- replicate(nsim,
                      popCorrelation(replicate(npop,
                                               sex(c1edF1, c1edF1),
                                               simplify = FALSE)))
c1edEffDimSim <- apply(c1edCors, 3,
                       function(cors) {
                           sapply(effDimFuncs, do.call,
                                  args = list(lambdas = eigen(cors)$values))
                       })
## plot all this together
boxplot(t(c1edEffDimSim))
points(x = 1:4, y = c1edEffDim, pch = 19, col = "steelblue")
