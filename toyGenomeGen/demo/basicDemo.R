## define the marker distribution by chromosome, each element gives
## the count of markers on the corresponding chromosome given by its
## name
markerDist <- c("1" = 5, "2" = 10, "3" = 7, "X" = 5)
## construct a genome based on this distribution
g <- simGenome(markerDist)
g # simple print method
plot(g, add.legend = TRUE) # and plot method
## changing the location arguments gives different structures
g <- simGenome(markerDist, markerFuns = markerPureDom,
               locFuns = locationUniform)
plot(g, add.legend = TRUE)

## by using list arguments, different locations and markers can be
## generated for each chromosome separately
g <- simGenome(markerDist,
               markerFuns = list(markerPureDom, markerPureRec,
                                 markerHybrid, markerHybrid),
               locFuns = list(locationUniform, locationRegular,
                              locationUniform, locationRegular))
plot(g, add.legend = TRUE)

## wrapping the defaults allows the user to change parameters
uni100500 <- function(n) locationUniform(n, min = 100, max = 500)
reg25 <- function(n) locationRegular(n, delta = 25)
g <- simGenome(markerDist,
               markerFuns = list(markerPureDom, markerPureRec,
                                 markerHybrid, markerHybrid),
               locFuns = list(uni100500, reg25,
                              locationUniform, locationRegular))
plot(g, add.legend = TRUE)

## changing the plot arguments controls the appearance
plot(g, add.legend = TRUE, chrLens = c(500, 400, 250, 400),
     epch = c(19, 14, 9, 1))
## under different encodings, epch must be set as well
formals(empiricalGenome:::plot.genome)

## more specific control is given by the makeGenome function, which
## constructs a genome given all of the slots, checking
## consistency at the same time
locs <- list("1" = c(10, 25, 60, 190, 200),
             "2" = c(1, 14, 20, 50, 200, 240, 250, 270, 300, 400),
             "3" = c(90, 100, 110, 150, 190, 240, 300),
             "X" = c(80, 160, 240, 300, 320))
chr <- factor(c(rep("1", 5), rep("2", 10), rep("3", 7), rep("X", 5)))
alls <- rep(list(c("A", "a")), length(chr))
enc <- markerRandom(length(chr))
g <- makeGenome(location = locs, alleles = alls, chromosome = chr,
                encoding = enc)
plot(g, add.legend = TRUE)

## data provided in the correct data.frame format can also be
## converted to a genome, without necessarily being in the form of
## an encoding
alls <- c("A", "a") # simple biallelic case
vals <- c(1, 0) # give encoded values for alleles
mv <- sample(alls, size = 27, replace = TRUE)
pv <- sample(alls, size = 27, replace = TRUE) # example annotations
gDF <- data.frame(mv = mv, pv = pv, chr = chr,
                  pos = unlist(locs)) # same structure as last example
head(gDF)
g <- asGenome(gDF, alleles = alls, values = vals)
plot(g, add.legend = TRUE)

## conformable genomes can be combined using sex
g1 <- simGenome(nmark = markerDist, markerFuns = markerHybrid)
g2 <- simGenome(nmark = markerDist, markerFuns = markerPureDom)
set.seed(31415) # reproducibility
offspring <- sex(g1, g2)
## this assumes the Haldane map/recombination probabilities by
## default, but we can change this behaviour by specifying
## crossFun of map (mapKosambi provided in the function)
set.seed(31415)
offspring2 <- sex(g1, g2, map = mapKosambi)
## we can also provide recombination probabilities for each marker
## directly
p1 <- lapply(markerDist, function(locs) runif(locs - 1))
p2 <- lapply(g2$location, function(ds) 1/diff(ds))
set.seed(31415)
offspring3 <- sex(g1, g2, probs1 = p1, probs2 = p2)
## plot to see the difference that results
plot(offspring, add.legend = TRUE, main = "Default")
dev.new()
plot(offspring2, add.legend = TRUE, main = "Kosambi probabilities")
dev.new()
plot(offspring3, add.legend = TRUE, main = "Custom probabilities")

## making a population is as simple as repeating a given cross many
## times
popsize <- 100
pop1 <- asPopulation(replicate(popsize, sex(g1, g2),
                               simplify = FALSE))
## we can then plot correlation
corObs <- popCorrelation(pop1)
corrImg(corObs, xaxt = "n", yaxt = "n")
addChromosomeLines(pop1)
## compare this to the theoretical correlation
corTh <- theoryCorrelation(g1)
corrImg(corTh, xaxt = "n", yaxt = "n")
addChromosomeLines(g1)
## could change the theoretical settings
corTh <- theoryCorrelation(g1, map = mapKosambi,
                           setting = "backcross")
corrImg(corTh, xaxt = "n", yaxt = "n")
addChromosomeLines(g1)

## some data is included in the function from real experiments
## crossing pure inbred mouse strains
data(package = "empiricalGenome")
data(ucla_bsb) # take the smallest for convenience
x1 <- selectGenome(ucla_bsb, ind = 1)
x1
plot(x1) # we can see where cross-overs occurred by the switch in pch
## compute population correlations
corUcla <- popCorrelation(ucla_bsb,
                          use = "pairwise.complete.obs")
corrImg(corUcla, xaxt = "n", yaxt = "n")
addChromosomeLines(ucla_bsb, lncol = "black")
## theoretical correlations
corUclaTh <- theoryCorrelation(ucla_bsb)
corrImg(corUclaTh, xaxt = "n", yaxt = "n")
addChromosomeLines(ucla_bsb, lncol = "black")
