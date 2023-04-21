## define the marker distribution by chromosome, each element gives
## the count of markers on the corresponding chromosome given by its
## name
markerDist <- c("1" = 5, "2" = 10, "3" = 7, "X" = 8)
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
formals(simpleGenome:::plot.genome)
