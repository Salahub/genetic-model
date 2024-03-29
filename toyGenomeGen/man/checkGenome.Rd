% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomes.R
\name{checkGenome}
\alias{checkGenome}
\title{Checking if an object is a valid genome}
\usage{
checkGenome(genome)
}
\arguments{
\item{genome}{object to be checked}
}
\value{
TRUE if the object is a genome, otherwise text outlining
what check failed.
}
\description{
Checking if an object is a valid genome
}
\details{
Performs checks of slot data types, dimensions, and
conformity to ensure `genome` is a valid instance of the S3 class
`genome`.
}
\author{
Chris Salahub
}
