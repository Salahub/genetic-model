% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomes.R
\name{simGenome}
\alias{simGenome}
\title{Simulate a genome object}
\usage{
simGenome(nmark, alleles, markerFuns = markerHybrid, locFuns = locationRegular)
}
\arguments{
\item{nmark}{vector of integers giving the count of measured
markers by chromosome, if it is named it is assumed the names
correspond to the chromosomes, otherwise chromosomes are numbered
in order}

\item{alleles}{vector of possible allele annotations, if missing
this is taken to be `c("A","a")` for every marker site}

\item{locFuns}{list of functions that generate locations for
markers by chromosome, recycled if it is not as long as nmark}

\item{markerFun}{list of functions that simulate marker
encodings by chromosome, recycled if it is not as long as nmark}
}
\value{
A genome object with encodings and locations generated
by the corresponding functions for each chromosome.
}
\description{
`simGenome` simulates a genome object based on the
distribution of markers in `nmark` and the marker and location
generating functions in `markerFuns` and `locFuns`
}
\details{
`alleles` must be a list of the same length as
`sum(nmark)`, with each element providing the possible alleles for
the corresponding marker. If it is not provided, the Mendelian
A/a biallelic annotation is assumed for every marker. If either of
`markerFuns` or `locFuns` provides a function rather than a list
of functions, this function is recycled across every chromosome.
Otherwise, the length of the list must match the number of
chromosomes to provide the different functions used to generate
encodings and locations on each chromosome.
}
\examples{
simGenome(c("1" = 10, "2" = 5, "X" = 19))
}
\author{
Chris Salahub
}
