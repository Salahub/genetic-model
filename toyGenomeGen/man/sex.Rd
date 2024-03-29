% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomeCombination.R
\name{sex}
\alias{sex}
\title{Crossing two genomes}
\usage{
sex(
  genome1,
  genome2,
  probs1 = NULL,
  probs2 = NULL,
  map = mapHaldane,
  crossFun = crossIndep
)
}
\arguments{
\item{genome1}{a `genome`}

\item{genome2}{a `genome` with markers observed at the same
locations as `genome1`}

\item{probs1}{vector of probabilities to be passed to `meiose`
for `genome1`}

\item{probs2}{vector of probabilities to be passed to `meiose`
for `genome2`}

\item{map}{function which converts map distances to probabilities}

\item{crossFun}{function dictating crossing over in `meiosis` for
both genomes}
}
\value{
A `genome` with the same marker locations as `genome1`
and `genome2` with encodings recombined according to `crossFun`.
}
\description{
Crossing two genomes
}
\details{
As its name suggests, `sex` simulates sexual reproduction
between two `genome`s. It first calls `meiose` on the objects
independently, passing in the optional arguments to `meiose`.
Next, independent assortment randomly assigns one copy from each
of the recombined gametes from the meiosis calls into a new
`genome` and returns it.
}
\author{
Chris Salahub
}
