# Exploring a simple genetic model

This repository presents and explores a simple genetic model.

## toyGenomeGen

First, this model is implemented using the S3 class system of R in the
package contained in the toyGenomeGenR directory. Numerous functions to
generate, cross, subset, and otherwise manipulate genomes are all
provided in this package, with demos to guide users. In addition, this
package provides processed versions of almost a dozen data sets
recording genetic measurements in backcross experiments in mice. The
package can be downloaded via the command
`install_github("Salahub/genetic-model", subdir="toyGenomeGenR")`
using `install_github` from the `devtools` package in R.

## Scripts

The package is applied in the numerous scripts contained in the other
directory, and several supporting scripts to extract and process the
package data are provided as well. Described broadly, these scripts
are

- **mgiConvert.R**: use API calls to https://phenome.jax.org/ and
      data extracted from https://www.informatics.jax.org/ to pull
	  phenotype and genotype data from the Jackson Lab databases and
      convert them to population and genome objects
- **mgiTests.R**: perform simulations based on the Jackson Lab
	  data, corresponds with Section 5.6 of the thesis and contains
	  some additional investigations not explored in the thesis
- **modelScheme.R**: draws the model diagram in Figure 5.1
- **literatureSettings.R**: uses toyGenomeGenR to replicate
	  settings discussed in the literature of effective dimension
- **phenoTypeAnalysis.R**: the script used to complete the example
	  at the end of Chapter 6 of the thesis
- **adjustingPooledChi.R**: the script showing the early
	  exploration and fitting of Chapter 6

The data files either store output or are reference values for these
scripts.

## effectiveDim

This directory contains some short code files used to investigate the
computation of effective dimension in genetic data. This was covered
only briefly in my thesis, and generally seems like a bad idea.
