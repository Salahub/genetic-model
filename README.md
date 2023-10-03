# Exploring a simple genetic model

This repository presents and explores a simple genetic model.

## toyGenomeGenR

First, this model is implemented using the S3 class system of R in the
package contained in the toyGenomeGenR directory. Numerous functions to
generate, cross, subset, and otherwise manipulate genomes are all
provided in this package, with demos to guide users. In addition, this
package provides processed versions of almost a dozen data sets
recording genetic measurements in backcross experiments in mice.

## R

The package is applied in the numerous scripts contained in the other
directory, and several supporting scripts to extract and process the
package data are provided as well.
For example, scripts **mgiFunctions.R** and **mgiTests.R** are used to
extract, process, and display data from the MGI database at
http://www.informatics.jax.org/, the source of the package data sets.
**simulationRedo.R** contains an implementation of the genetic model and a number of small simulations using the model, respectively.
**modelScheme.R** is rather irrelevant, it's simply an image that summarizes the model.
