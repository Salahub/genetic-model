# An implementation of a simple genetic model

Several different code files are used here for the gathering of experimental data, the implementation of the model, and the simulation of correlations under the model. **modelScheme.R** is not exceptionally relevant, it simply an image that summarizes the model.

The files **mgiFunctions.R** and **mgiTests.R** are used to extract, process, and display data from the MGI database at http://www.informatics.jax.org/. This is only used tangentially in the work, but could be relied upon more heavily for testing and verification in the future.

The files **simulationFunctions.R** and **simulationRedo.R** contain the implementation of the genetic model and a number of small simulations using the model, respectively.

Finally, **effectiveDim.R** gives the work required to compare different multiple testing approaches using the bounds of Grenander and Szego and the approximate circulant matrix.
