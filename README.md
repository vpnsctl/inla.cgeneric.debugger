# inla.cgeneric.debugger

`inla.cgeneric.debugger` is an R package designed to assist in debugging and working with INLA (Integrated Nested Laplace Approximation) CGeneric models. It provides utility functions for preparing data structures, calling dynamic CGeneric models, and building precision matrices.

## Installation

You can install the package directly from GitHub:

```r
devtools::install_github("vpnsctl/inla.cgeneric.debugger")
```

Loading the Package

Load the package into your R session:

```{r}
library(inla.cgeneric.debugger)
```

Example: Using the Package with a 2D Mesh and Anisotropic Model

Below is an illustrative example of how to use the package in combination with other packages like fmesher, inlabru, and rSPDE.

1. Generate a 2D Mesh

We first create a 2D mesh using the fm_mesh_2d function from the fmesher package.

```r
library(fmesher)
n_loc <- 20
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(
    loc = loc_2d_mesh,
    cutoff = 0.01,
    max.edge = c(0.1, 0.5)
)
```

2. Define an Anisotropic Model

We define an anisotropic SPDE model using the rspde.anistropic2d function from the rSPDE package.

```r
library(inlabru)
library(rSPDE)
model_aniso <- rspde.anistropic2d(mesh = mesh_2d)
```

3. Call the CGeneric Model

Using the call_dynamic_cgeneric_model function, we compute the precision matrix (Q) and graph structure for the SPDE model.

```r
Q <- call_dynamic_cgeneric_model(model_aniso, cmd = "Q", theta = c(-1, -1, 0, 0, 0))
graph <- call_dynamic_cgeneric_model(model_aniso, cmd = "graph", theta = c(-1, -1, 0, 0, 0))
```

4. Build the Full Precision Matrix

Finally, we use the build_Q_cgeneric function to construct the full precision matrix from the sparse representation.

```{r}
Q_full <- build_Q_cgeneric(Q, graph)
```

5. Result

The variable `Q_full` now contains the dense precision matrix that can be used for further computations or analyses.

Additional Information

This package is particularly useful for debugging CGeneric models in INLA, allowing for flexible data preparation and dynamic library integration.

Authors

	•	Alexandre Simas
	•	David Bolin