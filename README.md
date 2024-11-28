
# inla.cgeneric.debugger

**inla.cgeneric.debugger** is an R package designed to assist in debugging and working with INLA (Integrated Nested Laplace Approximation) CGeneric models. It provides utility functions for preparing data structures, calling dynamic CGeneric models, and building precision matrices.

## Installation

You can install the package directly from GitHub:

```R
devtools::install_github("vpnsctl/inla.cgeneric.debugger")
```

## Loading the Package

Load the package into your R session:

```R
library(inla.cgeneric.debugger)
```

## Overview

This package is particularly useful for debugging CGeneric models in INLA, allowing for flexible data preparation and dynamic library integration. It provides the following key functions:
- **call_dynamic_cgeneric_model()**: Dynamically calls a CGeneric model’s shared library with options for verbose output.
- **build_Q_cgeneric()**: Constructs the full precision matrix from a sparse representation.

## Function Details

### call_dynamic_cgeneric_model()

The `call_dynamic_cgeneric_model` function is used to dynamically call a CGeneric model’s shared library. It is especially helpful for debugging because it allows you to see the data passed to the C function and outputs the prints from the CGeneric functions.

#### Usage

```R
result <- call_dynamic_cgeneric_model(
  model,
  cmd,
  theta,
  verbose = FALSE
)
```

#### Arguments
- **model**: A model object returned by the `INLA::inla.cgeneric.define()` function. This object contains the CGeneric model definition, including the shared library path, model name, and data.
- **cmd**: A command string specifying the operation ("void", "Q", "graph", "mu", "initial", "log_norm_const", "log_prior", "quit").
- **theta**: A numeric vector of parameters for the CGeneric model.
- **verbose**: A logical value. If `TRUE`, the function will display the data passed to the C function, which is useful for debugging. If `FALSE`, it will return the command results and show the prints from the CGeneric functions.

#### Features
- **Debugging Aid**: By setting `verbose = TRUE`, you can inspect the data being passed to the C function, making it easier to identify issues.
- **Direct C Function Calls**: The function directly calls the C function without invoking INLA, which makes it faster and more efficient for testing and debugging.
- **Flexible Input**: Accepts a model object from `INLA::inla.cgeneric.define()`, ensuring compatibility with INLA-defined models.

#### Example

```R
# Load required packages
library(fmesher)
library(inlabru)
library(rSPDE)

# Generate a 2D mesh
n_loc <- 20
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(
  loc = loc_2d_mesh,
  cutoff = 0.01,
  max.edge = c(0.1, 0.5)
)

# Define an anisotropic model
model_aniso <- rspde.anistropic2d(mesh = mesh_2d)

# Call the CGeneric model with verbose output to see data passed to C
Q_inla_format <- call_dynamic_cgeneric_model(
  model_aniso,
  cmd = "Q",
  theta = c(-1, -1, 0, 0, 0),
  verbose = TRUE
)

# Call the CGeneric model without verbose output to get function prints
graph <- call_dynamic_cgeneric_model(
  model_aniso,
  cmd = "graph",
  theta = c(-1, -1, 0, 0, 0),
  verbose = FALSE
)
```

### build_Q_cgeneric()

Constructs the full precision matrix from a sparse representation returned by a CGeneric model.

#### Usage

```R
Q_sparse <- build_Q_cgeneric(Q, graph)
```

#### Arguments
- **Q**: Vector containing non-zero entries of the precision matrix (in INLA’s internal format), obtained from `call_dynamic_cgeneric_model()` with `cmd = "Q"`.
- **graph**: Graph representation (row, column indices, and counts), obtained from `call_dynamic_cgeneric_model()` with `cmd = "graph"`.

#### Example

Using the Q and graph computed previously:

```R
# Build the sparse precision matrix from the sparse representation
Q_sparse <- build_Q_cgeneric(Q_inla_format, graph)
```

#### Notes
- The function returns a sparse matrix (`dgCMatrix` class) representing the precision matrix Q.
- `call_dynamic_cgeneric_model()` with `cmd = "Q"` returns Q in INLA’s internal format, which `build_Q_cgeneric()` converts into a standard sparse matrix.

## Example Workflow

Below is an illustrative example of how to use the package in combination with other packages like `fmesher`, `inlabru`, and `rSPDE`.

### 1. Generate a 2D Mesh

```R
library(fmesher)
n_loc <- 20
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(
  loc = loc_2d_mesh,
  cutoff = 0.01,
  max.edge = c(0.1, 0.5)
)
```

### 2. Define an Anisotropic Model

```R
library(inlabru)
library(rSPDE)
model_aniso <- rspde.anistropic2d(mesh = mesh_2d)
```

### 3. Call the CGeneric Model with Verbose Output

```R
library(inla.cgeneric.debugger)

# Call with verbose output to see data passed to C
Q_inla_format <- call_dynamic_cgeneric_model(
  model_aniso,
  cmd = "Q",
  theta = c(-1, -1, 0, 0, 0),
  verbose = TRUE
)
```

### 4. Call the CGeneric Model without Verbose Output

```R
# Call without verbose output to get function prints
graph <- call_dynamic_cgeneric_model(
  model_aniso,
  cmd = "graph",
  theta = c(-1, -1, 0, 0, 0),
  verbose = FALSE
)
```

### 5. Build the Sparse Precision Matrix

```R
# Build the sparse precision matrix from the sparse representation
Q_sparse <- build_Q_cgeneric(Q_inla_format, graph)
```

### 6. Result

The variable `Q_sparse` now contains the sparse precision matrix that can be used for further computations or analyses.

## Additional Information

- **Compatibility**: The package is designed to work seamlessly with INLA and other related packages.
- **Efficiency**: Directly calls the C function, making it faster by avoiding unnecessary overhead.
- **Ease of Debugging**: The `verbose` argument in `call_dynamic_cgeneric_model` provides flexibility in debugging by controlling the level of output.

## Authors
- Alexandre Simas
- David Bolin

## License

This package is licensed under the GPL-3 license.

---

Feel free to explore the package and utilize its functions to enhance your debugging and modeling workflow with INLA CGeneric models.
