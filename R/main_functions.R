#' Prepare CGeneric Data for INLA
#' 
#' @description This function preprocesses data structures used in INLA's CGeneric models.
#' @param data A list containing `ints`, `doubles`, `smatrices`, and `matrices` as required by CGeneric models.
#' @return A list with properly formatted data elements.
#' @export
prepare_cgeneric_data <- function(data) {
  # Preprocess integers
  if (!is.null(data$ints)) {
    data$ints <- lapply(data$ints, function(int_vec) {
      if (length(int_vec) == 0 || is.null(int_vec)) {
        as.integer(NA) # Replace empty or NULL values with NA
      } else {
        as.integer(int_vec)
      }
    })
  } else {
    data$ints <- list() # Ensure ints is at least an empty list
  }

  # Preprocess doubles
  if (!is.null(data$doubles)) {
    data$doubles <- lapply(data$doubles, function(double_vec) {
      if (length(double_vec) == 0 || is.null(double_vec)) {
        as.numeric(NA) # Replace empty or NULL values with NA
      } else {
        as.numeric(double_vec)
      }
    })
  } else {
    data$doubles <- list() # Ensure doubles is at least an empty list
  }

  # Preprocess sparse matrices
  if (!is.null(data$smatrices)) {
    data$smatrices <- lapply(data$smatrices, function(sm) {
      nrow <- sm[1]
      ncol <- sm[2]
      n <- sm[3]
      if (n == 0) {
        return(NULL) # Handle completely empty sparse matrices
      }
      i <- sm[4:(4 + n - 1)]
      j <- sm[(4 + n):(4 + 2 * n - 1)]
      x <- sm[(4 + 2 * n):(4 + 3 * n - 1)]
      list(
        nrow = as.integer(nrow),
        ncol = as.integer(ncol),
        n = as.integer(n),
        i = as.integer(i),
        j = as.integer(j),
        x = as.numeric(x)
      )
    })
  } else {
    data$smatrices <- list() # Ensure smatrices is at least an empty list
  }

  # Preprocess dense matrices
  if (!is.null(data$matrices)) {
    data$matrices <- lapply(data$matrices, function(mat) {
      if (length(mat) < 2) {
        return(NULL) # Handle completely empty matrices
      }
      nrow <- mat[1]
      ncol <- mat[2]
      x <- mat[-(1:2)] # Exclude the first two elements (nrow, ncol)
      list(
        nrow = as.integer(nrow),
        ncol = as.integer(ncol),
        x = as.numeric(x)
      )
    })
  } else {
    data$matrices <- list() # Ensure matrices is at least an empty list
  }

  return(data)
}

#' Call Dynamic CGeneric Model
#' 
#' @description Dynamically calls a CGeneric model's shared library.
#' @param model A model list containing `f$cgeneric` with `shlib`, `model`, and `data`.
#' @param cmd Command string specifying the operation (`"void"`, `"Q"`, etc.).
#' @param theta Numeric vector of parameters for the CGeneric model.
#' @return Result of the dynamic function call.
#' @export
call_dynamic_cgeneric_model <- function(model, cmd, theta) {
  # Validate model structure
  if (!is.list(model) || is.null(model$f) || is.null(model$f$cgeneric)) {
    stop("The model must be a list containing an $f$cgeneric structure.")
  }

  cgeneric <- model$f$cgeneric

  if (!is.list(cgeneric) || is.null(cgeneric$shlib) || is.null(cgeneric$model) || is.null(cgeneric$data)) {
    stop("The model$f$cgeneric structure must contain 'shlib', 'model', and 'data' elements.")
  }

  # Preprocess data
  cgeneric$data <- prepare_cgeneric_data(cgeneric$data)

  # Convert cmd to integer
  cmd <- switch(
    cmd,
    "void" = 0,
    "Q" = 1,
    "graph" = 2,
    "mu" = 3,
    "initial" = 4,
    "log_norm_const" = 5,
    "log_prior" = 6,
    "quit" = 7,
    stop(sprintf("Invalid command '%s'", cmd))
  )

  # Validate theta
  if (!is.numeric(theta) || length(theta) == 0) {
    stop("The 'theta' argument must be a non-empty numeric vector.")
  }

  # Call the C function with PACKAGE specified
    result <- tryCatch(
      {
        .Call(
          "call_dynamic_inla_cgeneric",
          as.integer(cmd),
          as.numeric(theta),
          cgeneric$data,
          cgeneric$model,
          cgeneric$shlib
        )
      },
      error = function(e) {
        stop(
          sprintf(
            "Error calling function '%s' in shared library '%s': %s",
            cgeneric$model,
            cgeneric$shlib,
            e$message
          )
        )
      }
    )

  return(result)
}

#' Build Q Matrix from Sparse Representation
#' 
#' @description Converts a sparse representation of a precision matrix into a full dense matrix.
#' @param Q Vector containing non-zero entries of the precision matrix.
#' @param graph Graph representation (row, column indices, and counts).
#' @return Dense matrix representation of the precision matrix.
#' @export
build_Q_cgeneric <- function(Q, graph) {
  # Extract nrow (and ncol, since it's square) from the first element of graph
  n <- graph[1]
  
  # Extract the number of nonzero elements from the second element of graph
  num_nonzero <- graph[2]
  
  # Extract the row (i) and column (j) indices from the graph
  i <- graph[3:(2 + num_nonzero)]
  j <- graph[(3 + num_nonzero):(2 + 2 * num_nonzero)]
  
  # Extract the nonzero values from Q (skipping the first two elements)
  Q_values <- Q[-c(1, 2)]
  
  # Ensure the lengths match the number of nonzero elements
  if (length(Q_values) != num_nonzero) {
    stop("Mismatch between the number of nonzero entries and the provided values in Q.")
  }
  if (length(i) != num_nonzero || length(j) != num_nonzero) {
    stop("Mismatch between the graph indices and the number of nonzero entries.")
  }
  
  # Create a sparse matrix using the Matrix package
  Q_sparse <- Matrix::sparseMatrix(
    i = i+1,
    j = j+1,
    x = Q_values,
    dims = c(n, n),
    symmetric = TRUE  # Assuming Q is symmetric
  )
  return(Q_sparse)
}

#' @useDynLib inla_cgeneric_debugger, .registration = TRUE
NULL