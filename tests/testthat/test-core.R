# tests/testthat/test-core.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("Package loads without error", {
  library(topolow)
  expect_true(TRUE)
}
)

# Helper function to create test matrices
create_test_matrix <- function(size = 4, with_threshold = FALSE) {
  if (with_threshold) {
    matrix(c(0, ">2", 3, ">2", 0, 4, 3, 4, 0), nrow = 3)
  } else {
    mat <- as.matrix(dist(matrix(runif(size * 2), ncol = 2)))
    colnames(mat) <- rownames(mat) <- paste0("Point", 1:size)
    mat
  }
}

test_that("euclidean_embedding handles input validation correctly", {
  # Basic matrix validation
  expect_error(
    euclidean_embedding(dissimilarity_matrix = "not a matrix", ndim = 2, mapping_max_iter = 10,
                     k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01, write_positions_to_csv = FALSE),
    "dissimilarity_matrix must be a matrix"
  )

  # Non-square matrix
  expect_error(
    euclidean_embedding(dissimilarity_matrix = matrix(1:6, nrow = 2), ndim = 2, mapping_max_iter = 10,
                     k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01, write_positions_to_csv = FALSE),
    "dissimilarity_matrix must be square"
  )

  # Parameter validation
  params_to_test <- list(
    list(ndim = -1, msg = "ndim must be a positive integer"),
    list(k0 = -1, msg = "k0 must be a positive number"),
    list(cooling_rate = 1.5, msg = "cooling_rate must be between 0 and 1"),
    list(c_repulsion = 0, msg = "c_repulsion must be a positive number"),
    list(relative_epsilon = -1, msg = "relative_epsilon must be a positive number"),
    list(convergence_counter = 0.5, msg = "convergence_counter must be a positive integer")
  )

  test_mat <- create_test_matrix()
  for (param in params_to_test) {
    args <- list(
      dissimilarity_matrix = test_mat,
      ndim = 2, mapping_max_iter = 10,
      k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
      relative_epsilon = 1e-4, convergence_counter = 5
    )
    args[[names(param)[1]]] <- param[[1]]
    expect_error(do.call(euclidean_embedding, args), param$msg)
  }

  # Test high k0 warning
  expect_warning(
    euclidean_embedding(dissimilarity_matrix = test_mat, ndim = 2, mapping_max_iter = 10,
                     k0 = 35, cooling_rate = 0.01, c_repulsion = 0.01, write_positions_to_csv = FALSE),
    "High k0 value"
  )
})

test_that("euclidean_embedding handles initial positions correctly", {
  test_mat <- create_test_matrix()
  n <- nrow(test_mat)
  ndim <- 2

  # Test valid initial positions
  init_pos <- matrix(runif(n * ndim), nrow = n, ncol = ndim)
  rownames(init_pos) <- rownames(test_mat)
  result <- euclidean_embedding(test_mat, ndim = ndim, mapping_max_iter = 10,
                             k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
                             initial_positions = init_pos, write_positions_to_csv = FALSE)
  expect_equal(dim(result$positions), c(n, ndim))

  # Test invalid dimensions
  wrong_pos <- matrix(runif((n + 1) * ndim), nrow = n + 1, ncol = ndim)
  expect_error(
    euclidean_embedding(test_mat, ndim = ndim, mapping_max_iter = 10,
                     k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
                     initial_positions = wrong_pos, write_positions_to_csv = FALSE),
    "initial_positions must have same number of rows"
  )
})

test_that("euclidean_embedding handles missing values and thresholds", {
  test_mat <- create_test_matrix(with_threshold = TRUE)
  test_mat[1, 3] <- test_mat[3, 1] <- NA # Add a missing value

  result <- euclidean_embedding(test_mat, ndim = 2, mapping_max_iter = 10,
                             k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
                             write_positions_to_csv = FALSE)

  # Check that NA values result in a calculated distance (repulsive force)
  expect_true(!is.na(result$est_distances[1, 3]))
  expect_equal(result$est_distances[1, 3], result$est_distances[3, 1]) # Symmetry

  # Check that threshold values are handled
  expect_true(is.numeric(result$est_distances))
})

test_that("euclidean_embedding force calculations preserve distance relationships", {
  # Create triangle with known properties
  test_mat <- matrix(c(0,1,2, 1,0,1, 2,1,0), nrow=3)
  colnames(test_mat) <- rownames(test_mat) <- c("A","B","C")
  
  result <- euclidean_embedding(dissimilarity_matrix = test_mat, ndim=2, mapping_max_iter=10,
                        k0=1.0, cooling_rate=0.01, c_repulsion=0.01, 
                        write_positions_to_csv=FALSE)
  
  # Get pairwise distances from result
  get_dist <- function(p1, p2) {
    sqrt(sum((result$positions[p1,] - result$positions[p2,])^2))
  }
  
  # Check triangle inequality and relative distances
  d_ab <- get_dist("A","B")
  d_bc <- get_dist("B","C")
  d_ac <- get_dist("A","C")
  
  expect_true(d_ac > d_ab)  # AC should be longest
  expect_true(d_ac < (d_ab + d_bc))  # Triangle inequality
})

test_that("euclidean_embedding returns correct object structure", {
  result <- euclidean_embedding(create_test_matrix(), ndim = 2, mapping_max_iter = 10,
                             k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01, write_positions_to_csv = FALSE)

  expected_elements <- c("positions", "est_distances", "mae", "iter",
                         "parameters", "convergence")
  expect_true(all(expected_elements %in% names(result)))
  expect_s3_class(result, "topolow")
  expect_true(is.numeric(result$mae))
  expect_true(is.logical(result$convergence$achieved))
})


###########################################################################
################## Testing the wizard function Euclidify ##################


# Helper function to create test dissimilarity matrices
create_test_dissimilarity_matrix <- function(size = 20, with_missing = FALSE, with_thresholds = FALSE) {
  points <- paste0("Point", 1:size)
  
  test_data <- expand.grid(object = points, reference = points)
  
  if (with_thresholds) {
    test_data$score <- sample(c(1, 2, 4, 8, 16, 32, 64, "<1", ">64"), 
                              nrow(test_data), replace = TRUE)
  } else {
    test_data$score <- sample(c(1, 2, 4, 8, 16, 32, 64), 
                              nrow(test_data), replace = TRUE)
  }
  
  # Convert to matrix
  dist_mat <- list_to_matrix(
    data = test_data,
    object_col = "object",
    reference_col = "reference",
    value_col = "score",
    is_similarity = TRUE
  )
  
  # Add missing values if requested
  if (with_missing) {
    n_missing <- round(0.2 * size^2)  # 20% missing
    missing_indices <- sample(size^2, n_missing)
    dist_mat[missing_indices] <- NA
  }
  
  return(dist_mat)
}


test_that("Euclidify input validation works correctly", {
  test_mat <- create_test_dissimilarity_matrix(5)
  temp_dir <- tempdir()
  
  # Test missing output_dir
  expect_error(
    Euclidify(dissimilarity_matrix = test_mat, max_cores = 1),
    "output_dir must be specified"
  )
  
  # Test non-matrix input
  expect_error(
    Euclidify(dissimilarity_matrix = "not_a_matrix", output_dir = temp_dir,max_cores = 1),
    "dissimilarity_matrix must be a matrix"
  )
  
  # Test non-square matrix
  non_square <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(
    Euclidify(dissimilarity_matrix = non_square, output_dir = temp_dir, max_cores = 1),
    "dissimilarity_matrix must be square"
  )
  
  # Test invalid ndim_range
  expect_error(
    Euclidify(dissimilarity_matrix = test_mat, output_dir = temp_dir, 
              ndim_range = c(5, 3),
              max_cores = 1),  # min > max
    "ndim_range must be a vector of length 2 with min <= max"
  )
  
  expect_error(
    Euclidify(dissimilarity_matrix = test_mat, output_dir = temp_dir, 
              max_cores = 1,
              ndim_range = c(2)),  # wrong length
    "ndim_range must be a vector of length 2 with min <= max"
  )
  
  # Test invalid verbose parameter
  expect_error(
    Euclidify(dissimilarity_matrix = test_mat, output_dir = temp_dir, 
              max_cores = 1,
              verbose = "invalid"),
    "verbose must be 'off', 'standard', or 'full'"
  )
})


test_that("Euclidify works with basic valid input", {
  test_mat <- matrix(c(0, 1, 2, 3, 1, 0, 2.5, 3.5, 2, 2.5, 0, 4, 3, 3.5, 4, 0), 4, 4)
  rownames(test_mat) <- colnames(test_mat) <- paste0("Point", 1:4)
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Test with minimal parameters for speed
  result <- Euclidify(
    dissimilarity_matrix = test_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 6),
    n_initial_samples = 20,      # Increased slightly for stability
    n_adaptive_samples = 20,    # Increased slightly for stability
    max_cores = 1,
    folds = 10,
    mapping_max_iter = 100,      # Low for speed
    verbose = "off",
    fallback_to_defaults = TRUE
  )
  
  # Check return structure
  expected_elements <- c("positions", "est_distances", "mae", "optimal_params", 
                         "optimization_summary", "data_characteristics", "runtime")
  expect_true(all(expected_elements %in% names(result)))
  
  # Check positions matrix
  expect_true(is.matrix(result$positions))
  expect_equal(nrow(result$positions), nrow(test_mat))
  expect_true(ncol(result$positions) >= 2)
  expect_true(ncol(result$positions) <= 6)
  
  # Check estimated distances
  expect_true(is.matrix(result$est_distances))
  expect_equal(dim(result$est_distances), dim(test_mat))
  
  # Check optimal parameters
  expect_true(is.list(result$optimal_params))
  param_names <- c("ndim", "k0", "cooling_rate", "c_repulsion")
  expect_true(all(param_names %in% names(result$optimal_params)))
  
  # Check data characteristics
  expect_true(is.list(result$data_characteristics))
  expect_true("n_objects" %in% names(result$data_characteristics))
  expect_equal(result$data_characteristics$n_objects, nrow(test_mat))
  
  # Check runtime
  expect_true(result$runtime > 0)
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})


test_that("Euclidify handles missing data correctly", {
  test_mat <- create_test_dissimilarity_matrix(10, with_missing = TRUE)
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  result <- Euclidify(
    dissimilarity_matrix = test_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 6),
    n_initial_samples = 8,
    n_adaptive_samples = 12,
    max_cores = 1,
    folds = 5,
    mapping_max_iter = 50,
    verbose = "off",
    fallback_to_defaults = TRUE
  )
  
  # Should complete successfully even with missing data
  expect_true(is.matrix(result$positions))
  expect_true(result$data_characteristics$missing_percentage > 0)
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})


test_that("Euclidify handles threshold indicators correctly", {
  test_mat <- create_test_dissimilarity_matrix(60, with_thresholds = TRUE)
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  result <- Euclidify(
    dissimilarity_matrix = test_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 4),
    n_initial_samples = 15,
    n_adaptive_samples = 15,
    max_cores = 1,
    folds = 5,
    mapping_max_iter = 5,
    verbose = "off",
    fallback_to_defaults = TRUE
  )
    
  # Should complete successfully with threshold indicators
  expect_true(is.matrix(result$positions))
  expect_true(is.numeric(result$mae))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})


test_that("Euclidify verbose modes work correctly", {
  test_mat <- create_test_dissimilarity_matrix(10)
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Test "off" mode (should not error)
  expect_no_error({
    result_off <- Euclidify(
      dissimilarity_matrix = test_mat,
      output_dir = temp_dir,
      ndim_range = c(2, 2),
      n_initial_samples = 2,
      n_adaptive_samples = 2,
      max_cores = 1,
      folds = 5,
      mapping_max_iter = 20,
      verbose = "off",
      fallback_to_defaults = TRUE
    )
  })
  
  # Test "standard" mode (should not error)
  expect_no_error({
    result_standard <- Euclidify(
      dissimilarity_matrix = test_mat,
      output_dir = file.path(temp_dir, "standard"),
      ndim_range = c(2, 2),
      n_initial_samples = 2,
      n_adaptive_samples = 2,
      max_cores = 1,
      folds = 5,
      mapping_max_iter = 20,
      verbose = "standard",
      fallback_to_defaults = TRUE
    )
  })
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})


test_that("Euclidify fallback mechanism works", {
  # Create a very small matrix that might cause optimization issues
  small_mat <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), nrow = 3)
  rownames(small_mat) <- colnames(small_mat) <- paste0("P", 1:3)
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Test with fallback enabled (default)
  result <- Euclidify(
    dissimilarity_matrix = small_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 2),
    n_initial_samples = 4,    # Keep small to potentially trigger fallback
    n_adaptive_samples = 4,
    max_cores = 1,
    folds = 3,                # Reduced folds for small matrix
    mapping_max_iter = 20,
    verbose = "off",
    fallback_to_defaults = TRUE
  )
  
  # Should still return a valid result
  expect_true(is.matrix(result$positions))
  expect_true(is.list(result$optimal_params))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})


test_that("Euclidify save_results parameter works", {
  test_mat <- create_test_dissimilarity_matrix(10)
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  result <- Euclidify(
    dissimilarity_matrix = test_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 2),
    n_initial_samples = 2,
    n_adaptive_samples = 2,
    max_cores = 1,
    folds = 5,
    mapping_max_iter = 20,
    verbose = "off",
    save_results = TRUE,
    fallback_to_defaults = TRUE
  )
  
  # Check that positions file was created
  positions_file <- file.path(temp_dir, "euclidify_positions.csv")
  expect_true(file.exists(positions_file))
  
  # Verify file content
  saved_positions <- read.csv(positions_file, row.names = 1)
  expect_equal(dim(saved_positions), dim(result$positions))
  expect_equal(rownames(saved_positions), rownames(result$positions))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("Euclidify clean_intermediate parameter works", {
  test_mat <- create_test_dissimilarity_matrix(10)
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Test with clean_intermediate = FALSE
  result <- Euclidify(
    dissimilarity_matrix = test_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 2),
    n_initial_samples = 2,
    n_adaptive_samples = 2,
    max_cores = 1,
    folds = 5,
    mapping_max_iter = 20,
    verbose = "off",
    clean_intermediate = FALSE,
    fallback_to_defaults = TRUE
  )
  
  # Optimization directory should exist
  optimization_dir <- file.path(temp_dir, "optimization")
  expect_true(dir.exists(optimization_dir))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("Euclidify handles custom parameter ranges", {
  test_mat <- create_test_dissimilarity_matrix(10)
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Test with custom dimension range
  result <- Euclidify(
    dissimilarity_matrix = test_mat,
    output_dir = temp_dir,
    ndim_range = c(3, 6),      # Higher dimensions
    n_initial_samples = 10,
    n_adaptive_samples = 6,
    max_cores = 1,
    folds = 5,
    mapping_max_iter = 50,
    verbose = "off",
    fallback_to_defaults = TRUE
  )
  
  # Should respect the dimension range
  expect_true(ncol(result$positions) >= 3)
  expect_true(ncol(result$positions) <= 6)
  expect_true(result$optimal_params$ndim >= 3)
  expect_true(result$optimal_params$ndim <= 6)
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("Euclidify data characteristics analysis works", {
  # Create a matrix with known properties
  test_mat <- create_test_dissimilarity_matrix(5, with_missing = TRUE)
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  result <- Euclidify(
    dissimilarity_matrix = test_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 3),
    n_initial_samples = 4,
    n_adaptive_samples = 6,
    max_cores = 1,
    folds = 5,
    mapping_max_iter = 50,
    verbose = "off",
    fallback_to_defaults = TRUE
  )
  
  # Check data characteristics
  data_char <- result$data_characteristics
  expect_equal(data_char$n_objects, 5)
  expect_true(data_char$missing_percentage >= 0)
  expect_true(data_char$missing_percentage <= 100)
  expect_true(is.numeric(data_char$non_euclidean_score))
  expect_true(is.numeric(data_char$eigenvalue_dims_90))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})


test_that("Euclidify handles edge case matrices", {
  # Very simple 3x3 matrix
  simple_mat <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), nrow = 3)
  rownames(simple_mat) <- colnames(simple_mat) <- c("A", "B", "C")
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  result <- Euclidify(
    dissimilarity_matrix = simple_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 2),
    n_initial_samples = 4,
    n_adaptive_samples = 4,
    max_cores = 1,
    folds = 3,               # Reduced for small matrix
    mapping_max_iter = 20,
    verbose = "off",
    fallback_to_defaults = TRUE
  )
  
  # Should handle simple case gracefully
  expect_true(is.matrix(result$positions))
  expect_equal(nrow(result$positions), 3)
  expect_equal(ncol(result$positions), 2)
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})


# ---------------------------------------------------------------------------
# Regression tests for the Rcpp-only optimizer backend (2.1.0)
# ---------------------------------------------------------------------------

# Recomputes, in R, the MAE the C++ backend reports in convergence$error:
# the mean absolute error over CONTRIBUTING edges only. An edge contributes
# when it is exact (no threshold), or when its threshold constraint is
# violated (">" with dist < target, "<" with dist > target).
# The optimizer reorders points, so the matrix is realigned to the row order
# of the returned positions before comparing.
recompute_active_mae <- function(dissimilarity_matrix, positions) {
  ord <- rownames(positions)
  mat <- dissimilarity_matrix[ord, ord]
  vals <- mat[upper.tri(mat)]
  idx <- which(upper.tri(mat), arr.ind = TRUE)

  target <- suppressWarnings(as.numeric(gsub("^[<>]", "", vals)))
  thresh <- ifelse(grepl("^>", vals), 1L, ifelse(grepl("^<", vals), -1L, 0L))

  keep <- !is.na(target)
  idx <- idx[keep, , drop = FALSE]
  target <- target[keep]
  thresh <- thresh[keep]

  d <- sqrt(rowSums((positions[idx[, 2], , drop = FALSE] -
                     positions[idx[, 1], , drop = FALSE])^2))
  contributes <- (thresh == 0) |
                 (thresh == 1 & d < target) |
                 (thresh == -1 & d > target)

  list(mae = mean(abs(target - d)[contributes]), n_contributing = sum(contributes))
}

test_that("euclidean_embedding does not modify the caller's initial_positions", {
  # The C++ backend copies the position matrix before optimizing. A missed copy
  # would update the caller's R matrix in place.
  set.seed(101)
  n <- 6
  test_mat <- as.matrix(dist(matrix(runif(n * 2), ncol = 2)))
  dimnames(test_mat) <- list(paste0("Pt", 1:n), paste0("Pt", 1:n))

  init_pos <- matrix(runif(n * 2), nrow = n, ncol = 2)
  rownames(init_pos) <- rownames(test_mat)
  init_pos_before <- init_pos

  euclidean_embedding(test_mat, ndim = 2, mapping_max_iter = 50,
                      k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
                      initial_positions = init_pos, write_positions_to_csv = FALSE)

  expect_identical(init_pos, init_pos_before)
})

test_that("convergence$error matches the MAE recomputed from the returned positions", {
  # The reported error is recorded when the best snapshot is saved; the returned
  # positions are that same snapshot. If the snapshot aliased the live position
  # buffer instead of being deep-copied, later iterations would overwrite it and
  # the two would no longer agree.
  set.seed(102)
  n <- 8
  test_mat <- as.matrix(dist(matrix(runif(n * 2, 0, 10), ncol = 2)))
  dimnames(test_mat) <- list(paste0("Pt", 1:n), paste0("Pt", 1:n))

  result <- euclidean_embedding(test_mat, ndim = 2, mapping_max_iter = 300,
                                k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
                                write_positions_to_csv = FALSE)

  expect_true(is.finite(result$convergence$error))
  expect_equal(recompute_active_mae(test_mat, result$positions)$mae,
               result$convergence$error, tolerance = 1e-10)
})

test_that("thresholded edges contribute only when their constraint is violated", {
  # Same self-consistency check on a fixture mixing "<" and ">" entries. It
  # passes only if the C++ contribution rule for censored edges is unchanged;
  # averaging over all edges instead of the contributing ones gives a different
  # number.
  set.seed(103)
  n <- 7
  num <- as.matrix(dist(matrix(runif(n * 2, 0, 10), ncol = 2)))
  test_mat <- matrix(as.character(round(num, 3)), n, n)
  diag(test_mat) <- "0"
  for (p in list(c(1, 3), c(2, 5), c(4, 7))) {
    test_mat[p[1], p[2]] <- paste0(">", round(num[p[1], p[2]], 3))
    test_mat[p[2], p[1]] <- test_mat[p[1], p[2]]
  }
  for (p in list(c(1, 6), c(3, 5))) {
    test_mat[p[1], p[2]] <- paste0("<", round(num[p[1], p[2]], 3))
    test_mat[p[2], p[1]] <- test_mat[p[1], p[2]]
  }
  dimnames(test_mat) <- list(paste0("V", 1:n), paste0("V", 1:n))

  result <- euclidean_embedding(test_mat, ndim = 2, mapping_max_iter = 300,
                                k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
                                write_positions_to_csv = FALSE)

  active <- recompute_active_mae(test_mat, result$positions)
  expect_true(active$n_contributing < choose(n, 2))   # some constraints satisfied
  expect_equal(active$mae, result$convergence$error, tolerance = 1e-10)
})

test_that("diverging positions degrade gracefully rather than always erroring", {
  # The MAE is computed by multiplying each edge error by its contribution mask,
  # so a non-finite distance propagates NaN into the total even on edges that do
  # not contribute. A NaN MAE compares false against both convergence bounds,
  # lands in the "worsening" branch, and lets the optimizer restore its last
  # good snapshot and return. With a short patience window that happens before
  # the every-10-iterations stability guard can fire.
  set.seed(104)
  base <- as.matrix(dist(matrix(runif(5 * 2, 0, 5), ncol = 2)))
  test_mat <- base * 1e250          # spring forces overflow to non-finite
  dimnames(test_mat) <- list(paste0("A", 1:5), paste0("A", 1:5))

  result <- suppressWarnings(
    euclidean_embedding(test_mat, ndim = 2, mapping_max_iter = 40,
                        k0 = 5, cooling_rate = 0.001, c_repulsion = 0.01,
                        convergence_counter = 3, write_positions_to_csv = FALSE)
  )
  expect_true(result$convergence$achieved)
  expect_true(all(is.finite(result$positions)))

  # With a longer patience window the stability guard is reached first, and the
  # run stops with an error instead. Both outcomes are intentional.
  expect_error(
    suppressWarnings(
      euclidean_embedding(test_mat, ndim = 2, mapping_max_iter = 40,
                          k0 = 5, cooling_rate = 0.001, c_repulsion = 0.01,
                          convergence_counter = 5, write_positions_to_csv = FALSE)
    ),
    "Numerical instability"
  )
})

test_that("returned positions carry no column names from initial_positions", {
  # The backend returns a bare matrix; R/core.R then sets the row names. Column
  # names supplied on initial_positions must not survive into the result.
  set.seed(105)
  n <- 5
  test_mat <- as.matrix(dist(matrix(runif(n * 2), ncol = 2)))
  dimnames(test_mat) <- list(paste0("Pt", 1:n), paste0("Pt", 1:n))

  init_pos <- matrix(runif(n * 2), nrow = n, ncol = 2)
  rownames(init_pos) <- rownames(test_mat)
  colnames(init_pos) <- c("dim_x", "dim_y")

  result <- euclidean_embedding(test_mat, ndim = 2, mapping_max_iter = 50,
                                k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
                                initial_positions = init_pos,
                                write_positions_to_csv = FALSE)

  expect_null(colnames(result$positions))
  expect_equal(sort(rownames(result$positions)), sort(rownames(test_mat)))
})
