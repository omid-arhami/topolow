# =============================================================================
# Tests for incremental_embedding()
# =============================================================================
# File: tests/testthat/test-incremental.R

library(testthat)
library(topolow)

test_that("incremental_embedding adds single new point correctly", {
  # Create initial map with 4 points
  initial_dist <- matrix(c(
    0, 2, 3, 4,
    2, 0, 2.5, 3.5,
    3, 2.5, 0, 2,
    4, 3.5, 2, 0
  ), nrow = 4, byrow = TRUE)
  rownames(initial_dist) <- colnames(initial_dist) <- c("A", "B", "C", "D")
  
  initial_map <- euclidean_embedding(
    dissimilarity_matrix = initial_dist,
    ndim = 2,
    mapping_max_iter = 500,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  # Add one new point with measurements to existing points
  new_data <- data.frame(
    object = c("E", "E", "E"),
    reference = c("A", "B", "C"),
    value = c(2.0, 2.5, 3.0)
  )
  
  updated_map <- incremental_embedding(
    fixed_positions = initial_map$positions,
    new_measurements = new_data,
    mapping_max_iter = 500,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  # Check that we have 5 points now

  expect_equal(nrow(updated_map$positions), 5)
  
  # Check that fixed points haven't moved significantly
  fixed_names <- c("A", "B", "C", "D")
  for (name in fixed_names) {
    original_pos <- initial_map$positions[name, ]
    updated_pos <- updated_map$positions[name, ]
    expect_equal(original_pos, updated_pos, tolerance = 1e-10)
  }
  
  # Check that new point exists
  expect_true("E" %in% rownames(updated_map$positions))
  
  # Check incremental_info
  expect_equal(updated_map$incremental_info$n_fixed_points, 4)
  expect_equal(updated_map$incremental_info$n_new_points, 1)
  expect_equal(updated_map$incremental_info$new_point_names, "E")
})


test_that("incremental_embedding adds multiple new points correctly", {
  # Create initial map
  initial_dist <- matrix(c(
    0, 2, 3,
    2, 0, 4,
    3, 4, 0
  ), nrow = 3, byrow = TRUE)
  rownames(initial_dist) <- colnames(initial_dist) <- c("A", "B", "C")
  
  initial_map <- euclidean_embedding(
    dissimilarity_matrix = initial_dist,
    ndim = 2,
    mapping_max_iter = 500,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  # Add two new points
  new_data <- data.frame(
    object = c("D", "D", "E", "E", "D"),
    reference = c("A", "B", "A", "C", "E"),
    value = c(1.5, 2.0, 2.5, 1.8, 3.0)
  )
  
  updated_map <- incremental_embedding(
    fixed_positions = initial_map$positions,
    new_measurements = new_data,
    mapping_max_iter = 500,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  expect_equal(nrow(updated_map$positions), 5)
  expect_equal(updated_map$incremental_info$n_new_points, 2)
  expect_true(all(c("D", "E") %in% updated_map$incremental_info$new_point_names))
})


test_that("incremental_embedding warns when no new points", {
  # Create initial map
  initial_dist <- matrix(c(0, 2, 2, 0), nrow = 2)
  rownames(initial_dist) <- colnames(initial_dist) <- c("A", "B")
  
  initial_map <- euclidean_embedding(
    dissimilarity_matrix = initial_dist,
    ndim = 2,
    mapping_max_iter = 100,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  # Measurements only between existing points
  new_data <- data.frame(
    object = c("A"),
    reference = c("B"),
    value = c(2.5)
  )
  
  expect_warning(
    updated_map <- incremental_embedding(
      fixed_positions = initial_map$positions,
      new_measurements = new_data,
      verbose = FALSE
    ),
    "No new points found"
  )
  
  # Should return original positions unchanged
  expect_equal(updated_map$positions, initial_map$positions)
})


test_that("incremental_embedding warns when new point has insufficient measurements", {
  # Create initial map
  initial_dist <- matrix(c(
    0, 2, 3,
    2, 0, 4,
    3, 4, 0
  ), nrow = 3)
  rownames(initial_dist) <- colnames(initial_dist) <- c("A", "B", "C")
  
  initial_map <- euclidean_embedding(
    dissimilarity_matrix = initial_dist,
    ndim = 2,  # 2 dimensions
    mapping_max_iter = 200,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  # New point with only 1 measurement (less than ndim=2)
  new_data <- data.frame(
    object = c("D"),
    reference = c("A"),
    value = c(2.0)
  )
  
  expect_warning(
    updated_map <- incremental_embedding(
      fixed_positions = initial_map$positions,
      new_measurements = new_data,
      verbose = FALSE
    ),
    "fewer measurements"
  )
  
  # Check warning is recorded
  expect_true(length(updated_map$incremental_info$warnings) > 0)
})


test_that("incremental_embedding warns when new point has no anchor to fixed points", {
  # Create initial map
  initial_dist <- matrix(c(0, 2, 2, 0), nrow = 2)
  rownames(initial_dist) <- colnames(initial_dist) <- c("A", "B")
  
  initial_map <- euclidean_embedding(
    dissimilarity_matrix = initial_dist,
    ndim = 2,
    mapping_max_iter = 100,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  # Two new points that only connect to each other, not to fixed points
  # Plus one that does connect
  new_data <- data.frame(
    object = c("C", "D", "E"),
    reference = c("A", "E", "D"),  # C->A is anchored, but D and E only connect to each other
    value = c(2.0, 1.5, 1.5)
  )
  
  expect_warning(
    updated_map <- incremental_embedding(
      fixed_positions = initial_map$positions,
      new_measurements = new_data,
      verbose = FALSE
    ),
    "no measurements to fixed points"
  )
})


test_that("incremental_embedding handles threshold indicators", {
  # Create initial map
  initial_dist <- matrix(c(0, 2, 2, 0), nrow = 2)
  rownames(initial_dist) <- colnames(initial_dist) <- c("A", "B")
  
  initial_map <- euclidean_embedding(
    dissimilarity_matrix = initial_dist,
    ndim = 2,
    mapping_max_iter = 100,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  # New measurements with threshold indicators
  new_data <- data.frame(
    object = c("C", "C"),
    reference = c("A", "B"),
    value = c(">1.5", "<3.0")  # C is at least 1.5 from A, at most 3.0 from B
  )
  
  updated_map <- incremental_embedding(
    fixed_positions = initial_map$positions,
    new_measurements = new_data,
    mapping_max_iter = 500,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  expect_equal(nrow(updated_map$positions), 3)
  expect_true("C" %in% rownames(updated_map$positions))
})


test_that("incremental_embedding returns topolow class", {
  initial_dist <- matrix(c(0, 2, 2, 0), nrow = 2)
  rownames(initial_dist) <- colnames(initial_dist) <- c("A", "B")
  
  initial_map <- euclidean_embedding(
    dissimilarity_matrix = initial_dist,
    ndim = 2,
    mapping_max_iter = 100,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01
  )
  
  new_data <- data.frame(
    object = c("C"),
    reference = c("A"),
    value = c(2.0)
  )
  
  suppressWarnings({
    updated_map <- incremental_embedding(
      fixed_positions = initial_map$positions,
      new_measurements = new_data
    )
  })
  
  expect_s3_class(updated_map, "topolow")
  expect_true(all(c("positions", "est_distances", "mae", "iter", 
                    "parameters", "convergence") %in% names(updated_map)))
})


test_that("incremental_embedding validates inputs", {
  expect_error(
    incremental_embedding(
      fixed_positions = "not a matrix",
      new_measurements = data.frame(object = "A", reference = "B", value = 1)
    ),
    "must be a matrix"
  )
  
  expect_error(
    incremental_embedding(
      fixed_positions = matrix(1:4, 2, 2),  # No row names
      new_measurements = data.frame(object = "A", reference = "B", value = 1)
    ),
    "must have row names"
  )
  
  mat <- matrix(1:4, 2, 2)
  rownames(mat) <- c("A", "B")
  
  expect_error(
    incremental_embedding(
      fixed_positions = mat,
      new_measurements = data.frame(obj = "C", ref = "A", val = 1)  # Wrong column names
    ),
    "missing required columns"
  )
})