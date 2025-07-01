# tests/testthat/test-core.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("Package loads without error", {
  library(topolow)
  expect_true(TRUE)
}
)

# Helper function to create test matrices
create_test_matrix <- function(size=3, with_threshold=FALSE) {
  if(with_threshold) {
    matrix(c(0, ">2", 3, ">2", 0, 4, 3, 4, 0), nrow=size)
  } else {
    mat <- as.matrix(dist(matrix(runif(size*2), ncol=2)))
    colnames(mat) <- rownames(mat) <- paste0("Point", 1:size)
    mat
  }
}

test_that("create_topolow_map handles input validation correctly", {
  # Basic matrix validation
  expect_error(
    create_topolow_map(distance_matrix = "not a matrix", ndim = 2, mapping_max_iter = 10,
                 k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01, write_positions_to_csv=FALSE),
    "distance_matrix must be a matrix"
  )
  
  # Non-square matrix
  expect_error(
    create_topolow_map(distance_matrix = matrix(1:6, nrow=2), ndim = 2, mapping_max_iter = 10,
                 k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01, write_positions_to_csv=FALSE),
    "distance_matrix must be square"
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
  for(param in params_to_test) {
    args <- list(
      distance_matrix = test_mat,
      ndim = 2, mapping_max_iter = 10,
      k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01
    )
    args[[names(param)[1]]] <- param[[1]]
    expect_error(do.call(create_topolow_map, args), param$msg)
  }
  
  # Test high k0 warning
  expect_warning(
    create_topolow_map(distance_matrix = test_mat, ndim = 2, mapping_max_iter = 10,
                 k0 = 35, cooling_rate = 0.01, c_repulsion = 0.01, write_positions_to_csv=FALSE),
    "High k0 value"
  )
})

test_that("create_topolow_map handles initial positions correctly", {
  test_mat <- create_test_matrix()
  n <- nrow(test_mat)
  ndim <- 2
  
  # Test valid initial positions
  init_pos <- matrix(runif(n*ndim), nrow=n, ncol=ndim)
  rownames(init_pos) <- rownames(test_mat)
  result <- create_topolow_map(test_mat, ndim=ndim, mapping_max_iter=10,
                        k0=1.0, cooling_rate=0.01, c_repulsion=0.01,
                        initial_positions=init_pos, write_positions_to_csv=FALSE)
  expect_equal(dim(result$positions), c(n, ndim))
  
  # Test invalid dimensions
  wrong_pos <- matrix(runif((n+1)*ndim), nrow=n+1, ncol=ndim)
  expect_error(
    create_topolow_map(test_mat, ndim=ndim, mapping_max_iter=10,
                 k0=1.0, cooling_rate=0.01, c_repulsion=0.01,
                 initial_positions=wrong_pos, write_positions_to_csv=FALSE),
    "initial_positions must have same number of rows"
  )
})

test_that("create_topolow_map convergence behavior works correctly", {
  test_mat <- create_test_matrix()
  
  # Test convergence tracking
  result <- create_topolow_map(test_mat, ndim=2, mapping_max_iter=10,
                        k0=1.0, cooling_rate=0.01, c_repulsion=0.01,
                        relative_epsilon=1e-4, convergence_counter=5, 
                        write_positions_to_csv=FALSE)
  
  expect_true(!is.null(result$convergence))
  expect_true(is.logical(result$convergence$achieved))
  expect_true(is.numeric(result$convergence$error))
})

test_that("create_topolow_map force calculations preserve distance relationships", {
  # Create triangle with known properties
  test_mat <- matrix(c(0,1,2, 1,0,1, 2,1,0), nrow=3)
  colnames(test_mat) <- rownames(test_mat) <- c("A","B","C")
  
  result <- create_topolow_map(test_mat, ndim=2, mapping_max_iter=10,
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

test_that("create_topolow_map handles missing values appropriately", {
  test_mat <- create_test_matrix()
  test_mat[1,2] <- test_mat[2,1] <- NA
  
  result <- create_topolow_map(test_mat, ndim=2, mapping_max_iter=10,
                        k0=1.0, cooling_rate=0.01, c_repulsion=0.01, 
                        write_positions_to_csv=FALSE)
  
  # Missing values should use repulsive forces
  expect_true(!is.na(result$est_distances[1,2]))
  expect_equal(result$est_distances[1,2], result$est_distances[2,1])  # Symmetry
})

test_that("create_topolow_map returns correct object structure", {
  result <- create_topolow_map(create_test_matrix(), ndim=2, mapping_max_iter=10,
                        k0=1.0, cooling_rate=0.01, c_repulsion=0.01, write_positions_to_csv=FALSE)
  
  expected_elements <- c("positions", "est_distances", "mae", "iter",
                         "parameters", "convergence")
  expect_true(all(expected_elements %in% names(result)))
  expect_s3_class(result, "topolow")
  expect_true(is.numeric(result$mae))
})