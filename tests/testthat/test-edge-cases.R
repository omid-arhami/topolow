# tests/testthat/test-edge-cases.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("euclidean_embedding handles degenerate matrices", {
  # All points at the same location (zero distances) - suppress expected warning
  zero_mat <- matrix(0, 3, 3)
  rownames(zero_mat) <- colnames(zero_mat) <- paste0("Point", 1:3)

  suppressWarnings({
    result_zero <- euclidean_embedding(
      zero_mat, ndim = 2, mapping_max_iter = 20,
      k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
      verbose = FALSE
    )
  })

  expect_s3_class(result_zero, "topolow")
  expect_true(is.finite(result_zero$mae))

  # Matrix with very large values
  large_mat <- matrix(runif(9, 1000, 10000), 3, 3)
  large_mat[lower.tri(large_mat)] <- t(large_mat)[lower.tri(large_mat)]
  diag(large_mat) <- 0

  expect_no_error(
    result_large <- euclidean_embedding(
      large_mat, ndim = 2, mapping_max_iter = 20,
      k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
      verbose = FALSE
    )
  )

  # Matrix with very small values
  small_mat <- matrix(runif(9, 1e-6, 1e-3), 3, 3)
  small_mat[lower.tri(small_mat)] <- t(small_mat)[lower.tri(small_mat)]
  diag(small_mat) <- 0

  expect_no_error(
    result_small <- euclidean_embedding(
      small_mat, ndim = 2, mapping_max_iter = 20,
      k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
      verbose = FALSE
    )
  )
})

test_that("euclidean_embedding handles matrices with only thresholds", {
  # Matrix with only threshold values - suppress expected warning
  threshold_mat <- matrix(c("0", ">5", "<10", ">5", "0", ">20", "<10", ">20", "0"), 3, 3)
  rownames(threshold_mat) <- colnames(threshold_mat) <- paste0("Point", 1:3)

  suppressWarnings({
    result_thresh <- euclidean_embedding(
      threshold_mat, ndim = 2, mapping_max_iter = 30,
      k0 = 2.0, cooling_rate = 0.01, c_repulsion = 0.05,
      verbose = FALSE
    )
  })

  expect_s3_class(result_thresh, "topolow")
  expect_true(is.matrix(result_thresh$positions))
  expect_true(is.matrix(result_thresh$est_distances))
})

test_that("euclidean_embedding handles extremely sparse matrices", {
  # Matrix with only one measurement
  sparse_mat <- matrix(NA, 4, 4)
  sparse_mat[1, 2] <- sparse_mat[2, 1] <- 5
  diag(sparse_mat) <- 0
  rownames(sparse_mat) <- colnames(sparse_mat) <- paste0("Point", 1:4)

  result_sparse <- euclidean_embedding(
    sparse_mat, ndim = 2, mapping_max_iter = 50,
    k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.1,
    verbose = FALSE
  )

  expect_s3_class(result_sparse, "topolow")
  # All points should be positioned (repulsive forces should spread them out)
  expect_true(all(is.finite(result_sparse$positions)))
})

test_that("error_calculator_comparison handles edge cases", {
  # Test with all identical values
  identical_true <- matrix(5, 3, 3)
  identical_pred <- matrix(5, 3, 3)
  diag(identical_true) <- diag(identical_pred) <- 0

  errors_identical <- error_calculator_comparison(identical_pred, identical_true)
  expect_equal(errors_identical$Completeness, 1)

  # Test with matrices containing Inf values
  inf_true <- matrix(c(0, 1, Inf, 1, 0, 2, Inf, 2, 0), 3, 3)
  inf_pred <- matrix(c(0, 1.1, 3, 1.1, 0, 2.1, 3, 2.1, 0), 3, 3)

  expect_no_error(
    errors_inf <- error_calculator_comparison(inf_pred, inf_true)
  )

  # Test with matrices where all predictions are NA
  na_pred <- matrix(NA, 3, 3)
  valid_true <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3)

  errors_na <- error_calculator_comparison(na_pred, valid_true)
  # When no predictions can be made, completeness should be 0
  expect_equal(errors_na$Completeness, 0)
})


test_that("functions handle single-element inputs", {
  # Test extract_numeric_values with single element
  single_val <- extract_numeric_values("<5")
  expect_equal(single_val, 5)

  # Test clean_data with single element
  single_clean <- clean_data(c(10))
  expect_equal(single_clean, 10)

  # Test with minimal viable matrix (3x3 - smallest that reliably works)
  minimal_mat <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), 3, 3)
  rownames(minimal_mat) <- colnames(minimal_mat) <- paste0("Point", 1:3)

  expect_no_error(
    minimal_result <- euclidean_embedding(
      minimal_mat, ndim = 1, mapping_max_iter = 10,
      k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
      verbose = FALSE
    )
  )
})

test_that("visualization functions handle empty or minimal data", {
  skip_if_not_installed("ggplot2")

  # Test network plot with minimal connectivity
  minimal_adj <- matrix(FALSE, 3, 3)
  minimal_adj[1, 2] <- minimal_adj[2, 1] <- TRUE
  rownames(minimal_adj) <- colnames(minimal_adj) <- paste0("Node", 1:3)

  minimal_network <- list(
    adjacency = minimal_adj,
    connectivity = data.frame(
      node = paste0("Node", 1:3),
      degree = c(1, 1, 0),
      completeness = c(0.5, 0.5, 0)
    ),
    summary = list(
      n_points = 3,
      n_measurements = 1,
      completeness = 0.33
    )
  )

  expect_no_error(
    p_minimal <- plot_network_structure(minimal_network)
  )
  expect_s3_class(p_minimal, "ggplot")
})

test_that("parameter optimization handles extreme parameter ranges", {
  # Test with very small parameter ranges
  small_test_mat <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), 3, 3)

  expect_no_error(
    small_range_results <- initial_parameter_optimization(
      dissimilarity_matrix = small_test_mat,
      mapping_max_iter = 20,
      relative_epsilon = 1e-2,
      convergence_counter = 2,
      scenario_name = "edge_test",
      N_min = 2, N_max = 2,  # Fixed dimension (should be allowed)
      k0_min = 0.9, k0_max = 1.1,  # Very narrow range
      c_repulsion_min = 0.009, c_repulsion_max = 0.011,
      cooling_rate_min = 0.009, cooling_rate_max = 0.011,
      num_samples = 2,
      folds = 2,
      max_cores = 1,
      write_files = FALSE
    )
  )

  expect_true(is.data.frame(small_range_results))
})

test_that("adaptive sampling handles problematic samples", {
  # Create samples with extreme NLL values
  temp_samples_file <- tempfile(fileext = ".csv")
  extreme_samples <- data.frame(
    log_N = log(c(2, 3, 4)),
    log_k0 = log(c(1, 2, 3)),
    log_cooling_rate = log(c(0.01, 0.02, 0.03)),
    log_c_repulsion = log(c(0.1, 0.2, 0.3)),
    NLL = c(1e6, 1e-6, Inf),  # Extreme values
    Holdout_MAE = c(100, 0.001, 50)
  )
  write.csv(extreme_samples, temp_samples_file, row.names = FALSE)

  # Test matrix for likelihood calculation
  test_mat <- matrix(runif(16, 1, 5), 4, 4)
  diag(test_mat) <- 0

  # Should handle extreme values gracefully
  expect_no_error(
    result <- adaptive_MC_sampling(
      samples_file = temp_samples_file,
      dissimilarity_matrix = test_mat,
      iterations = 1,
      mapping_max_iter = 10,
      relative_epsilon = 1e-2,
      folds = 2,
      scenario_name = "edge_test",
      verbose = FALSE
    )
  )

  unlink(temp_samples_file)
})

test_that("file I/O edge cases are handled", {
  skip_if_not_installed("filelock")

  # Test with non-existent directory for output
  non_existent_dir <- file.path(tempdir(), "non_existent_subdir")

  test_mat <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3)

  # Function should create directory if needed
  expect_no_error(
    result <- euclidean_embedding(
      test_mat, ndim = 2, mapping_max_iter = 10,
      k0 = 1.0, cooling_rate = 0.01, c_repulsion = 0.01,
      write_positions_to_csv = TRUE,
      output_dir = non_existent_dir,
      verbose = FALSE
    )
  )

  expect_true(dir.exists(non_existent_dir))
  unlink(non_existent_dir, recursive = TRUE)
})

test_that("numerical stability with extreme inputs", {
  # Test with very high dimensional embedding
  test_mat <- matrix(runif(25, 1, 10), 5, 5)
  test_mat[lower.tri(test_mat)] <- t(test_mat)[lower.tri(test_mat)]
  diag(test_mat) <- 0

  # High dimensional embedding might be unstable but should not crash
  expect_no_error(
    high_dim_result <- euclidean_embedding(
      test_mat, ndim = 4,  # High dimension relative to data
      mapping_max_iter = 30,
      k0 = 0.1,  # Low spring constant
      cooling_rate = 0.001,  # Slow cooling
      c_repulsion = 0.001,  # Low repulsion
      verbose = FALSE
    )
  )

  # Check that result is still valid
  expect_true(all(is.finite(high_dim_result$positions)))
})

test_that("error metrics handle pathological cases", {
  # Test prediction interval with perfect predictions
  true_mat <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3)
  perfect_pred <- true_mat  # Perfect predictions

  expect_no_error(
    perfect_interval <- calculate_prediction_interval(true_mat, perfect_pred)
  )
  expect_true(perfect_interval >= 0)

  # Test with insufficient data points - matrices with all NA values
  tiny_true <- matrix(NA, 2, 2)
  tiny_pred <- matrix(NA, 2, 2)

  expect_error(
    calculate_prediction_interval(tiny_true, tiny_pred),
    "Fewer than 2 valid paired observations"
  )
})
