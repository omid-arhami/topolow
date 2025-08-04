# tests/testthat/test-diagnostics.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("error_calculator_comparison works correctly", {
  true_mat <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3)
  pred_mat <- true_mat + 0.1
  input_mat <- true_mat
  input_mat[1, 3] <- input_mat[3, 1] <- NA # Create a holdout set

  errors <- error_calculator_comparison(pred_mat, true_mat, input_mat)

  expect_true(is.list(errors))
  expect_true(all(c("report_df", "Completeness") %in% names(errors)))
  expect_true(is.data.frame(errors$report_df))
  # Check that out-of-sample error is calculated only for the NA part
  expect_equal(sum(!is.na(errors$report_df$OutSampleError)), 2)
  # Check that in-sample error is calculated for the non-NA part
  expect_equal(sum(!is.na(errors$report_df$InSampleError)), 7) # Includes diagonal
})

test_that("analyze_network_structure returns correct stats", {
  dist_mat <- matrix(runif(25), 5, 5)
  diag(dist_mat) <- 0
  dist_mat[1, 3] <- dist_mat[3, 1] <- NA

  metrics <- analyze_network_structure(dist_mat)

  expect_true(is.list(metrics))
  expect_true(all(c("adjacency", "connectivity", "summary") %in% names(metrics)))
  expect_equal(metrics$summary$n_measurements, (25 - 5 - 2) / 2) # (all - diag - 2 NAs) / 2
  expect_equal(metrics$connectivity$degree[1], 3) # Node 1 connects to 2, 4, 5
})

test_that("check_gaussian_convergence identifies convergence", {
  # Converged data
  conv_data <- as.data.frame(matrix(rnorm(1000 * 2, mean = 5, sd = 0.001), ncol = 2))
  conv_results <- check_gaussian_convergence(conv_data, window_size = 200, tolerance = 0.3)
  expect_true(conv_results$converged)
  expect_s3_class(conv_results, "topolow_convergence")

  # Non-converged data (trending)
  non_conv_data <- as.data.frame(matrix(1:(1000 * 2), ncol = 2))
  non_conv_results <- check_gaussian_convergence(non_conv_data, window_size = 100, tolerance = 0.001)
  expect_false(non_conv_results$converged)
})
