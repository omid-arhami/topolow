# tests/testthat/test-s3-methods.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

# Helper function to create a test topolow object
create_test_topolow_object <- function() {
  test_mat <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3)
  rownames(test_mat) <- colnames(test_mat) <- paste0("Point", 1:3)

  euclidean_embedding(
    test_mat,
    ndim = 2,
    mapping_max_iter = 10,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    verbose = FALSE
  )
}

test_that("print.topolow works correctly", {
  result <- create_test_topolow_object()

  expect_s3_class(result, "topolow")

  # Test that print method produces expected output
  expect_output(print(result), "topolow optimization result:")
  expect_output(print(result), "Dimensions:")
  expect_output(print(result), "Iterations:")
  expect_output(print(result), "MAE:")
  expect_output(print(result), "Convergence achieved:")
  expect_output(print(result), "Final convergence error:")

  # Test that print returns the object invisibly
  returned <- capture.output(returned_obj <- print(result))
  expect_identical(returned_obj, result)
})

test_that("summary.topolow works correctly", {
  result <- create_test_topolow_object()

  # Test summary method output
  expect_output(summary(result), "topolow optimization result:")
  expect_output(summary(result), "Parameters:")
  expect_output(summary(result), "k0:")
  expect_output(summary(result), "cooling_rate:")
  expect_output(summary(result), "c_repulsion:")
})

test_that("print.topolow_convergence works correctly", {
  # Create test convergence data
  test_data <- data.frame(
    param1 = rnorm(500, mean = 1, sd = 0.01),
    param2 = rnorm(500, mean = 2, sd = 0.01)
  )

  conv_result <- check_gaussian_convergence(test_data, window_size = 100, tolerance = 0.1)

  expect_s3_class(conv_result, "topolow_convergence")

  # Test print output
  expect_output(print(conv_result), "topolow Convergence Diagnostics")
  expect_output(print(conv_result), "Overall Convergence Achieved:")
  expect_output(print(conv_result), "Mean Vector Converged:")
  expect_output(print(conv_result), "Covariance Matrix Converged:")
  expect_output(print(conv_result), "Final Parameter Means:")
})

test_that("plot.topolow_convergence creates valid plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")

  # Create test convergence data
  test_data <- data.frame(
    param1 = rnorm(200, mean = 1, sd = 0.1),
    param2 = rnorm(200, mean = 2, sd = 0.1)
  )

  conv_result <- check_gaussian_convergence(test_data, window_size = 50)

  # Test that plot method works
  expect_no_error(plot(conv_result))

  # Test with custom parameter names
  expect_no_error(plot(conv_result, param_names = c("Alpha", "Beta")))
})

test_that("topolow_diagnostics S3 methods work", {
  skip_if_not_installed("coda")

  # Create temporary chain files
  temp_dir <- tempdir()
  chain_files <- character(2)

  # Create sample chain data - ensure no missing values
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  set.seed(123)  # For reproducible data
  sample_data <- data.frame(
    log_N = rnorm(100, mean = 1, sd = 0.1),
    log_k0 = rnorm(100, mean = 0, sd = 0.1),
    log_cooling_rate = rnorm(100, mean = -3, sd = 0.1),
    log_c_repulsion = rnorm(100, mean = -2, sd = 0.1),
    NLL = runif(100, 50, 150),
    Holdout_MAE = runif(100, 0.5, 2.0)
  )

  for (i in 1:2) {
    chain_files[i] <- file.path(temp_dir, paste0("chain", i, ".csv"))
    write.csv(sample_data, chain_files[i], row.names = FALSE)
  }

  # Calculate diagnostics - suppress potential warnings about effective sample size
  suppressWarnings({
    diag_results <- calculate_diagnostics(chain_files, mutual_size = 50)
  })

  expect_s3_class(diag_results, "topolow_diagnostics")

  # Test print method
  expect_output(print(diag_results), "topolow Adaptive Sampling Diagnostics")
  expect_output(print(diag_results), "R-hat values")
  expect_output(print(diag_results), "Effective Sample Sizes")

  # Clean up
  unlink(chain_files)
})

test_that("profile_likelihood S3 methods work", {
  # Create test samples
  set.seed(123)  # For reproducible tests
  central_log_N <- 1.5
  test_samples <- data.frame(
    log_N = rnorm(100, mean = central_log_N, sd = 0.2),  # More concentrated
    log_k0 = rnorm(100, mean = 1.0, sd = 0.2),
    log_cooling_rate = rnorm(100, mean = -3.0, sd = 0.3),
    log_c_repulsion = rnorm(100, mean = -1.0, sd = 0.3),
    NLL = runif(100, 20, 100)
  )

  # Suppress expected warnings for small datasets
  suppressWarnings({
    pl_result <- profile_likelihood("log_N", test_samples,
                                    grid_size = 6,
                                    bandwidth_factor = 0.4,  # Larger bandwidth
                                    min_samples = 2)        # Fewer minimum samples
  })
  expect_s3_class(pl_result, "profile_likelihood")

  # Test print method
  expect_output(print(pl_result), "Profile Likelihood Analysis")
  expect_output(print(pl_result), "Parameter: log_N")
  expect_output(print(pl_result), "Grid Points Evaluated:")
  expect_output(print(pl_result), "Bandwidth Used:")
  expect_output(print(pl_result), "Sample Counts per Window")
})

test_that("plot.profile_likelihood creates valid plots", {
  skip_if_not_installed("ggplot2")

  # Create test samples
  set.seed(123)
  test_samples <- data.frame(
    log_N = rnorm(100, mean = 1.5, sd = 0.2),
    log_k0 = rnorm(100, mean = 1.0, sd = 0.2),
    log_cooling_rate = rnorm(100, mean = -3.0, sd = 0.3),
    log_c_repulsion = rnorm(100, mean = -1.0, sd = 0.3),
    NLL = runif(100, 20, 100)
  )

  suppressWarnings({
    pl_result <- profile_likelihood("log_N", test_samples,
                                    grid_size = 6,
                                    bandwidth_factor = 0.4,
                                    min_samples = 2)
  })
  LL_max <- max(-test_samples$NLL)

  # Test plot creation - suppress ggplot warnings about missing values
  suppressWarnings({
    p <- plot(pl_result, LL_max)
  })
  expect_s3_class(p, "ggplot")

  # Test plot with saving (to temporary file)
  temp_dir <- tempdir()
  suppressWarnings({
    expect_no_error(
      plot(pl_result, LL_max, save_plot = TRUE, output_dir = temp_dir)
    )
  })

  # Check that file was created
  expected_file <- file.path(temp_dir, "profile_likelihood_log_N.pdf")
  expect_true(file.exists(expected_file))
  unlink(expected_file)
})

test_that("parameter_sensitivity S3 methods work", {
  # Create test samples
  test_samples <- data.frame(
    log_N = log(runif(100, 2, 6)),
    log_k0 = log(runif(100, 1, 5)),
    Holdout_MAE = runif(100, 0.5, 3.0)
  )

  # Calculate parameter sensitivity
  sens_result <- parameter_sensitivity_analysis("log_N", test_samples, bins = 10)

  expect_s3_class(sens_result, "parameter_sensitivity")

  # Test print method
  expect_output(print(sens_result), "Parameter Sensitivity Analysis")
  expect_output(print(sens_result), "Parameter Analyzed: log_N")
  expect_output(print(sens_result), "Number of Bins:")
  expect_output(print(sens_result), "Minimum MAE Found:")
  expect_output(print(sens_result), "Performance Threshold")
  expect_output(print(sens_result), "Sample Counts per Bin")
})

test_that("plot.parameter_sensitivity creates valid plots", {
  skip_if_not_installed("ggplot2")

  # Create test samples
  test_samples <- data.frame(
    log_k0 = log(runif(80, 1, 5)),
    Holdout_MAE = runif(80, 0.5, 3.0)
  )

  sens_result <- parameter_sensitivity_analysis("log_k0", test_samples, bins = 8)

  # Test plot creation
  p <- plot(sens_result)
  expect_s3_class(p, "ggplot")

  # Test plot with custom y-limit factor
  p_custom <- plot(sens_result, y_limit_factor = 1.2)
  expect_s3_class(p_custom, "ggplot")

  # Test plot with saving
  temp_dir <- tempdir()
  expect_no_error(
    plot(sens_result, save_plot = TRUE, output_dir = temp_dir)
  )

  expected_file <- file.path(temp_dir, "parameter_sensitivity_log_k0.pdf")
  expect_true(file.exists(expected_file))
  unlink(expected_file)
})

test_that("S3 method error handling works", {
  # Test plot.profile_likelihood with missing output_dir
  set.seed(123)  # For reproducible tests
  central_log_N <- 1.5
  test_samples <- data.frame(
    log_N = rnorm(100, mean = central_log_N, sd = 0.3),  # More concentrated around mean
    log_k0 = rnorm(100, mean = 1.0, sd = 0.2),
    log_cooling_rate = rnorm(100, mean = -3.0, sd = 0.3),
    log_c_repulsion = rnorm(100, mean = -1.0, sd = 0.3),
    NLL = runif(100, 20, 100)
  )
  suppressWarnings({
  pl_result <- profile_likelihood("log_N", test_samples,
                                  grid_size = 5,
                                  bandwidth_factor = 0.7,  # Even larger bandwidth
                                  min_samples = 2)
  })
  LL_max <- max(-test_samples$NLL)

  expect_error(
    plot(pl_result, LL_max, save_plot = TRUE),
    "'output_dir' must be provided when save_plot is TRUE"
  )

  # Test plot.parameter_sensitivity with missing output_dir
  sens_result <- parameter_sensitivity_analysis("log_N", test_samples, bins = 5, mae_col = "NLL")

  expect_error(
    plot(sens_result, save_plot = TRUE),
    "'output_dir' must be provided when save_plot is TRUE"
  )
})
