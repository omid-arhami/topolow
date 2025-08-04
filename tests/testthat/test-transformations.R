# tests/testthat/test-transformations.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("coordinates_to_matrix works correctly", {
  # Test with matrix input
  coords_matrix <- matrix(rnorm(20), 5, 4)
  dist_mat <- coordinates_to_matrix(coords_matrix)

  expect_true(is.matrix(dist_mat))
  expect_equal(dim(dist_mat), c(5, 5))
  expect_true(isSymmetric(dist_mat))
  expect_equal(as.numeric(diag(dist_mat)), rep(0, 5))
  expect_true(all(dist_mat >= 0))

  # Test with data frame input
  coords_df <- as.data.frame(coords_matrix)
  dist_mat_df <- coordinates_to_matrix(coords_df)
  expect_equal(dist_mat, dist_mat_df)

  # Test with row names
  rownames(coords_matrix) <- paste0("Point", 1:5)
  dist_mat_named <- coordinates_to_matrix(coords_matrix)
  expect_equal(rownames(dist_mat_named), paste0("Point", 1:5))
  expect_equal(colnames(dist_mat_named), paste0("Point", 1:5))
})

test_that("coordinates_to_matrix handles edge cases", {
  # Test with 2D coordinates
  # Create points: (0, 0) and (1, 1) - distance should be sqrt(2)
  coords_2d <- matrix(c(0, 0, 1, 1), 2, 2, byrow = TRUE)
  dist_mat_2d <- coordinates_to_matrix(coords_2d)
  expect_equal(dist_mat_2d[1, 2], sqrt(2), tolerance = 1e-10)

  # Test with single point (should work but give 1x1 matrix)
  coords_single <- matrix(c(1, 2, 3), 1, 3)
  dist_mat_single <- coordinates_to_matrix(coords_single)
  expect_equal(dim(dist_mat_single), c(1, 1))
  expect_equal(dist_mat_single[1, 1], 0)
})

test_that("coordinates_to_matrix input validation", {
  # Test invalid inputs
  expect_error(coordinates_to_matrix("not a matrix"),
               "positions must be a matrix or a data frame")
  expect_error(coordinates_to_matrix(list(a = 1, b = 2)),
               "positions must be a matrix or a data frame")
})

test_that("log_transform_parameters works with valid input", {
  # Create temporary input file
  temp_input <- tempfile(fileext = ".csv")
  test_params <- data.frame(
    N = c(2, 3, 4),
    k0 = c(1.0, 1.5, 2.0),
    cooling_rate = c(0.01, 0.02, 0.03),
    c_repulsion = c(0.1, 0.2, 0.3),
    Holdout_MAE = c(0.5, 0.6, 0.7),
    NLL = c(100, 110, 120)
  )
  write.csv(test_params, temp_input, row.names = FALSE)

  # Test transformation without output file
  result <- log_transform_parameters(temp_input)

  expect_true(is.data.frame(result))
  expect_true(all(c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion") %in% names(result)))
  expect_false(any(c("N", "k0", "cooling_rate", "c_repulsion") %in% names(result)))
  expect_true(all(c("Holdout_MAE", "NLL") %in% names(result)))

  # Check that log transformation was applied correctly
  expect_equal(result$log_N[1], log(2))
  expect_equal(result$log_k0[2], log(1.5))

  unlink(temp_input)
})

test_that("log_transform_parameters handles output file", {
  # Create temporary input and output files
  temp_input <- tempfile(fileext = ".csv")
  temp_output <- tempfile(fileext = ".csv")

  test_params <- data.frame(
    N = c(2, 3),
    k0 = c(1.0, 1.5),
    cooling_rate = c(0.01, 0.02),
    c_repulsion = c(0.1, 0.2),
    Holdout_MAE = c(0.5, 0.6),
    NLL = c(100, 110)
  )
  write.csv(test_params, temp_input, row.names = FALSE)

  # Test with output file
  expect_message(
    result <- log_transform_parameters(temp_input, temp_output),
    "Log transformed parameters"
  )

  expect_true(file.exists(temp_output))

  # Read the output file and verify
  output_data <- read.csv(temp_output)
  expect_equal(names(output_data), names(result))
  expect_equal(nrow(output_data), 2)

  unlink(c(temp_input, temp_output))
})

test_that("log_transform_parameters handles edge cases", {
  temp_input <- tempfile(fileext = ".csv")

  # Test with no transformable parameters
  no_transform_data <- data.frame(
    other_param = c(1, 2, 3),
    Holdout_MAE = c(0.5, 0.6, 0.7)
  )
  write.csv(no_transform_data, temp_input, row.names = FALSE)

  expect_message(
    result <- log_transform_parameters(temp_input),
    "No parameters found to transform"
  )
  expect_equal(result, no_transform_data)

  # Test with NA values
  na_data <- data.frame(
    N = c(2, NA, 4),
    k0 = c(1.0, 1.5, 2.0),
    Holdout_MAE = c(0.5, 0.6, 0.7)
  )
  write.csv(na_data, temp_input, row.names = FALSE)

  result_na <- log_transform_parameters(temp_input)
  expect_equal(nrow(result_na), 2) # NA row should be removed

  unlink(temp_input)
})

test_that("log_transform_parameters input validation", {
  # Test with non-existent file
  expect_error(log_transform_parameters("nonexistent.csv"),
               "Input file not found")

  # Test with invalid output file parameter
  temp_input <- tempfile(fileext = ".csv")
  write.csv(data.frame(N = 1, k0 = 1), temp_input, row.names = FALSE)

  expect_error(log_transform_parameters(temp_input, output_file = c("file1", "file2")),
               "'output_file' must be a single character string")

  # Test with non-positive values
  invalid_data <- data.frame(
    N = c(2, -1, 4),  # Negative value
    k0 = c(1.0, 1.5, 0),  # Zero value
    Holdout_MAE = c(0.5, 0.6, 0.7)
  )
  write.csv(invalid_data, temp_input, row.names = FALSE)

  expect_error(log_transform_parameters(temp_input),
               "Non-positive values found")

  unlink(temp_input)
})
