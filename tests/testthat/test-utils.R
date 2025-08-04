# tests/testthat/test-utils.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("extract_numeric_values handles mixed data correctly", {
  # Test with various threshold indicators
  mixed_data <- c(10, 20, "<5", ">100", 50, "25")
  result <- extract_numeric_values(mixed_data)

  expect_equal(result[1], 10)
  expect_equal(result[2], 20)
  expect_equal(result[3], 5)   # "<5" should become 5
  expect_equal(result[4], 100) # ">100" should become 100
  expect_equal(result[5], 50)
  expect_equal(result[6], 25)  # "25" should become 25
})

test_that("extract_numeric_values handles edge cases", {
  # Test with NA values
  data_with_na <- c(10, NA, "<5", ">100")
  result_na <- extract_numeric_values(data_with_na)
  expect_true(is.na(result_na[2]))
  expect_equal(result_na[3], 5)

  # Test with empty vector
  empty_result <- extract_numeric_values(c())
  expect_length(empty_result, 0)

  # Test with all numeric
  all_numeric <- c(1, 2, 3, 4, 5)
  result_numeric <- extract_numeric_values(all_numeric)
  expect_equal(result_numeric, all_numeric)

  # Test with all character (no threshold indicators)
  all_char <- c("10", "20", "30")
  result_char <- extract_numeric_values(all_char)
  expect_equal(result_char, c(10, 20, 30))

  # Test with invalid character values - suppress expected warnings
  invalid_char <- c("abc", "def", "10")
  suppressWarnings({
    result_invalid <- extract_numeric_values(invalid_char)
  })
  expect_true(is.na(result_invalid[1]))
  expect_true(is.na(result_invalid[2]))
  expect_equal(result_invalid[3], 10)
})

test_that("extract_numeric_values handles complex thresholds", {
  # Test with decimal thresholds
  decimal_thresholds <- c("<0.5", ">10.25", "5.75")
  result_decimal <- extract_numeric_values(decimal_thresholds)
  expect_equal(result_decimal[1], 0.5)
  expect_equal(result_decimal[2], 10.25)
  expect_equal(result_decimal[3], 5.75)

  # Test with scientific notation
  sci_notation <- c("<1e-3", ">1E5", "2.5e2")
  result_sci <- extract_numeric_values(sci_notation)
  expect_equal(result_sci[1], 0.001)
  expect_equal(result_sci[2], 100000)
  expect_equal(result_sci[3], 250)
})

test_that("create_cv_folds creates proper fold structure", {
  # Create test dissimilarity matrix
  test_mat <- matrix(runif(25, 1, 10), 5, 5)
  diag(test_mat) <- 0
  rownames(test_mat) <- colnames(test_mat) <- paste0("Point", 1:5)

  # Test basic fold creation
  folds <- create_cv_folds(test_mat, n_folds = 3, random_seed = 123)

  expect_length(folds, 3)
  expect_true(all(sapply(folds, is.list)))
  expect_true(all(sapply(folds, function(f) all(c("truth", "train") %in% names(f)))))

  # Check that each fold has the correct structure
  for (i in 1:3) {
    expect_true(is.matrix(folds[[i]]$truth))
    expect_true(is.matrix(folds[[i]]$train))
    expect_equal(dim(folds[[i]]$truth), dim(test_mat))
    expect_equal(dim(folds[[i]]$train), dim(test_mat))

    # Training matrix should have some NA values (holdout set)
    expect_true(sum(is.na(folds[[i]]$train)) > sum(is.na(test_mat)))
  }
})

test_that("create_cv_folds handles ground truth matrix", {
  # Create noisy input and clean ground truth
  input_mat <- matrix(runif(16, 1, 10), 4, 4)
  diag(input_mat) <- 0

  ground_truth_mat <- input_mat * 1.1  # Slightly different ground truth
  diag(ground_truth_mat) <- 0

  folds <- create_cv_folds(input_mat, ground_truth_mat, n_folds = 2, random_seed = 456)

  expect_length(folds, 2)

  # Truth should be the ground truth matrix, not the input
  for (i in 1:2) {
    expect_equal(folds[[i]]$truth, ground_truth_mat)
    # Training should be based on input matrix
    expect_true(sum(is.na(folds[[i]]$train)) > sum(is.na(input_mat)))
  }
})

test_that("create_cv_folds maintains matrix symmetry", {
  test_mat <- matrix(runif(16), 4, 4)
  test_mat[lower.tri(test_mat)] <- t(test_mat)[lower.tri(test_mat)]  # Make symmetric
  diag(test_mat) <- 0

  folds <- create_cv_folds(test_mat, n_folds = 2)

  for (i in 1:2) {
    train_mat <- folds[[i]]$train
    # Check that symmetry is maintained (if value is NA, its symmetric counterpart should also be NA)
    na_indices <- which(is.na(train_mat))
    for (idx in na_indices) {
      row <- (idx - 1) %/% nrow(train_mat) + 1
      col <- (idx - 1) %% ncol(train_mat) + 1
      if (row != col) {  # Skip diagonal
        expect_true(is.na(train_mat[col, row]))
      }
    }
  }
})

test_that("create_cv_folds input validation", {
  # Test with non-matrix input
  expect_error(create_cv_folds("not a matrix"),
               "`dissimilarity_matrix` must be a matrix")

  # Test with mismatched dimensions
  mat1 <- matrix(1:9, 3, 3)
  mat2 <- matrix(1:12, 3, 4)
  expect_error(create_cv_folds(mat1, mat2),
               "must have the same dimensions")

  # Test with invalid n_folds
  test_mat <- matrix(runif(9), 3, 3)
  expect_error(create_cv_folds(test_mat, n_folds = 1),
               "`n_folds` must be an integer greater than or equal to 2")
  expect_error(create_cv_folds(test_mat, n_folds = 5),
               "`n_folds` cannot be larger than the number of rows")

  # Test with invalid random seed
  expect_error(create_cv_folds(test_mat, random_seed = "abc"),
               "`random_seed` must be an integer")
})

test_that("create_cv_folds with reproducible results", {
  test_mat <- matrix(runif(25), 5, 5)
  diag(test_mat) <- 0

  # Test reproducibility with same seed
  folds1 <- create_cv_folds(test_mat, n_folds = 3, random_seed = 123)
  folds2 <- create_cv_folds(test_mat, n_folds = 3, random_seed = 123)

  # Should be identical
  expect_equal(folds1, folds2)

  # Test different results with different seeds
  folds3 <- create_cv_folds(test_mat, n_folds = 3, random_seed = 456)
  expect_false(identical(folds1, folds3))
})

test_that("ggsave_white_bg works as wrapper", {
  skip_if_not_installed("ggplot2")

  # Create a simple plot
  p <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1, y = 1))
  temp_file <- tempfile(fileext = ".png")

  # Test that function works without error
  expect_no_error(
    ggsave_white_bg(temp_file, p, width = 2, height = 2, dpi = 72)
  )

  # Check that file was created
  expect_true(file.exists(temp_file))

  unlink(temp_file)
})

test_that("color palette c25 is properly defined", {
  expect_true(exists("c25"))
  expect_true(is.character(c25))
  expect_length(c25, 20)  # Based on the definition in utils.R
  expect_true(all(nchar(c25) > 0))  # All elements should be non-empty strings

  # Test that colors are valid (basic check)
  expect_true(all(c25 %in% colors() | grepl("^#", c25)))
})
