# tests/testthat/test-data-preprocessing.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("titers_list_to_matrix creates correct matrix", {
  # Based on the original failing test data
  test_data <- data.frame(
    object = c("A1", "A2", "S1"),
    reference = c("S1", "S1", "A2"),
    dissimilarity = c(1.0, 2.0, 3.0)
  )

  result <- titers_list_to_matrix(
    test_data,
    chnames = "object",
    rnames = "reference", 
    values_column = "dissimilarity"
    # rc = FALSE is the default, so V/ and S/ prefixes will be added
  )

  expect_true(is.matrix(result))
  # With rc=FALSE, we get prefixed names creating 5 unique points:
  # V/A1, V/A2, V/S1, S/S1, S/A2
  expect_equal(dim(result), c(5, 5))
  expect_true(isSymmetric(result))
  
  # Check values with prefixed names
  expect_equal(result["V/A1", "S/S1"], 1.0)
  expect_equal(result["V/A2", "S/S1"], 2.0)  
  expect_equal(result["V/S1", "S/A2"], 3.0)
  
  # Check symmetry
  expect_equal(result["S/S1", "V/A1"], 1.0)
  expect_equal(result["S/S1", "V/A2"], 2.0)
  expect_equal(result["S/A2", "V/S1"], 3.0)
})


test_that("clean_data removes outliers", {
  x <- c(1, 1.1, 1.2, 10, 0.9, 0.8, -5)
  cleaned_x <- clean_data(x, k = 2)
  # Expect outliers (10 and -5) to be replaced by NA
  expect_true(is.na(cleaned_x[4]))
  expect_true(is.na(cleaned_x[7]))
  expect_equal(sum(is.na(cleaned_x)), 2)
})