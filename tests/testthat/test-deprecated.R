# tests/testthat/test-deprecated.R
# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("create_topolow_map gives deprecation warning", {
  # Create simple test matrix
  test_matrix <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow = 3)
  
  # Test that deprecation warning is issued
  expect_warning(
    result <- create_topolow_map(
      test_matrix,
      ndim = 2,
      mapping_max_iter = 10,
      k0 = 1.0,
      cooling_rate = 0.001,
      c_repulsion = 0.01,
      verbose = FALSE
    ),
    "was deprecated"
  )
  
  # Test that function still works
  expect_s3_class(result, "topolow")
  expect_true("est_distances" %in% names(result))
})

test_that("deprecated and new functions both work and return valid results", {
  test_matrix <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow = 3)
  
  # Set seed for reproducibility
  set.seed(123)
  suppressWarnings({
    result_old <- create_topolow_map(
      test_matrix,
      ndim = 2,
      mapping_max_iter = 50,
      k0 = 1.0,
      cooling_rate = 0.001,
      c_repulsion = 0.01,
      verbose = FALSE
    )
  })
  
  set.seed(123)  # Same seed for comparison
  result_new <- euclidean_embedding(
    test_matrix,
    ndim = 2,
    mapping_max_iter = 50,
    k0 = 1.0,
    cooling_rate = 0.001,
    c_repulsion = 0.01,
    verbose = FALSE
  )
  
  # Test that both functions return valid topolow objects
  expect_s3_class(result_old, "topolow")
  expect_s3_class(result_new, "topolow")
  
  # Test that both have finite, reasonable results
  expect_true(is.finite(result_old$mae))
  expect_true(is.finite(result_new$mae))
  expect_true(all(is.finite(result_old$est_distances)))
  expect_true(all(is.finite(result_new$est_distances)))
  
  # With same seed, results should be almost identical
  expect_equal(result_old$mae, result_new$mae, tolerance = 1e-2)
  expect_equal(result_old$est_distances, result_new$est_distances, tolerance = 1e-2)
})
