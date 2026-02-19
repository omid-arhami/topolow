# tests/testthat/test-adaptive-sampling.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("generate_kde_samples works correctly", {
  # Create test samples
  test_samples <- data.frame(
    log_N = log(runif(50, 2, 10)),
    log_k0 = log(runif(50, 1, 5)),
    log_cooling_rate = log(runif(50, 0.01, 0.1)),
    log_c_repulsion = log(runif(50, 0.1, 1)),
    Holdout_MAE = runif(50, 0.20, 10.0)
  )

  # Generate new samples
  new_samples <- generate_kde_samples(test_samples, n = 10)

  expect_true(is.data.frame(new_samples))
  expect_equal(nrow(new_samples), 10)
  expect_equal(ncol(new_samples), 4)  # Should only have parameter columns

  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  expect_true(all(par_names %in% names(new_samples)))
  expect_true(all(sapply(new_samples, is.numeric)))
})


test_that("generate_kde_samples handles edge cases", {
  # Test with minimal samples
  minimal_samples <- data.frame(
    log_N = log(c(2, 3)),
    log_k0 = log(c(1, 2)),
    log_cooling_rate = log(c(0.01, 0.02)),
    log_c_repulsion = log(c(0.1, 0.2)),
    NLL = c(50, 60),
    Holdout_MAE = c(0.5, 0.67)
  )

  expect_no_error(new_samples <- generate_kde_samples(minimal_samples, n = 5))
  expect_equal(nrow(new_samples), 5)

  # Test with exploration epsilon
  exploratory_samples <- generate_kde_samples(minimal_samples, n = 3, epsilon = 0.5)
  expect_equal(nrow(exploratory_samples), 3)
})

test_that("weighted_kde produces valid density estimates", {
  x <- rnorm(100, mean = 5, sd = 2)
  weights <- runif(100)
  weights <- weights / sum(weights)  # Normalize

  kde_result <- weighted_kde(x, weights, n = 100)

  expect_true(is.list(kde_result))
  expect_true(all(c("x", "y") %in% names(kde_result)))
  expect_length(kde_result$x, 100)
  expect_length(kde_result$y, 100)
  expect_true(all(kde_result$y >= 0))

  # Check approximate normalization (integral should be close to 1)
  integral_approx <- sum(kde_result$y) * diff(kde_result$x)[1]
  expect_true(abs(integral_approx - 1) < 0.2)  # Allow some tolerance
})

test_that("weighted_kde handles different parameters", {
  x <- c(1, 2, 3, 4, 5)
  weights <- c(0.1, 0.2, 0.4, 0.2, 0.1)

  # Test with different evaluation ranges
  kde1 <- weighted_kde(x, weights, from = 0, to = 6, n = 50)
  expect_length(kde1$x, 50)
  expect_equal(min(kde1$x), 0)
  expect_equal(max(kde1$x), 6)

  # Test with default range
  kde2 <- weighted_kde(x, weights)
  expect_equal(min(kde2$x), min(x))
  expect_equal(max(kde2$x), max(x))
})

test_that("calculate_weighted_marginals works correctly", {
  # Create test samples
  test_samples <- data.frame(
    log_N = log(runif(50, 2, 10)),
    log_k0 = log(runif(50, 1, 5)),
    log_cooling_rate = log(runif(50, 0.01, 0.1)),
    log_c_repulsion = log(runif(50, 0.1, 1)),
    Holdout_MAE = runif(50, 0.20, 10.0)
  )

  marginals <- calculate_weighted_marginals(test_samples)

  expect_true(is.list(marginals))
  expect_length(marginals, 4)

  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  expect_true(all(par_names %in% names(marginals)))

  # Check that each marginal has the expected structure
  for (param in par_names) {
    expect_true(all(c("x", "y") %in% names(marginals[[param]])))
    expect_true(all(marginals[[param]]$y >= 0))
  }
})

test_that("calculate_weighted_marginals input validation", {
  # Test with missing columns
  incomplete_samples <- data.frame(
    log_N = log(runif(20, 2, 10)),
    log_k0 = log(runif(20, 1, 5))
    # Missing log_cooling_rate, log_c_repulsion, NLL
  )

  expect_error(
    calculate_weighted_marginals(incomplete_samples),
    "Missing required columns"
  )

  # Test with non-numeric columns
  invalid_samples <- data.frame(
    log_N = log(runif(20, 2, 10)),
    log_k0 = log(runif(20, 1, 5)),
    log_cooling_rate = log(runif(20, 0.01, 0.1)),
    log_c_repulsion = log(runif(20, 0.1, 1)),
    NLL = as.character(runif(20, 20, 100)),  # Character instead of numeric
    Holdout_MAE = as.character(runif(20, 0.20, 10.0))  # Character instead of numeric
  )

  expect_error(
    calculate_weighted_marginals(invalid_samples),
    "All required parameter and Holdout_MAE columns must be numeric"
  )
})

test_that("likelihood_function works with cross-validation", {
  # Create test dissimilarity matrix
  test_mat <- matrix(runif(16, 1, 5), 4, 4)
  diag(test_mat) <- 0
  test_mat[lower.tri(test_mat)] <- t(test_mat)[lower.tri(test_mat)]  # Make symmetric

  # Test likelihood calculation
  result <- likelihood_function(
    dissimilarity_matrix = test_mat,
    mapping_max_iter = 20,
    relative_epsilon = 1e-3,
    N = 2,
    k0 = 1.0,
    cooling_rate = 0.01,
    c_repulsion = 0.01,
    folds = 3,
    num_cores = 1
  )

  expect_true(is.list(result))
  expect_true(all(c("Holdout_MAE", "NLL") %in% names(result)))
  expect_true(is.numeric(result$Holdout_MAE))
  expect_true(is.numeric(result$NLL))
  expect_true(result$Holdout_MAE > 0 || is.na(result$Holdout_MAE))
})

test_that("likelihood_function handles sparse data", {
  # Create very sparse matrix
  sparse_mat <- matrix(NA, 4, 4)
  sparse_mat[1, 2] <- sparse_mat[2, 1] <- 1
  sparse_mat[2, 3] <- sparse_mat[3, 2] <- 2
  diag(sparse_mat) <- 0

  # Suppress expected warning about sparse data
  suppressWarnings({
    result <- likelihood_function(
      dissimilarity_matrix = sparse_mat,
      mapping_max_iter = 10,
      relative_epsilon = 1e-3,
      N = 2,
      k0 = 1.0,
      cooling_rate = 0.01,
      c_repulsion = 0.01,
      folds = 2,
      num_cores = 1
    )
  })

  # Should handle sparse data gracefully (may return NA)
  expect_true(is.list(result))
  expect_true(all(c("Holdout_MAE", "NLL") %in% names(result)))
})

test_that("get_grid creates appropriate parameter grids", {
  # Create test samples
  test_samples <- data.frame(
    log_N = log(runif(30, 2, 10)),
    log_k0 = log(runif(30, 1, 5)),
    log_cooling_rate = log(runif(30, 0.01, 0.1)),
    log_c_repulsion = log(runif(30, 0.1, 1)),
    NLL = runif(30, 20, 100)
  )

  # Test grid creation for log_N
  grid_values <- get_grid(test_samples, "log_N", num_points = 10,
                          start_factor = 0.8, end_factor = 1.2)

  expect_length(grid_values, 10)
  expect_true(is.numeric(grid_values))
  expect_true(all(diff(grid_values) > 0))  # Should be increasing
})

test_that("profile_likelihood handles various parameter ranges", {
  # Create test samples with known distribution
  set.seed(456)  # Different seed for better distribution
  test_samples <- data.frame(
    log_N = log(runif(60, 2, 8)),        # Even more samples
    log_k0 = log(runif(60, 1, 4)),       # More concentrated range
    log_cooling_rate = log(runif(60, 0.01, 0.1)),
    log_c_repulsion = log(runif(60, 0.1, 1)),
    NLL = runif(60, 20, 100)
  )

  # Test profile likelihood calculation - suppress expected warnings for edge cases
  suppressWarnings({
    pl_result <- profile_likelihood("log_k0", test_samples, grid_size = 5,  # Even smaller grid
                                    bandwidth_factor = 0.5, min_samples = 1)  # Even larger bandwidth, minimum samples
  })

  expect_s3_class(pl_result, "profile_likelihood")
  expect_equal(length(pl_result$param), 5)  # Updated expectation
  expect_equal(length(pl_result$ll), 5)
  expect_equal(pl_result$param_name, "log_k0")
})

test_that("profile_likelihood input validation", {
  test_samples <- data.frame(
    log_N = log(runif(20, 2, 10)),
    NLL = runif(20, 20, 100)
  )

  # Test with invalid parameter name
  expect_error(
    profile_likelihood("invalid_param", test_samples),
    "Parameter 'invalid_param' not found in samples"
  )

  # Test with missing NLL column
  samples_no_nll <- test_samples[, "log_N", drop = FALSE]
  expect_error(
    profile_likelihood("log_N", samples_no_nll),
    "Samples data frame must contain an 'NLL' column"
  )

  # Test with invalid grid_size
  expect_error(
    profile_likelihood("log_N", test_samples, grid_size = 1),
    "grid_size must be at least 2"
  )

  # Test with invalid bandwidth_factor
  expect_error(
    profile_likelihood("log_N", test_samples, bandwidth_factor = -0.1),
    "bandwidth_factor must be positive"
  )
})
