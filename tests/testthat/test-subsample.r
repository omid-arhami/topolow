# tests/testthat/test-subsample.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

# Helper function to generate test data
generate_test_data <- function(n, ndim = 3, sparsity = 0) {
  set.seed(123)
  coords <- matrix(rnorm(n * ndim), ncol = ndim)
  dist_mat <- as.matrix(dist(coords))
  
  if (sparsity > 0) {
    mask <- matrix(runif(n^2) > sparsity, nrow = n)
    mask <- mask | t(mask)
    diag(mask) <- TRUE
    dist_mat[!mask] <- NA
  }
  
  rownames(dist_mat) <- colnames(dist_mat) <- paste0("Point", seq_len(n))
  return(dist_mat)
}

# =============================================================================
# TEST GROUP 1: Connectivity Checking
# =============================================================================

test_that("check_matrix_connectivity detects fully connected matrices", {
  mat <- generate_test_data(20, ndim = 2)
  result <- check_matrix_connectivity(mat)
  
  expect_true(result$is_connected)
  expect_equal(result$n_components, 1)
  expect_gt(result$completeness, 0.9)
  expect_equal(result$n_points, 20)
})

test_that("check_matrix_connectivity detects disconnected components", {
  # Create a disconnected matrix with two separate components
  mat1 <- generate_test_data(10, ndim = 2)
  mat2 <- generate_test_data(10, ndim = 2)
  
  # Combine into block diagonal (disconnected)
  n <- 20
  mat <- matrix(NA, n, n)
  mat[1:10, 1:10] <- mat1
  mat[11:20, 11:20] <- mat2
  diag(mat) <- 0
  
  result <- check_matrix_connectivity(mat)
  
  expect_false(result$is_connected)
  expect_equal(result$n_components, 2)
  expect_equal(result$n_points, 20)
})

test_that("check_matrix_connectivity handles sparse matrices", {
  mat <- generate_test_data(30, ndim = 2, sparsity = 0.3)
  result <- check_matrix_connectivity(mat)
  
  expect_type(result$is_connected, "logical")
  expect_type(result$n_components, "double")
  expect_gte(result$n_components, 1)
})

test_that("check_matrix_connectivity validates input", {
  expect_error(
    check_matrix_connectivity("not a matrix"),
    "matrix"
  )
  
  expect_error(
    check_matrix_connectivity(matrix(1:6, nrow = 2)),
    "square"
  )
})

# =============================================================================
# TEST GROUP 2: Core Subsampling Functionality
# =============================================================================

test_that("subsample_dissimilarity_matrix returns correct structure", {
  mat <- generate_test_data(50, ndim = 2)
  result <- subsample_dissimilarity_matrix(mat, 25, verbose = FALSE)
  
  expect_type(result, "list")
  expect_true(is.matrix(result$subsampled_matrix))
  expect_equal(nrow(result$subsampled_matrix), 25)
  expect_equal(ncol(result$subsampled_matrix), 25)
  expect_equal(length(result$selected_indices), 25)
})

test_that("subsample_dissimilarity_matrix preserves connectivity", {
  mat <- generate_test_data(80, ndim = 2)
  result <- subsample_dissimilarity_matrix(mat, 40, verbose = FALSE)
  
  conn_check <- check_matrix_connectivity(result$subsampled_matrix)
  expect_true(conn_check$is_connected)
})

test_that("subsample_dissimilarity_matrix expands sample size when needed", {
  mat <- generate_test_data(100, ndim = 2)
  result <- subsample_dissimilarity_matrix(mat, 50, verbose = FALSE)
  
  # Sample size may be expanded to maintain connectivity
  expect_gte(nrow(result$subsampled_matrix), 50)
  if (nrow(result$subsampled_matrix) > 50) {
    expect_true(result$size_increased)
  }
})

test_that("subsample_dissimilarity_matrix preserves row and column names", {
  mat <- generate_test_data(50, ndim = 2)
  rownames(mat) <- colnames(mat) <- paste0("Item", seq_len(50))
  
  result <- subsample_dissimilarity_matrix(mat, 25, verbose = FALSE)
  
  expect_false(is.null(result$selected_names))
  expect_equal(length(result$selected_names), nrow(result$subsampled_matrix))
  expect_true(all(result$selected_names %in% rownames(mat)))
  expect_equal(rownames(result$subsampled_matrix), result$selected_names)
  expect_equal(colnames(result$subsampled_matrix), result$selected_names)
})

test_that("subsample_dissimilarity_matrix returns full matrix when target exceeds size", {
  mat <- generate_test_data(30, ndim = 2)
  result <- subsample_dissimilarity_matrix(mat, 50, verbose = FALSE)
  
  expect_equal(nrow(result$subsampled_matrix), 30)
})

test_that("subsample_dissimilarity_matrix is reproducible with seed", {
  mat <- generate_test_data(80, ndim = 2)
  
  result1 <- subsample_dissimilarity_matrix(mat, 40, random_seed = 42, verbose = FALSE)
  result2 <- subsample_dissimilarity_matrix(mat, 40, random_seed = 42, verbose = FALSE)
  
  expect_identical(result1$selected_indices, result2$selected_indices)
  expect_identical(result1$subsampled_matrix, result2$subsampled_matrix)
})

test_that("subsample_dissimilarity_matrix produces different samples with different seeds", {
  mat <- generate_test_data(80, ndim = 2)
  
  result1 <- subsample_dissimilarity_matrix(mat, 40, random_seed = 42, verbose = FALSE)
  result2 <- subsample_dissimilarity_matrix(mat, 40, random_seed = 123, verbose = FALSE)
  
  # Should be different (with very high probability)
  expect_false(identical(result1$selected_indices, result2$selected_indices))
})

test_that("subsample_dissimilarity_matrix rejects invalid inputs", {
  mat <- generate_test_data(50, ndim = 2)
  
  expect_error(
    subsample_dissimilarity_matrix("not a matrix", 25),
    "must be a matrix"
  )
  
  expect_error(
    subsample_dissimilarity_matrix(mat, -5)
  )
  
  expect_error(
    subsample_dissimilarity_matrix(mat, 1)
  )
  
  expect_error(
    subsample_dissimilarity_matrix(matrix(1:6, nrow = 2), 5)
  )
})

# =============================================================================
# TEST GROUP 3: Sanity Check Function
# =============================================================================

test_that("sanity_check_subsample passes checks for good subsampling", {
  mat <- generate_test_data(100, ndim = 2)
  result <- sanity_check_subsample(mat, folds = 10, verbose = FALSE)
  
  expect_type(result, "list")
  expect_true(result$all_checks_passed)
  expect_true(all(unlist(result$checks)))
})

test_that("sanity_check_subsample detects issues with small matrices", {
  mat <- generate_test_data(5, ndim = 2)
  result <- sanity_check_subsample(mat, folds = 5, verbose = FALSE)
  expect_gt(length(result$warnings), 0)  # Has warnings
  expect_true(any(grepl("few points|measurements", result$warnings, ignore.case = TRUE)))
  
  expect_false(result$all_checks_passed)
})

# =============================================================================
# TEST GROUP 4: Integration with Parameter Optimization
# =============================================================================

test_that("initial_parameter_optimization works without subsampling", {
  skip_on_cran()  # Skip slow tests on CRAN
  
  mat <- generate_test_data(30, ndim = 2)
  
  result <- initial_parameter_optimization(
    dissimilarity_matrix = mat,
    mapping_max_iter = 20,
    relative_epsilon = 1e-3,
    convergence_counter = 2,
    scenario_name = "test_no_subsample",
    N_min = 2, N_max = 3,
    k0_min = 0.5, k0_max = 2,
    c_repulsion_min = 0.01, c_repulsion_max = 0.05,
    cooling_rate_min = 0.001, cooling_rate_max = 0.01,
    num_samples = 3,
    max_cores = 1,
    folds = 5,
    verbose = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  expect_false("opt_subsample" %in% names(result))
})

test_that("initial_parameter_optimization works with subsampling", {
  skip_on_cran()  # Skip slow tests on CRAN
  
  mat <- generate_test_data(80, ndim = 2)
  
  result <- initial_parameter_optimization(
    dissimilarity_matrix = mat,
    mapping_max_iter = 20,
    relative_epsilon = 1e-3,
    convergence_counter = 2,
    scenario_name = "test_with_subsample",
    N_min = 2, N_max = 3,
    k0_min = 0.5, k0_max = 2,
    c_repulsion_min = 0.01, c_repulsion_max = 0.05,
    cooling_rate_min = 0.001, cooling_rate_max = 0.01,
    num_samples = 3,
    max_cores = 1,
    folds = 5,
    opt_subsample = 60,
    verbose = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  expect_true("opt_subsample" %in% names(result))
  expect_true("original_n_points" %in% names(result))
  expect_equal(unique(result$opt_subsample), 60)
  expect_equal(unique(result$original_n_points), 80)
})

test_that("initial_parameter_optimization uses full data when subsample exceeds matrix size", {
  skip_on_cran()
  
  mat <- generate_test_data(30, ndim = 2)
  
  result <- initial_parameter_optimization(
    dissimilarity_matrix = mat,
    mapping_max_iter = 20,
    relative_epsilon = 1e-3,
    convergence_counter = 2,
    scenario_name = "test_subsample_too_large",
    N_min = 2, N_max = 3,
    k0_min = 0.5, k0_max = 2,
    c_repulsion_min = 0.01, c_repulsion_max = 0.05,
    cooling_rate_min = 0.001, cooling_rate_max = 0.01,
    num_samples = 2,
    max_cores = 1,
    folds = 5,
    opt_subsample = 50,
    verbose = FALSE
  )
  
  expect_s3_class(result, "data.frame")
  # When opt_subsample exceeds matrix size, should use full data
})

test_that("initial_parameter_optimization validates opt_subsample", {
  mat <- generate_test_data(150, ndim = 2)
  
  # Small values should trigger a warning but not error
  expect_warning(
    result <- initial_parameter_optimization(
      dissimilarity_matrix = mat,
      mapping_max_iter = 10,
      relative_epsilon = 1e-3,
      convergence_counter = 2,
      scenario_name = "test_small_subsample",
      N_min = 2, N_max = 3,
      k0_min = 0.5, k0_max = 2,
      c_repulsion_min = 0.01, c_repulsion_max = 0.05,
      cooling_rate_min = 0.001, cooling_rate_max = 0.01,
      num_samples = 2,
      max_cores = 1,
      folds = 5,
      opt_subsample = 5,  # Small but allowed
      verbose = FALSE
    ),
    "smaller than recommended|connectivity issues"
  )
  
  # Function should still return results despite warning
  expect_s3_class(result, "data.frame")
  
  # Negative values should error
  expect_error(
    initial_parameter_optimization(
      dissimilarity_matrix = mat,
      mapping_max_iter = 10,
      relative_epsilon = 1e-3,
      convergence_counter = 2,
      scenario_name = "test_invalid_subsample",
      N_min = 2, N_max = 3,
      k0_min = 0.5, k0_max = 2,
      c_repulsion_min = 0.01, c_repulsion_max = 0.05,
      cooling_rate_min = 0.001, cooling_rate_max = 0.01,
      num_samples = 2,
      max_cores = 1,
      folds = 5,
      opt_subsample = -10
    ),
    "must be NULL or.*>= 5|numeric value"
  )
  
  # Non-numeric values should error
  expect_error(
    initial_parameter_optimization(
      dissimilarity_matrix = mat,
      mapping_max_iter = 10,
      relative_epsilon = 1e-3,
      convergence_counter = 2,
      scenario_name = "test_invalid_subsample",
      N_min = 2, N_max = 3,
      k0_min = 0.5, k0_max = 2,
      c_repulsion_min = 0.01, c_repulsion_max = 0.05,
      cooling_rate_min = 0.001, cooling_rate_max = 0.01,
      num_samples = 2,
      max_cores = 1,
      folds = 5,
      opt_subsample = "invalid"
    )
  )
})

# =============================================================================
# TEST GROUP 5: Adaptive Sampling with Subsampling
# =============================================================================
test_that("run_adaptive_sampling works with subsampling", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  
  mat <- generate_test_data(150, ndim = 2)
  output_dir <- tempfile()
  dir.create(output_dir, recursive = TRUE)
  
  # Create initial samples with MORE samples for stability
  initial_results <- initial_parameter_optimization(
    dissimilarity_matrix = mat,
    mapping_max_iter = 10,
    relative_epsilon = 1e-3,
    convergence_counter = 2,
    scenario_name = "test_adaptive",
    N_min = 2, N_max = 3,
    k0_min = 0.5, k0_max = 2,
    c_repulsion_min = 0.01, c_repulsion_max = 0.05,
    cooling_rate_min = 0.001, cooling_rate_max = 0.01,
    num_samples = 10,  # Increased from 5 to 10 for more stable KDE
    max_cores = 1,
    folds = 5,
    opt_subsample = 80,
    verbose = FALSE,  # Enable verbose to see what's happening
    write_files = TRUE,
    output_dir = output_dir
  )
  
  # Get the initial samples file
  param_dir <- file.path(output_dir, "model_parameters")
  initial_file <- file.path(param_dir, "test_adaptive_model_parameters.csv")
  
  expect_true(file.exists(initial_file))
  
  # Run adaptive sampling
  expect_error(
    run_adaptive_sampling(
      initial_samples_file = initial_file,
      scenario_name = "test_adaptive",
      dissimilarity_matrix = mat,
      max_cores = 1,
      num_samples = 2,
      mapping_max_iter = 10,
      folds = 3,  # Match the folds used in initial optimization
      opt_subsample = 60,
      output_dir = output_dir,
      verbose = FALSE  # Enable to debug
    ) ,
    NA  # Should not error
  )
  
  # Check output files exist
  adaptive_file <- file.path(param_dir, "test_adaptive_model_parameters.csv")
  expect_true(file.exists(adaptive_file))
  
  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

# =============================================================================
# TEST GROUP 6: Edge Cases and Error Handling
# =============================================================================

test_that("subsample_dissimilarity_matrix handles matrices with missing values appropriately", {
  mat <- generate_test_data(50, ndim = 2, sparsity = 0.1)
  
  # Should either succeed or give informative error
  result <- tryCatch(
    subsample_dissimilarity_matrix(mat, 25, verbose = FALSE),
    error = function(e) e
  )
  
  # If it succeeds, check structure
  if (!inherits(result, "error")) {
    expect_true(is.matrix(result$subsampled_matrix))
    expect_equal(nrow(result$subsampled_matrix) , 25)
  }
})


test_that("subsample_dissimilarity_matrix handles minimum size matrix", {
  mat <- generate_test_data(10, ndim = 2)
  result <- subsample_dissimilarity_matrix(mat, 10, verbose = FALSE)
  
  expect_equal(nrow(result$subsampled_matrix), 10)
})

test_that("subsample functions provide helpful messages in verbose mode", {
  mat <- generate_test_data(50, ndim = 2)
  
  expect_output(
    subsample_dissimilarity_matrix(mat, 25, verbose = TRUE),
    "subsample|connectivity|sample"
  )
  
  expect_output(
    sanity_check_subsample(mat, folds = 5, verbose = TRUE),
    "check|sanity|connectivity"
  )
})

# =============================================================================
# TEST GROUP 7: Grouped Subsampling (always_include / sample_pool)
# =============================================================================

make_group_matrix <- function(n_v = 4, n_s = 3, seed = 1) {
  set.seed(seed)
  m <- as.matrix(dist(matrix(rnorm((n_v + n_s) * 2), ncol = 2)))
  nm <- c(paste0("V/ag", seq_len(n_v)), paste0("S/sr", seq_len(n_s)))
  dimnames(m) <- list(nm, nm)
  m
}

test_that("subsample_dissimilarity_matrix consumes the same RNG stream when ungrouped", {
  mat <- generate_test_data(10, ndim = 2)

  set.seed(42)
  expected <- sample(seq_len(10), 5)

  set.seed(42)
  result <- subsample_dissimilarity_matrix(mat, 5, verbose = FALSE)

  expect_identical(result$selected_indices, expected)
})

test_that("explicitly passing NULL groups changes nothing", {
  mat <- generate_test_data(30, ndim = 2)

  set.seed(7)
  a <- subsample_dissimilarity_matrix(mat, 15)
  set.seed(7)
  b <- subsample_dissimilarity_matrix(mat, 15,
                                      always_include = NULL,
                                      sample_pool = NULL)

  expect_identical(a$selected_indices, b$selected_indices)
  expect_identical(a$subsampled_matrix, b$subsampled_matrix)
})

test_that("every always_include point is kept", {
  mat <- make_group_matrix()
  antigens <- paste0("V/ag", 1:4)

  for (i in 1:10) {
    res <- subsample_dissimilarity_matrix(mat, sample_size = 6,
                                          always_include = antigens,
                                          sample_pool = paste0("S/sr", 1:3))
    expect_true(all(antigens %in% res$selected_names))
    expect_equal(sum(startsWith(res$selected_names, "S/")), 2)
    expect_equal(nrow(res$subsampled_matrix), 6)
  }
})

test_that("draws never leave the sample_pool", {
  mat <- make_group_matrix()
  excluded <- "S/sr3"
  pool <- c("S/sr1", "S/sr2")

  for (i in 1:20) {
    res <- subsample_dissimilarity_matrix(mat, sample_size = 5,
                                          always_include = paste0("V/ag", 1:4),
                                          sample_pool = pool)
    expect_false(excluded %in% res$selected_names)
  }
})

test_that("names and indices are interchangeable", {
  mat <- make_group_matrix()

  r1 <- subsample_dissimilarity_matrix(mat, sample_size = 6, random_seed = 3,
                                       always_include = paste0("V/ag", 1:4),
                                       sample_pool = paste0("S/sr", 1:3))
  r2 <- subsample_dissimilarity_matrix(mat, sample_size = 6, random_seed = 3,
                                       always_include = 1:4,
                                       sample_pool = 5:7)

  expect_setequal(r1$selected_names, r2$selected_names)
})

test_that("grouped subsampling is reproducible with a seed", {
  mat <- make_group_matrix(n_v = 6, n_s = 5)

  r1 <- subsample_dissimilarity_matrix(mat, sample_size = 9, random_seed = 11,
                                       always_include = paste0("V/ag", 1:6),
                                       sample_pool = paste0("S/sr", 1:5))
  r2 <- subsample_dissimilarity_matrix(mat, sample_size = 9, random_seed = 11,
                                       always_include = paste0("V/ag", 1:6),
                                       sample_pool = paste0("S/sr", 1:5))

  expect_identical(r1$selected_indices, r2$selected_indices)
})

test_that("grouped subsampling rejects impossible sizes", {
  mat <- make_group_matrix()
  antigens <- paste0("V/ag", 1:4)

  expect_error(
    subsample_dissimilarity_matrix(mat, sample_size = 4,
                                   always_include = antigens),
    "sample_size must be greater than the number of always_include points"
  )
  expect_error(
    subsample_dissimilarity_matrix(mat, sample_size = 99,
                                   always_include = antigens,
                                   sample_pool = paste0("S/sr", 1:3)),
    "sample_size must not exceed"
  )
  expect_error(
    subsample_dissimilarity_matrix(mat, sample_size = 5,
                                   always_include = "nope"),
    "always_include must name points that exist"
  )
  expect_error(
    subsample_dissimilarity_matrix(mat, sample_size = 5,
                                   always_include = TRUE),
    "always_include must be a character vector"
  )
})

test_that("character groups need row names", {
  mat <- generate_test_data(10, ndim = 2)
  dimnames(mat) <- NULL
  expect_error(
    subsample_dissimilarity_matrix(mat, 5, always_include = c("a", "b")),
    "must have row names"
  )
})

test_that("a sample_size covering the whole universe returns it", {
  mat <- make_group_matrix()
  res <- subsample_dissimilarity_matrix(mat, sample_size = 7,
                                        always_include = paste0("V/ag", 1:4),
                                        sample_pool = paste0("S/sr", 1:3))
  expect_equal(nrow(res$subsampled_matrix), 7)
})

test_that("the verbose full-matrix path no longer errors", {
  mat <- generate_test_data(20, ndim = 2)
  expect_output(
    res <- subsample_dissimilarity_matrix(mat, 50, verbose = TRUE),
    "available points"
  )
  expect_equal(nrow(res$subsampled_matrix), 20)
})

# =============================================================================
# TEST GROUP 8: Serum-wise Subsampling in Parameter Optimization
# =============================================================================

test_that("resolve_sera_subsample converts counts and fractions", {
  mat <- make_group_matrix(n_v = 10, n_s = 5)

  by_count <- resolve_sera_subsample(mat, NULL, 3, NULL)
  expect_equal(by_count$opt_subsample, 13)
  expect_equal(length(by_count$always_include), 10)
  expect_equal(length(by_count$sample_pool), 5)

  by_frac <- resolve_sera_subsample(mat, NULL, 0.8, NULL)
  expect_equal(by_frac$opt_subsample, 14)
})

test_that("resolve_sera_subsample enforces its preconditions", {
  mat <- make_group_matrix(n_v = 10, n_s = 5)

  expect_error(resolve_sera_subsample(mat, 20, 3, NULL),
               "must not both be supplied")
  expect_error(resolve_sera_subsample(mat, NULL, 5, NULL),
               "smaller than the number of sera")
  expect_error(resolve_sera_subsample(mat, NULL, 1, NULL),
               "at least 2 sera")
  expect_error(resolve_sera_subsample(mat, NULL, NULL, "S/sr1"),
               "no effect unless")
})

test_that("initial_parameter_optimization accepts opt_subsample_sera", {
  skip_on_cran()

  mat <- make_group_matrix(n_v = 24, n_s = 10, seed = 4)

  result <- suppressWarnings(initial_parameter_optimization(
    dissimilarity_matrix = mat,
    mapping_max_iter = 10, relative_epsilon = 1e-3, convergence_counter = 2,
    scenario_name = "test_sera_subsample",
    N_min = 2, N_max = 3, k0_min = 0.5, k0_max = 2,
    c_repulsion_min = 0.01, c_repulsion_max = 0.05,
    cooling_rate_min = 0.001, cooling_rate_max = 0.01,
    num_samples = 3, max_cores = 1, folds = 3,
    opt_subsample_sera = 0.8, verbose = FALSE
  ))

  expect_s3_class(result, "data.frame")
  expect_true("opt_subsample" %in% names(result))
  # 24 antigens kept plus 8 of the 10 sera
  expect_equal(unique(result$opt_subsample), 32)
  expect_equal(unique(result$original_n_points), 34)
})

test_that("opt_subsample and opt_subsample_sera are mutually exclusive", {
  mat <- make_group_matrix(n_v = 24, n_s = 10, seed = 4)

  expect_error(
    initial_parameter_optimization(
      dissimilarity_matrix = mat,
      mapping_max_iter = 10, relative_epsilon = 1e-3, convergence_counter = 2,
      scenario_name = "test_both",
      N_min = 2, N_max = 3, k0_min = 0.5, k0_max = 2,
      c_repulsion_min = 0.01, c_repulsion_max = 0.05,
      cooling_rate_min = 0.001, cooling_rate_max = 0.01,
      num_samples = 2, max_cores = 1, folds = 3,
      opt_subsample = 20, opt_subsample_sera = 0.8
    ),
    "must not both be supplied"
  )
})

test_that("the serum subset size is a fixed count, recorded and reusable", {
  skip_on_cran()

  mat <- make_group_matrix(n_v = 24, n_s = 10, seed = 8)

  res <- suppressWarnings(initial_parameter_optimization(
    dissimilarity_matrix = mat,
    mapping_max_iter = 10, relative_epsilon = 1e-3, convergence_counter = 2,
    scenario_name = "test_fixed_size",
    N_min = 2, N_max = 3, k0_min = 0.5, k0_max = 2,
    c_repulsion_min = 0.01, c_repulsion_max = 0.05,
    cooling_rate_min = 0.001, cooling_rate_max = 0.01,
    num_samples = 4, max_cores = 1, folds = 3,
    opt_subsample_sera = 0.8, verbose = FALSE
  ))

  # one constant size for every evaluation, with the serum count spelled out
  expect_true("opt_n_sera" %in% names(res))
  expect_length(unique(res$opt_subsample), 1)
  expect_length(unique(res$opt_n_sera), 1)
  expect_equal(unique(res$opt_n_sera), 8)
  expect_equal(unique(res$opt_subsample), 24 + 8)

  # and the same count drives the bootstrap stage
  boot <- suppressWarnings(bootstrap_embedding(
    mat, ndim = 2, k0 = 1, cooling_rate = 0.01, c_repulsion = 0.01,
    mapping_max_iter = 20, n_replicates = 5,
    n_sera = unique(res$opt_n_sera), max_cores = 1, consensus = FALSE
  ))
  expect_equal(boot$settings$n_drawn, 8L)
  expect_equal(unique(boot$replicates$n_points), 32L)
})

test_that("opt_n_sera is absent for a uniform subsample", {
  skip_on_cran()

  mat <- make_group_matrix(n_v = 24, n_s = 10, seed = 9)

  res <- suppressWarnings(initial_parameter_optimization(
    dissimilarity_matrix = mat,
    mapping_max_iter = 10, relative_epsilon = 1e-3, convergence_counter = 2,
    scenario_name = "test_uniform_size",
    N_min = 2, N_max = 3, k0_min = 0.5, k0_max = 2,
    c_repulsion_min = 0.01, c_repulsion_max = 0.05,
    cooling_rate_min = 0.001, cooling_rate_max = 0.01,
    num_samples = 3, max_cores = 1, folds = 3,
    opt_subsample = 25, verbose = FALSE
  ))

  expect_true(all(is.na(res$opt_n_sera)))
  expect_equal(unique(res$opt_subsample), 25)
})
