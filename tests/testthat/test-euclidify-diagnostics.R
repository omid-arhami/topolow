# tests/testthat/test-euclidify-diagnostics.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("Euclidify creates diagnostic plots when requested", {
  dist_mat <- matrix(runif(400), 20, 20)
  temp_dir <- tempfile()
  
  result <- Euclidify(
    dissimilarity_matrix = dist_mat,
    output_dir = temp_dir,
    ndim_range = c(2, 2),
    n_initial_samples = 10,
    n_adaptive_samples = 10,
    max_cores = 1,
    folds = 5,
    mapping_max_iter = 10,
    create_diagnostic_plots = TRUE,
    verbose = "off"
  )
  
  # Check plots created
  diag_dir <- file.path(temp_dir, "diagnostics")
  expect_true(dir.exists(diag_dir))
  expect_true(file.exists(file.path(diag_dir, "parameter_search_trajectory.png")))
  expect_true(file.exists(file.path(diag_dir, "embedding_quality.png")))
  expect_true(file.exists(file.path(diag_dir, "diagnostic_report.txt")))
  
  # Check result structure
  expect_true("all_samples" %in% names(result))
  expect_true("diagnostic_plots" %in% names(result))
  
  unlink(temp_dir, recursive = TRUE)
})
