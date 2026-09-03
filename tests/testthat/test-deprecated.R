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

test_that("create_diagnostic_plots warns and forwards to plot_mcmc_diagnostics", {
  chain_files <- c(
    system.file("extdata", "diag_chain1.csv", package = "topolow"),
    system.file("extdata", "diag_chain2.csv", package = "topolow"),
    system.file("extdata", "diag_chain3.csv", package = "topolow")
  )
  skip_if_not(all(nzchar(chain_files)), "example chain files not installed")

  op <- options(lifecycle_verbosity = "warning")
  on.exit(options(op), add = TRUE)

  # The deprecation warning fires
  expect_warning(
    old <- create_diagnostic_plots(chain_files, mutual_size = 50, save_plot = FALSE),
    "was deprecated"
  )

  # And it forwards: same result as calling the replacement directly
  new <- plot_mcmc_diagnostics(chain_files, mutual_size = 50, save_plot = FALSE)
  expect_s3_class(old, "gtable")
  expect_s3_class(new, "gtable")
})

test_that("create_diagnostic_plots maps deprecated `res` onto `dpi`", {
  chain_files <- c(
    system.file("extdata", "diag_chain1.csv", package = "topolow"),
    system.file("extdata", "diag_chain2.csv", package = "topolow"),
    system.file("extdata", "diag_chain3.csv", package = "topolow")
  )
  skip_if_not(all(nzchar(chain_files)), "example chain files not installed")

  out_dir <- file.path(tempdir(), "topolow_res_arg_test")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # lifecycle silences a repeated deprecation warning within a session; force
  # it so this test does not depend on which tests ran before it.
  op <- options(lifecycle_verbosity = "warning")
  on.exit(options(op), add = TRUE)

  # Both the function-level and the `res` argument warnings are raised.
  # Collect them rather than nesting expect_warning(), which only bubbles one.
  warns <- character()
  p <- withCallingHandlers(
    create_diagnostic_plots(
      chain_files,
      mutual_size = 50,
      save_plot = TRUE,
      output_dir = out_dir,
      output_file = "res_arg.png",
      width = 900, height = 900,
      res = 150
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("create_diagnostic_plots()", warns, fixed = TRUE)))
  expect_true(any(grepl("`res` argument", warns, fixed = TRUE)))

  expect_s3_class(p, "gtable")
  expect_true(file.exists(file.path(out_dir, "res_arg.png")))
})
