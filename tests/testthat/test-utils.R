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
  expect_true(all(nchar(c25) > 0))  # All elements should be non-empty strings

  # Test that colors are valid (basic check)
  expect_true(all(c25 %in% colors() | grepl("^#", c25)))
})

# =============================================================================
# identify_antigens_sera
# =============================================================================

make_vs_mat <- function(nms) {
  m <- as.matrix(dist(matrix(rnorm(length(nms) * 2), ncol = 2)))
  dimnames(m) <- list(nms, nms)
  m
}

test_that("identify_antigens_sera splits on the V/ and S/ convention", {
  set.seed(1)
  nms <- c("V/ag1", "V/ag2", "V/ag3", "V/ag4", "S/sr1", "S/sr2", "S/sr3")
  m <- make_vs_mat(nms)

  g <- identify_antigens_sera(m)

  expect_equal(g$antigens, nms[1:4])
  expect_equal(g$sera, nms[5:7])
  expect_equal(g$antigen_indices, 1:4)
  expect_equal(g$serum_indices, 5:7)
  expect_length(g$unassigned, 0)
  expect_equal(rownames(m)[g$antigen_indices], g$antigens)
})

test_that("identify_antigens_sera infers the other group from an explicit one", {
  set.seed(2)
  nms <- c("a1", "a2", "a3", "s1", "s2")
  m <- make_vs_mat(nms)

  g <- identify_antigens_sera(m, sera = c("s1", "s2"))
  expect_equal(g$sera, c("s1", "s2"))
  expect_equal(g$antigens, c("a1", "a2", "a3"))

  g2 <- identify_antigens_sera(m, antigens = c("a1", "a2", "a3"))
  expect_equal(g2$sera, c("s1", "s2"))
})

test_that("identify_antigens_sera accepts custom patterns", {
  set.seed(3)
  nms <- c("ag1", "ag2", "ag3", "sr1", "sr2")
  m <- make_vs_mat(nms)

  g <- identify_antigens_sera(m, antigen_pattern = "^ag", serum_pattern = "^sr")
  expect_equal(g$antigens, c("ag1", "ag2", "ag3"))
  expect_equal(g$sera, c("sr1", "sr2"))
})

test_that("identify_antigens_sera reports unassigned points", {
  set.seed(4)
  nms <- c("V/ag1", "V/ag2", "S/sr1", "X/other")
  m <- make_vs_mat(nms)

  g <- identify_antigens_sera(m)
  expect_equal(g$unassigned, "X/other")
})

test_that("identify_antigens_sera errors when nothing matches", {
  set.seed(5)
  nms <- paste0("Point", 1:5)
  m <- make_vs_mat(nms)

  expect_error(
    identify_antigens_sera(m),
    "antigen_pattern or serum_pattern must match some row names"
  )
})

test_that("identify_antigens_sera rejects overlapping groups and bad input", {
  set.seed(6)
  nms <- c("a1", "a2", "s1", "s2")
  m <- make_vs_mat(nms)

  expect_error(
    identify_antigens_sera(m, antigens = c("a1", "s1"), sera = c("s1", "s2")),
    "must not share any point"
  )
  expect_error(identify_antigens_sera("not a matrix"),
               "dissimilarity_matrix must be a matrix")

  unnamed <- m
  dimnames(unnamed) <- NULL
  expect_error(identify_antigens_sera(unnamed),
               "must have row names")

  dup <- m
  rownames(dup)[2] <- rownames(dup)[1]
  colnames(dup) <- rownames(dup)
  expect_error(identify_antigens_sera(dup), "unique row names")
})

# =============================================================================
# resolve_point_indices
# =============================================================================

test_that("resolve_point_indices handles names, indices and NULL", {
  set.seed(7)
  nms <- c("p1", "p2", "p3", "p4")
  m <- make_vs_mat(nms)

  expect_null(resolve_point_indices(NULL, m, "x"))
  expect_equal(resolve_point_indices(c("p2", "p4"), m, "x"), c(2L, 4L))
  expect_equal(resolve_point_indices(c(2, 4), m, "x"), c(2L, 4L))
  expect_equal(resolve_point_indices(c(2, 2, 3), m, "x"), c(2L, 3L))

  expect_error(resolve_point_indices("nope", m, "x"),
               "x must name points that exist")
  expect_error(resolve_point_indices(c(0, 1), m, "x"),
               "x must be row indices between 1 and")
  expect_error(resolve_point_indices(character(0), m, "x"),
               "x must not be empty")
})

test_that(".onUnload releases the compiled shared library", {
  hook <- topolow:::.onUnload
  expect_true(is.function(hook))
  expect_identical(names(formals(hook)), "libpath")
  expect_true(any(grepl("library.dynam.unload", deparse(body(hook)),
                        fixed = TRUE)))
})
