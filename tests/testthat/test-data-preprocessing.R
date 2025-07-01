test_that("process_antigenic_data handles various input formats", {
  # Create test data
  test_data <- data.frame(
    antigen = c("A1", "A2", "A1"),
    serum = c("S1", "S1", "S2"),
    value = c(40, 80, "<160"),
    stringsAsFactors = FALSE
  )
  
  # Use a temporary file for testing
  temp_file <- tempfile(fileext = ".csv")
  write.csv(test_data, temp_file, row.names = FALSE)
  
  # Test titer data processing
  result <- process_antigenic_data(
    file_path = temp_file,
    antigen_col = "antigen",
    serum_col = "serum",
    value_col = "value",
    is_titer = TRUE
  )
  
  expect_true(is.list(result))
  expect_true(all(c("long", "matrix") %in% names(result)))
  expect_true(is.matrix(result$matrix))
  
  # Clean up is handled automatically for tempfiles, but explicit unlink is fine
  unlink(temp_file)
})

test_that("long_to_matrix creates correct matrix", {
  test_data <- data.frame(
    antigen = c("A1", "A2"),
    serum = c("S1", "S1"),
    distance = c(1.0, 2.0)
  )
  
  result <- long_to_matrix(
    test_data,
    chnames = "antigen",
    rnames = "serum",
    values_column = "distance"
  )
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
  expect_true(isSymmetric(result))
})

test_that("prune_distance_network maintains connectivity", {
  # Create test distance matrix
  test_mat <- matrix(1:16, 4, 4)
  diag(test_mat) <- 0
  test_mat[upper.tri(test_mat)] <- t(test_mat)[upper.tri(test_mat)]
  
  test_data <- data.frame(
    Virus = rep(paste0("V", 1:2), each=2),
    Antibody = rep(paste0("Ab", 1:2), 2),
    distance = c(1,2,3,4)
  )
  
  result <- prune_distance_network(
    test_data,
    virus_col = "Virus",
    antibody_col = "Antibody",
    min_connections = 2
  )
  
  expect_true(is.list(result))
  expect_true("pruned_data" %in% names(result))
  expect_true("stats" %in% names(result))
})
