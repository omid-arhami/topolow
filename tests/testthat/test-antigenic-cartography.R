test_that("process_antigenic_data handles various input formats", {
  # Create test data
  test_data <- data.frame(
    antigen = c("A1", "A2", "A1"),
    serum = c("S1", "S1", "S2"),
    value = c(40, 80, "<160"),
    stringsAsFactors = FALSE
  )

  # Test titer data processing
  result <- process_antigenic_data(
    test_data,
    antigen_col = "antigen",
    serum_col = "serum",
    value_col = "value",
    is_similarity = TRUE
  )
  
  expect_true(is.list(result))
  expect_true(all(c("long", "matrix") %in% names(result)))
  expect_true(is.matrix(result$matrix))
})
