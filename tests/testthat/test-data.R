test_that("example_positions dataset loads correctly", {
  # Check data loads
  data(example_positions)
  
  # Check dimensions 
  expect_equal(ncol(example_positions), 11)
  expect_true(nrow(example_positions) > 0)
  
  # Check required columns exist
  expect_true(all(paste0("V", 1:5) %in% names(example_positions)))
  expect_true(all(c("name", "antigen", "antiserum", "cluster", 
                    "color", "year") %in% names(example_positions)))
  
  # Check data types
  expect_true(is.logical(example_positions$antigen))
  expect_true(is.logical(example_positions$antiserum))
  expect_true(is.factor(example_positions$cluster))
  
  # Check no missing values
  expect_true(!any(is.na(example_positions)))
})

test_that("h3n2_data loads correctly", {
  data(h3n2_data)
  expect_true(all(c("virusStrain", "serumStrain", "titer", 
                    "virusYear", "serumYear", "cluster", "color") %in% 
                    names(h3n2_data)))
  expect_true(is.factor(h3n2_data$cluster))
})

test_that("HIV data loads correctly", {
  data(hiv_viruses)
  data(hiv_titers)
  
  expect_true(all(c("Virus.name", "Country", "Subtype", "Year") %in% 
                    names(hiv_viruses)))
  expect_true(all(c("Antibody", "Virus", "IC50") %in% 
                    names(hiv_titers)))
})