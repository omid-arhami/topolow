# tests/testthat/test-integration.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("full data processing and optimization workflow executes on general data", {
  # 1. Create temporary data file
  test_data <- data.frame(
    object = rep(paste0("V", 1:4), each = 4),
    reference = rep(paste0("S", 1:4), 4),
    score = sample(c(10, 20, 40, 80, 160, 320, 640, "<10", ">1280"), 16, replace = TRUE)
  )
  # 2.  process it
  processed_matrix <- list_to_matrix(
    data = test_data,  # Pass the data frame, not file path
    object_col = "object",
    reference_col = "reference",
    value_col = "score",
    is_similarity = TRUE
  )
  expect_true(is.matrix(processed_matrix))

  # 3. Run optimization
  result <- euclidean_embedding(
    dissimilarity_matrix = processed_matrix,
    ndim = 2,
    mapping_max_iter = 50, # Keep low for testing
    k0 = 5.0,
    cooling_rate = 0.05,
    c_repulsion = 0.01,
    write_positions_to_csv = FALSE
  )
  expect_s3_class(result, "topolow")

  # 4. Run a diagnostic plot
  plots <- scatterplot_fitted_vs_true(
    dissimilarity_matrix = processed_matrix,
    p_dissimilarity_mat = result$est_distances,
    save_plot = FALSE
  )
  expect_s3_class(plots$scatter_plot, "ggplot")
})

test_that("full workflow executes correctly on antigenic data", {
  # Create a temporary file path
  temp_csv_path <- tempfile(fileext = ".csv")
  # Create test data
  test_data <- data.frame(
    antigen = rep(paste0("V", 1:3), each=3),
    serum = rep(paste0("S", 1:3), 3),
    titer = c(40, 80, 160, 320, "<10", ">640", 40, 80, 160),
    virusYear = rep(2000:2002, each=3),
    serumYear = rep(2000:2002, 3),
    cluster = rep(c("A", "B", "C"), 3),
    color = rep(c("red", "blue", "green"), 3)
  )
  
  # Add year to antigen/serum names before writing
  test_data$antigen <- paste0(test_data$antigen, "/", test_data$virusYear)
  test_data$serum <- paste0(test_data$serum, "/", test_data$serumYear)
  
  # Process data using the temporary file path
  results <- process_antigenic_data(
    test_data,
    antigen_col = "antigen",
    serum_col = "serum",
    value_col = "titer",
    is_similarity = TRUE,
    metadata_cols = c("cluster", "color", "virusYear", "serumYear")
  )
  
  # Run optimization
  topo_result <- euclidean_embedding(dissimilarity_matrix =  results$matrix,
    ndim = 2,
    mapping_max_iter = 100,
    k0 = 3.0,
    cooling_rate = 0.1,
    c_repulsion = 0.001,
    write_positions_to_csv = FALSE
  )
  
  # Create visualization
  positions <- as.data.frame(topo_result$positions)
  positions$name <- rownames(positions)
  positions$antigen <- grepl("^V/", rownames(positions))
  positions$antiserum <- grepl("^S/", rownames(positions))
  
  # Extract year based on point type
  positions$year <- sapply(rownames(positions), function(x) {
    if (grepl("^V/", x)) {
      # Extract year from virus name
      year <- as.numeric(sub(".*/(\\d{4})$", "\\1", x))
    } else {
      # Extract year from serum name
      year <- as.numeric(sub(".*/(\\d{4})$", "\\1", x))
    }
    return(year)
  })
  
  # Verify year extraction worked
  expect_true(!any(is.na(positions$year)))
  expect_true(all(positions$year %in% 2000:2002))
  
  plot <- plot_temporal_mapping(positions, ndim = 2)
  expect_s3_class(plot, "ggplot")

})

test_that("parameter optimization workflow works", {
  # Create test matrix with enough variation
  test_mat <- matrix(c(0, 1, 2, 3, 1, 0, 2.5, 3.5, 2, 2.5, 0, 4, 3, 3.5, 4, 0), 4, 4)
  rownames(test_mat) <- colnames(test_mat) <- paste0("Point", 1:4)

  # Run parameter optimization with minimal settings for testing
  results <- initial_parameter_optimization(
    dissimilarity_matrix = test_mat,
    mapping_max_iter = 50,
    relative_epsilon = 1e-3,
    convergence_counter = 3,
    scenario_name = "test_opt",
    N_min = 2,
    N_max = 3,
    k0_min = 0.5,
    k0_max = 5,
    c_repulsion_min = 0.001,
    c_repulsion_max = 0.01,
    cooling_rate_min = 0.001,
    cooling_rate_max = 0.05,
    num_samples = 2, # Reduced for testing
    folds = 2,       # Reduced for testing
    max_cores = 1,
    write_files = FALSE
  )

  expect_true(is.data.frame(results))
  expect_true(all(c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion", "Holdout_MAE", "NLL") %in%
                    names(results)))
})

test_that("adaptive sampling workflow executes", {
  # Create a temporary file path for the samples
  temp_samples_path <- tempfile(fileext = ".csv")

  # Create initial samples with log-transformed names
  samples <- data.frame(
    log_N = log(c(3, 4)),
    log_k0 = log(c(1, 1.5)),
    log_cooling_rate = log(c(0.001, 0.002)),
    log_c_repulsion = log(c(0.001, 0.002)),
    NLL = c(100, 90),
    Holdout_MAE = c(2, 1.8)
  )
  write.csv(samples, temp_samples_path, row.names = FALSE)

  # Create test dissimilarity matrix
  test_mat <- as.matrix(dist(matrix(rnorm(10 * 3), ncol = 3)))

  # Run adaptive sampling using the temporary file
  # We test the internal function as it's the core of the logic
  result <- adaptive_MC_sampling(
    samples_file = temp_samples_path,
    dissimilarity_matrix = test_mat,
    iterations = 1, # Just one iteration for a quick test
    mapping_max_iter = 10,
    relative_epsilon = 1e-3,
    folds = 2,
    scenario_name = "test_amc",
    verbose = FALSE
  )

  expect_true(is.data.frame(result))
  expect_true(nrow(result) > nrow(samples)) # Check that samples were added
  expect_true(all(names(samples) %in% names(result)))

  unlink(temp_samples_path)
})



test_that("adaptive sampling workflow executes", {
  # Create a temporary file path for the samples
  temp_samples_path <- tempfile(fileext = ".csv")

  # Create initial samples with log-transformed names
  samples <- data.frame(
    log_N = log(c(3, 4)),
    log_k0 = log(c(1, 1.5)),
    log_cooling_rate = log(c(0.001, 0.002)),
    log_c_repulsion = log(c(0.001, 0.002)),
    NLL = c(100, 90),
    Holdout_MAE = c(2, 1.8)
  )
  write.csv(samples, temp_samples_path, row.names = FALSE)

  # Create test dissimilarity matrix
  test_mat <- as.matrix(dist(matrix(rnorm(10 * 3), ncol = 3)))

  # Run adaptive sampling using the temporary file
  # We test the internal function as it's the core of the logic
  result <- adaptive_MC_sampling(
    samples_file = temp_samples_path,
    dissimilarity_matrix = test_mat,
    iterations = 1, # Just one iteration for a quick test
    mapping_max_iter = 10,
    relative_epsilon = 1e-3,
    folds = 2,
    scenario_name = "test_amc",
    verbose = FALSE
  )

  expect_true(is.data.frame(result))
  expect_true(nrow(result) > nrow(samples)) # Check that samples were added
  expect_true(all(names(samples) %in% names(result)))

  unlink(temp_samples_path)
})
