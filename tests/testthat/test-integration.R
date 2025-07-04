test_that("full workflow executes correctly", {
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
  
  # Write to the temporary path
  write.csv(test_data, temp_csv_path, row.names = FALSE)
  
  # Process data using the temporary file path
  results <- process_antigenic_data(
    temp_csv_path,
    antigen_col = "antigen",
    serum_col = "serum",
    value_col = "titer",
    is_titer = TRUE,
    metadata_cols = c("cluster", "color", "virusYear", "serumYear")
  )
  
  # Run optimization
  topo_result <- create_topolow_map(
    distance_matrix = results$matrix,
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
  test_mat <- matrix(c(0,1,2,3, 1,0,2.5,3.5, 2,2.5,0,4, 3,3.5,4,0), 4, 4)
  rownames(test_mat) <- colnames(test_mat) <- paste0("Point", 1:4)
  
  # Run parameter optimization with minimal settings for testing
  results <- initial_parameter_optimization(
    distance_matrix = test_mat,
    mapping_max_iter = 100,
    relative_epsilon = 1e-3,
    convergence_counter = 3,
    scenario_name = "test_opt",
    N_min = 2,
    N_max = 4,
    k0_min = 0.5,
    k0_max = 5,
    c_repulsion_min = 0.001,
    c_repulsion_max = 0.01,
    cooling_rate_min = 0.001,
    cooling_rate_max = 0.05,
    num_samples = 3, # Reduced for testing
    folds = 3,      # Reduced for testing
    max_cores = 1,
    write_files = FALSE
  )
  
  expect_true(is.data.frame(results))
  expect_true(all(c("N", "k0", "cooling_rate", "c_repulsion", "Holdout_MAE", "NLL") %in% 
                    names(results)))
})


test_that("adaptive sampling workflow executes", {
  # Create a temporary file path for the samples
  temp_samples_path <- tempfile(fileext = ".csv")
  # Create initial samples
  samples <- data.frame(
    log_N = log(c(5,3,4,3, 3.5, 4.5, 2.5, 5.5, 3.2)),
    log_k0 = log(c(1,1.5,2,0.5, 1.2, 1.8, 0.8, 2.2, 1.1)),
    log_cooling_rate = log(c(0.001, 0.002, 0.003, 0.0005, 0.0015, 0.0025, 0.0008, 0.0032, 0.0011)),
    log_c_repulsion = log(c(0.001, 0.002, 0.003, 0.0001, 0.0012, 0.0022, 0.0002, 0.0035, 0.0015)),
    NLL = c(100, 90, 95, 120, 92, 98, 110, 85, 105),
    Holdout_MAE = c(2, 1.8, 1.9, 2.8, 1.85, 1.95, 2.6, 1.75, 2.1)
  )
  
  # Write the samples to the temporary file
  write.csv(samples, temp_samples_path, row.names = FALSE)
  
  # Create test distance matrix
  test_mat <- matrix(0, 10, 10)
  
  for (i in 1:10) {
    for (j in 1:10) {
      if (i != j) {
        test_mat[i, j] <- abs(i - j) + sample(0:3, 1) # Introducing slight variation
      }
    }
  }
  
  # Run adaptive sampling using the temporary file
  result <- adaptive_MC_sampling(
    samples_file = temp_samples_path,
    distance_matrix = test_mat,
    iterations = 1,
    mapping_max_iter = 10,
    relative_epsilon = 1e-3,
    folds = 3,
    scenario_name = "test_amc"
  )
  
  expect_true(is.data.frame(result))
  expect_true(all(c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion", 
                    "NLL", "Holdout_MAE") %in% names(result)))

})
