test_that("full workflow executes correctly", {
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
  
  write.csv(test_data, "test_workflow.csv", row.names = FALSE)
  
  # Process data
  results <- process_antigenic_data(
    "test_workflow.csv",
    antigen_col = "antigen",
    serum_col = "serum",
    value_col = "titer",
    is_titer = TRUE,
    metadata_cols = c("cluster", "color", "virusYear", "serumYear")
  )
  
  # Run optimization
  topo_result <- topolow_full(
    distance_matrix = results$matrix,
    ndim = 2,
    max_iter = 100,
    k0 = 3.0,
    k_decay = 0.1,
    cqq = 0.001
  )
  
  # Create visualization
  positions <- as.data.frame(topo_result$positions)
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
  
  # Clean up
  unlink("test_workflow.csv")
})


test_that("parameter optimization workflow works", {
  # Create test matrix with enough variation
  test_mat <- matrix(c(0,1,2,3, 1,0,2.5,3.5, 2,2.5,0,4, 3,3.5,4,0), 4, 4)
  rownames(test_mat) <- colnames(test_mat) <- paste0("Point", 1:4)
  
  # Run parameter optimization with minimal settings for testing
  results <- run_parameter_optimization(
    distance_matrix = test_mat,
    max_iter = 100,
    relative_epsilon = 1e-3,
    convergence_counter = 3,
    scenario_name = "test_opt",
    N_min = 2,
    N_max = 4,
    k0_min = 0.5,
    k0_max = 5,
    cqq_min = 0.001,
    cqq_max = 0.01,
    k_decay_min = 0.001,
    k_decay_max = 0.05,
    num_samples = 5, # Reduced for testing
    folds = 3,      # Reduced for testing
    write_files = FALSE,
    num_cores = 1   # Force single core for testing
  )
  
  expect_true(is.data.frame(results))
  expect_true(all(c("N", "k0", "k_decay", "cqq", "Holdout_MAE", "NLL") %in% 
                    names(results)))
})


test_that("adaptive sampling workflow executes", {
  # Create initial samples
  samples <- data.frame(
    log_N = log(c(2,3,4)),
    log_k0 = log(c(1,1.5,2)),
    log_k_decay = log(c(0.001, 0.002, 0.003)),
    log_cqq = log(c(0.001, 0.002, 0.003)),
    NLL = c(100, 90, 95),
    Holdout_MAE = c(2, 1.8, 1.9)
  )
  write.csv(samples, "test_samples.csv", row.names = FALSE)
  
  # Create test distance matrix
  test_mat <- matrix(c(0,1,2, 1,0,3, 2,3,0), 3, 3)
  
  # Run adaptive sampling
  result <- adaptive_MC_sampling(
    samples_file = "test_samples.csv",
    distance_matrix = test_mat,
    n_iter = 1,
    batch_size = 1,
    max_iter = 10,
    relative_epsilon = 1e-3,
    folds = 3,
    num_cores = 1,
    scenario_name = "test_amc"
  )
  
  expect_true(is.data.frame(result))
  expect_true(all(c("log_N", "log_k0", "log_k_decay", "log_cqq", 
                    "NLL", "Holdout_MAE") %in% names(result)))
  
  # Clean up
  unlink("test_samples.csv")
})