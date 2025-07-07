# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/adaptive_sampling.R

#' Latin Hypercube and Adaptive Monte Carlo Sampling Functions
#' 
#' @description
#' This file contains functions for performing Latin Hypercube and adaptive Monte Carlo 
#' sampling in parameter space. The AMC sampling adapts based on previous evaluations to focus
#' sampling in high-likelihood regions. The functions run locally using parallel processing.
#' 
#' Functions handle:
#' - A suite of functions to get an initial estimate of the likelihood space through LHS 
#' - Core adaptive sampling algorithm 
#' - Likelihood calculations with cross-validation
#' - Distribution updating and resampling
#' - Safe evaluation wrappers
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats sd
#' @importFrom parallel mclapply detectCores
#' @keywords internal
"_PACKAGE"


#' Run Parameter Optimization Via Latin Hypercube Sampling
#'
#' @description
#' Performs parameter optimization using Latin Hypercube Sampling (LHS) combined with
#' k-fold cross-validation. Parameters are sampled from specified ranges using maximin
#' LHS design to ensure good coverage of parameter space. Each parameter set is evaluated
#' using k-fold cross-validation to assess prediction accuracy. To calculate one NLL per set of
#' parameters, the function uses a pooled errors approach which combine all validation errors into 
#' one set, then calculate a single NLL. This approach has two main advantages:
#' 1- It treats all validation errors equally, respecting the underlying error distribution assumption
#' 2- It properly accounts for the total number of validation points
#'
#' @details
#' The function performs these steps:
#' 1. Generates LHS samples in parameter space 
#' 2. Creates k-fold splits of input data
#' 3. For each parameter set and fold:
#'    - Trains model on training set
#'    - Evaluates on validation set
#'    - Calculates MAE and negative log likelihood
#' 4. Computations are run locally in parallel.
#'
#' Parameters ranges are transformed to log scale where appropriate to handle
#' different scales effectively.
#'
#' @param distance_matrix Matrix or data frame. Input distance matrix. Must be square
#'        and symmetric. Can contain NA values for missing measurements.
#' @param mapping_max_iter Integer. Maximum number of optimization iterations.
#' @param relative_epsilon Numeric. Convergence threshold for relative change in error.
#' @param convergence_counter Integer. Number of iterations below threshold before 
#'        declaring convergence.
#' @param scenario_name Character. Name for output files and job identification.
#' @param N_min,N_max Integer. Range for number of dimensions parameter.
#' @param k0_min,k0_max Numeric. Range for initial spring constant parameter.
#' @param c_repulsion_min,c_repulsion_max Numeric. Range for repulsion constant parameter.
#' @param cooling_rate_min,cooling_rate_max Numeric. Range for spring decay parameter.
#' @param num_samples Integer. Number of LHS samples to generate (default: 20).
#' @param max_cores Integer. Maximum number of cores to use for parallel processing. If NULL,
#'        uses all available cores minus 1 (default: NULL).
#' @param folds Integer. Number of cross-validation folds. Default: 20.
#' @param verbose Logical. Whether to print progress messages. Default: FALSE.
#' @param write_files Logical. Whether to save results to CSV. Default: FALSE.
#' @param output_dir Character. Directory where output files will be saved. 
#'        Required if `write_files` is `TRUE`.
#' 
#' @return A `data.frame` containing the parameter sets and their performance metrics 
#' (`Holdout_MAE` and `NLL`). The columns of the data frame are `N`, `k0`, 
#' `cooling_rate`, `c_repulsion`, `Holdout_MAE`, and `NLL`.
#' If `write_files` is `TRUE`, this data frame is also saved to a CSV file as a side effect.
#'
#' @examples
#' \donttest{
#' # This example is wrapped in \donttest{} because it can exceed 5 seconds,
#' # 1. Create a structured, synthetic dataset for the example
#' # Generate coordinates for a more realistic test case
#' synth_coords <- generate_complex_data(n_points = 20, n_dim = 3)
#' # Convert coordinates to a distance matrix
#' dist_mat <- coordinates_to_matrix(synth_coords)
#'
#' # 2. Run the optimization on the synthetic data
#' # ensuring it passes CRAN's automated checks.
#' results <- initial_parameter_optimization(
#'   distance_matrix = dist_mat,
#'   mapping_max_iter = 100,
#'   relative_epsilon = 1e-3,
#'   convergence_counter = 2,
#'   scenario_name = "test_opt_synthetic",
#'   N_min = 2, N_max = 5,
#'   k0_min = 1, k0_max = 10,
#'   c_repulsion_min = 0.001, c_repulsion_max = 0.05,
#'   cooling_rate_min = 0.001, cooling_rate_max = 0.02,
#'   num_samples = 4,
#'   max_cores = 2,
#'   verbose = FALSE
#' )
#' print(results)
#' }
#' 
#' @seealso
#' \code{\link{create_topolow_map}} for the core optimization algorithm
#'
#' @importFrom lhs maximinLHS
#' @importFrom stats qunif complete.cases aggregate
#' @importFrom parallel detectCores mclapply makeCluster clusterExport clusterEvalQ parLapply stopCluster
#' @export
initial_parameter_optimization <- function(# Mapping related arguments:
                                          distance_matrix,
                                          mapping_max_iter = 1000,
                                          relative_epsilon,
                                          convergence_counter,
                                          scenario_name,
                                          N_min,
                                          N_max,
                                          k0_min,
                                          k0_max,
                                          c_repulsion_min,
                                          c_repulsion_max,
                                          cooling_rate_min,
                                          cooling_rate_max,
                                          # Sampling related arguments:
                                          num_samples = 20,
                                          max_cores = NULL,
                                          folds = 20,
                                          verbose = FALSE,
                                          write_files = FALSE,
                                          output_dir) {
  # Input validation
  if (!is.matrix(distance_matrix)) {
    stop("distance_matrix must be a matrix")
  }
  if (nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix must be square")
  }
  
  # Validate integer parameters
  integer_params <- list(
    mapping_max_iter = mapping_max_iter,
    convergence_counter = convergence_counter,
    N_min = N_min,
    N_max = N_max,
    num_samples = num_samples,
    folds = folds
  )
  
  for (param_name in names(integer_params)) {
    param_value <- integer_params[[param_name]]
    if (!is.numeric(param_value) || 
        param_value != round(param_value) || 
        param_value < 1) {
      stop(sprintf("%s must be a positive integer", param_name))
    }
  }
  
  # Validate numeric parameters
  numeric_params <- list(
    relative_epsilon = relative_epsilon,
    k0_min = k0_min,
    k0_max = k0_max,
    c_repulsion_min = c_repulsion_min,
    c_repulsion_max = c_repulsion_max,
    cooling_rate_min = cooling_rate_min,
    cooling_rate_max = cooling_rate_max
  )
  
  for (param_name in names(numeric_params)) {
    param_value <- numeric_params[[param_name]]
    if (!is.numeric(param_value) || param_value <= 0) {
      stop(sprintf("%s must be a positive number", param_name))
    }
  }
  
  # Validate ranges
  if (N_min >= N_max) {
    stop("N_min must be less than N_max")
  }
  if (k0_min >= k0_max) {
    stop("k0_min must be less than k0_max")
  }
  if (c_repulsion_min >= c_repulsion_max) {
    stop("c_repulsion_min must be less than c_repulsion_max")
  }
  if (cooling_rate_min >= cooling_rate_max) {
    stop("cooling_rate_min must be less than cooling_rate_max")
  }
  
  # Validate logical parameters
  logical_params <- list(
    verbose = verbose,
    write_files = write_files
  )
  
  for (param_name in names(logical_params)) {
    param_value <- logical_params[[param_name]]
    if (!is.logical(param_value) || length(param_value) != 1) {
      stop(sprintf("%s must be a single logical value", param_name))
    }
  }
  
  # Validate scenario name
  if (!is.character(scenario_name) || length(scenario_name) != 1) {
    stop("scenario_name must be a single character string")
  }
  
  # Check for output_dir if writing files
  if (write_files && missing(output_dir)) {
    stop("An 'output_dir' must be provided when 'write_files' is TRUE.", call. = FALSE)
  }
  
  # Determine maximum cores for parallel processing
  available_cores <- parallel::detectCores()
  if (is.null(max_cores)) {
    # Use all cores minus 1 to avoid bogging down the system
    max_cores <- max(1, available_cores - 1)
  } else {
    # Validate max_cores
    if (!is.numeric(max_cores) || max_cores < 1 || max_cores != round(max_cores)) {
      stop("max_cores must be a positive integer")
    }
    # Limit to available cores
    max_cores <- min(max_cores, available_cores)
  }
  
  if (verbose) {
    cat(sprintf("Processing %d samples with maximum %d cores\n", num_samples, max_cores))
  }
  
  # Only set up directories if writing files
  param_dir <- NULL
  run_topolow_dir <- NULL
  if (write_files) {
    param_dir <- file.path(output_dir, "model_parameters")
    run_topolow_dir <- file.path(output_dir, "init_param_optimization")
    for (dir_path in c(param_dir, run_topolow_dir)) {
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      }
    }
    # Verify write permissions
    test_file <- file.path(param_dir, "test_write.txt")
    tryCatch({
      write.table(data.frame(test=1), test_file)
      unlink(test_file)
    }, error = function(e) {
      stop("No write permission in directory: ", param_dir)
    })
  }

  # Generate LHS samples
  lhs_samples <- maximinLHS(n = num_samples, k = 4)
  lhs_params <- data.frame(
    N = floor(qunif(lhs_samples[,1], min = N_min, max = N_max + 1)),
    k0 = qunif(lhs_samples[,2], min = k0_min, max = k0_max),
    c_repulsion = qunif(lhs_samples[,3], min = c_repulsion_min, max = c_repulsion_max),
    cooling_rate = qunif(lhs_samples[,4], min = cooling_rate_min, max = cooling_rate_max)
  )
  
  # Create cross-validation folds
  # Initialize list for storing fold data
  matrix_list <- vector("list", folds)
  for (i in 1:folds) {
    matrix_list[[i]] <- list(
      truth_matrix = distance_matrix,  # Original matrix
      train_matrix = NULL  # Will hold training data with validation set masked
    )
  }
  
  # Calculate holdout size
  num_elements <- sum(!is.na(distance_matrix))
  holdout_size <- floor(num_elements/(folds*2)) # Factor of 2 accounts for symmetry
  
  # Create training matrices for each fold
  D_train <- distance_matrix  # Copy for progressive masking
  
  for(i in 1:folds) {
    if(verbose) {
      cat(sprintf("Creating fold %d/%d\n", i, folds))
    }
    
    # Randomly select elements for validation
    random_indices <- sample(which(!is.na(D_train)), size=holdout_size)
    
    # Create training matrix for this fold
    input_matrix <- distance_matrix
    for(index in random_indices) {
      # Convert linear index to row/column
      row <- (index - 1) %/% nrow(distance_matrix) + 1
      col <- (index - 1) %% ncol(distance_matrix) + 1
      # Mask validation elements
      input_matrix[row, col] <- NA
      input_matrix[col, row] <- NA  # Maintain symmetry
    }
    
    # Store training matrix
    matrix_list[[i]][[2]] <- input_matrix
    
    # Update D_train by masking used elements
    for(index in random_indices) {
      row <- (index - 1) %/% nrow(D_train) + 1
      col <- (index - 1) %% ncol(D_train) + 1
      D_train[row, col] <- NA
      D_train[col, row] <- NA
    }
  }

  # Process samples locally
  if(verbose) cat("Processing samples locally\n")
  
  # Determine the optimal number of cores to use
  local_cores <- min(num_samples, max_cores)
  if(verbose) cat(sprintf("Using %d cores for processing %d samples\n", local_cores, num_samples))
  
  # Process each sample and fold
  process_sample <- function(i) {
      sample_idx <- ((i - 1) %% num_samples) + 1
      fold_idx <- floor((i - 1) / num_samples) + 1
      
      N <- lhs_params$N[sample_idx]
      k0 <- lhs_params$k0[sample_idx]
      c_repulsion <- lhs_params$c_repulsion[sample_idx]
      cooling_rate <- lhs_params$cooling_rate[sample_idx]
      
      truth_matrix <- matrix_list[[fold_idx]][[1]]
      input_matrix <- matrix_list[[fold_idx]][[2]]
      
      tryCatch({
          res_train <- create_topolow_map(
          distance_matrix = input_matrix,
          ndim = N,
          mapping_max_iter = mapping_max_iter,
          k0 = k0,
          cooling_rate = cooling_rate,
          c_repulsion = c_repulsion,
          relative_epsilon = relative_epsilon,
          convergence_counter = convergence_counter,
          initial_positions = NULL,
          write_positions_to_csv = FALSE,
          verbose = FALSE
          )
          
          p_dist_mat <- as.matrix(res_train$est_distances)
          
          errors <- error_calculator_comparison(
          p_dist_mat = p_dist_mat,
          truth_matrix = truth_matrix,
          input_matrix = input_matrix
          )
          
          df <- errors$report_df
          
          # Store data needed for per-fold and pooled calculations
          out_sample_errors <- df$OutSampleError[!is.na(df$OutSampleError)]
          n_samples <- length(out_sample_errors)
          sum_abs_errors <- sum(abs(out_sample_errors))
          
          # Calculate fold-specific MAE
          mae_holdout <- if(n_samples > 0) sum_abs_errors / n_samples else NA
          
          # Return valid results with temporary additional columns for pooling
          if(is.finite(mae_holdout) && n_samples > 0) {
          # Include the pooling data as temporary columns
          result <- data.frame(
              N = N,
              k0 = k0,
              cooling_rate = cooling_rate,
              c_repulsion = c_repulsion,
              Holdout_MAE = mae_holdout,
              NLL = n_samples * (1 + log(2*mae_holdout)),
              # Temporary columns for pooling calculation
              temp_n_samples = n_samples,
              temp_sum_abs_errors = sum_abs_errors
          )
          
          # Save individual result if requested - save only the standard columns
          if(write_files) {
              result_file <- file.path(run_topolow_dir,
                                      sprintf("%d_params_%s.csv", i, scenario_name))
              write.csv(result[1:6], result_file, row.names = FALSE)
          }
          
          return(result)
          } else {
          if(verbose) {
              cat(sprintf("Sample %d produced invalid results (inf/NA)\n", i))
          }
          return(NULL)
          }
          
      }, error = function(e) {
          if(verbose) {
          cat(sprintf("Error processing sample %d: %s\n", i, e$message))
          }
          return(NULL)
      })
  }
  
  # Create batches if num_samples*folds exceeds what we can process at once
  total_combinations <- num_samples * folds
  if(verbose) cat(sprintf("Total cross validations to process: %d\n", total_combinations))
  
  # Determine batch size based on available cores and memory considerations
  # Process in reasonable sized batches to avoid memory issues
  batch_size <- min(local_cores * 10, total_combinations)
  num_batches <- ceiling(total_combinations / batch_size)
  
  all_results <- list()
  
  for(batch in 1:num_batches) {
    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, total_combinations)
    batch_indices <- batch_start:batch_end
    
    if(verbose) {
      cat(sprintf("Processing batch %d/%d (indices %d-%d)\n", 
                  batch, num_batches, batch_start, batch_end))
    }
    
    # Process all samples in current batch with appropriate parallel method
    if(local_cores > 1) {
      if(.Platform$OS.type == "windows") {
        if(verbose) cat("Using parallel cluster on Windows\n")
        # Create cluster
        cl <- parallel::makeCluster(local_cores)
        
        # Export required functions and data to cluster
        parallel::clusterExport(cl, c("matrix_list", "lhs_params", "mapping_max_iter",
                                      "relative_epsilon", "convergence_counter",
                                      "scenario_name", "write_files", "verbose",
                                      "run_topolow_dir", "param_dir", "num_samples"),
                                envir = environment())
        
        # Load required packages on each cluster node
        parallel::clusterEvalQ(cl, {
          library(topolow)
        })
        
        # Run processing
        batch_results <- parallel::parLapply(cl, batch_indices, process_sample)
        
        # Clean up
        parallel::stopCluster(cl)
      } else {
        # Use mclapply on Unix-like systems
        if(verbose) cat("Using mclapply on Unix-like system\n")
        batch_results <- parallel::mclapply(
          batch_indices,
          process_sample,
          mc.cores = local_cores
        )
      }
    } else {
      # Sequential processing
      if(verbose) cat("Running sequentially with single core\n")
      batch_results <- lapply(batch_indices, process_sample)
    }
    
    # Add batch results to all_results
    all_results <- c(all_results, batch_results)
    
    # Clean memory between batches
    gc()
  }
  
  # Remove NULL results
  res_list <- Filter(Negate(is.null), all_results)
  # Check if we have any valid results
  if (length(res_list) == 0) stop("No valid results obtained")
  # Combine results
  res_list <- do.call(rbind, res_list)
  # Remove any remaining invalid values
  res_list <- res_list[complete.cases(res_list) & 
                      apply(res_list, 1, function(x) all(is.finite(x))), ]
  if(nrow(res_list) == 0) {
  stop("All results were invalid after filtering infinities and NAs")
  }

  # Calculate pooled statistics using aggregate
  pooled_results <- aggregate(
          cbind(temp_sum_abs_errors, temp_n_samples) ~ N + k0 + cooling_rate + c_repulsion,
          data = res_list,
          FUN = sum
  )

  # Calculate pooled MAE and NLL
  pooled_results$Holdout_MAE <- pooled_results$temp_sum_abs_errors / pooled_results$temp_n_samples
  pooled_results$NLL <- pooled_results$temp_n_samples * (1 + log(2*pooled_results$Holdout_MAE))

  # Remove temporary columns for final output
  res_list_median <- pooled_results[, c("N", "k0", "cooling_rate", "c_repulsion", "Holdout_MAE", "NLL")]

  # --- RENAME & REORDER FOR OUTPUT ---
  final_df <- res_list_median
  colnames(final_df) <- c(
    "N",
    "k0",
    "cooling_rate",
    "c_repulsion",
    "Holdout_MAE",
    "NLL"
  )

  # Write aggregated results
  if (write_files) {
      file_name <- file.path(param_dir, paste0(scenario_name, "_model_parameters.csv"))
      if (file.exists(file_name)) {
        existing_data <- read.csv(file_name, stringsAsFactors = FALSE)
        combined      <- rbind(existing_data, final_df)
        write.csv(combined, file_name, row.names = FALSE)
      } else {
        write.csv(final_df, file_name, row.names = FALSE)
      }
  }
  
  return(final_df)
}




#' Run Adaptive Monte Carlo Sampling
#'
#' @description 
#' Performs adaptive Monte Carlo sampling to explore parameter space, running locally
#' in parallel. Samples are drawn adaptively based on previous evaluations
#' to focus sampling in high-likelihood regions. Results from all parallel jobs 
#' accumulate in a single output file. This function always writes to the file system
#' and therefore requires the `output_dir` argument.
#'
#' @param initial_samples_file Character. Path to CSV file containing initial samples.
#'        Must contain columns: log_N, log_k0, log_cooling_rate, log_c_repulsion, NLL
#' @param distance_matrix Matrix. Distance matrix of the input data.
#' @param mapping_max_iter Integer. Maximum iterations per map optimization.
#' @param relative_epsilon Numeric. Convergence threshold.
#' @param folds Integer. Number of CV folds (default: 20).
#' @param num_parallel_jobs Integer. Number of parallel local jobs (chains) to run.
#' @param max_cores Integer. Maximum number of cores to use for parallel processing across all jobs. 
#'        If NULL, uses all available cores minus 1 (default: NULL).
#' @param num_samples Integer. Number of new samples to be added to the CSV file containing initial samples 
#'        through Adaptive Monte Carlo sampling (default: 10).
#' @param scenario_name Character. Name for output files.
#' @param output_dir Character. Directory for output job files. The project's working directory is a 
#'        straightforward example. This argument is required.
#' @param verbose Logical. Whether to print progress messages. Default: FALSE.
#'
#' @return No return value, called for side effects. The function writes the 
#' results of the adaptive sampling to a CSV file within the specified `output_dir`.
#'
#' @examples
#' \donttest{
#' # 1. Locate the example initial samples file included with the package
#' initial_file <- system.file(
#'   "extdata", "initial_samples_example.csv", 
#'   package = "topolow"
#' )
#' 
#' # 2. Create a temporary directory for the function's output
#' # This function requires a writable directory for its results.
#' temp_out_dir <- tempdir()
#' 
#' # 3. Create a sample distance matrix for the function to use
#' dist_mat <- matrix(runif(100, 1, 10), 10, 10)
#' diag(dist_mat) <- 0
#' 
#' # 4. Run the adaptive sampling only if the example file is found
#' if (nzchar(initial_file)) {
#'   run_adaptive_sampling(
#'     initial_samples_file = initial_file,
#'     scenario_name = "adaptive_test_example",
#'     distance_matrix = dist_mat,
#'     output_dir = temp_out_dir,
#'     num_parallel_jobs = 2, # Use small values for a quick example
#'     num_samples = 2,
#'     verbose = FALSE
#'   )
#'   
#'   # 5. Verify output files were created
#'   print("Output files from adaptive sampling:")
#'   print(list.files(temp_out_dir, recursive = TRUE))
#'   
#'   # 6. Clean up the temporary directory
#'   unlink(temp_out_dir, recursive = TRUE)
#' }
#' }
#' @importFrom utils read.csv write.csv write.table
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ parLapply stopCluster mclapply
#' @export
run_adaptive_sampling <- function(initial_samples_file,
                                  scenario_name,
                                  distance_matrix,
                                  num_parallel_jobs = 5,
                                  max_cores = NULL,
                                  num_samples = 10,
                                  mapping_max_iter = 1000, 
                                  relative_epsilon = 1e-4,
                                  folds = 20,
                                  output_dir,
                                  verbose = FALSE) {
  # --- CAPTURE ORIGINAL COLUMN ORDER ---
  orig_cols <- names(read.csv(initial_samples_file, 
                              stringsAsFactors = FALSE, 
                              nrows = 1))
  # Parameter names
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  iterations <- ceiling(num_samples / num_parallel_jobs)

  # --- Input Validation ---
  if (missing(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("'output_dir' must be a single character string specifying the directory.", call. = FALSE)
  }
  # Validate num_samples
  if (!is.numeric(num_samples) || num_samples < 1 || num_samples != round(num_samples)) {
    stop("num_samples must be a positive integer")
  }
  for (p in c("mapping_max_iter", "folds", "num_parallel_jobs", "iterations")) {
    v <- get(p)
    if (!is.numeric(v) || v < 1 || v != round(v)) stop(sprintf("%s must be a positive integer", p))
  }
  
  if (verbose) {
    cat(sprintf("Calculating iterations per job: %d samples / %d jobs = %d iterations per job\n", 
                num_samples, num_parallel_jobs, iterations))
    cat(sprintf("This will produce approximately %d total new samples\n", 
                iterations * num_parallel_jobs))
  }
  
  if(!is.matrix(distance_matrix)) {
    stop("distance_matrix must be a matrix")
  }
  if (!is.numeric(relative_epsilon) || relative_epsilon <= 0) stop("relative_epsilon must be positive")
  if (!is.character(scenario_name) || length(scenario_name) != 1) stop("scenario_name must be a single string")

  # Determine available cores
  available_cores <- parallel::detectCores()
  if (is.null(max_cores)) {
    # Use all cores minus 1 to avoid bogging down the system
    max_cores <- max(1, available_cores - 1)
  } else {
    # Validate max_cores
    if (!is.numeric(max_cores) || max_cores < 1 || max_cores != round(max_cores)) {
      stop("max_cores must be a positive integer")
    }
    # Limit to available cores
    max_cores <- min(max_cores, available_cores)
  }
  
  # Setup directories
  adaptive_dir <- file.path(output_dir, "adaptive_sampling_jobs")
  param_dir    <- file.path(output_dir, "model_parameters")

  for (dir in c(adaptive_dir, param_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # --- "COMPREHENSIVE CLEANUP"  ---
  if (verbose) cat("Cleaning adaptive sampling jobs directory...\n")

  # List ALL files in the directory
  files_to_remove <- list.files(adaptive_dir, full.names = TRUE, recursive = TRUE)
  if (length(files_to_remove) > 0) file.remove(files_to_remove)

  # Validate initial samples and prepare master file
  # Check initial samples
  if (!file.exists(initial_samples_file)) stop("initial_samples_file not found: ", initial_samples_file)
  
  init <- read.csv(initial_samples_file, stringsAsFactors=FALSE)
  req <- c(par_names, "NLL")
  if (!all(req %in% names(init))) stop("Missing columns in initial samples: ", paste(setdiff(req, names(init)), collapse=", "))
  init <- init[complete.cases(init[, req]) & apply(init[, req], 1, function(x) all(is.finite(x))), ]
  if (nrow(init)==0) stop("No valid initial samples after filtering")
  results_file <- file.path(param_dir, paste0(scenario_name, "_model_parameters.csv"))

  # Copy initial samples to final destination if it doesn't exist
  if (!file.exists(results_file)) file.copy(initial_samples_file, results_file)

  # --- LOCAL EXECUTION PATH ---
  if (verbose) {
    cat(sprintf("Running %d parallel jobs with maximum %d cores\n", 
              num_parallel_jobs, max_cores))
  }

  # Create temporary files for each parallel job to write to
  make_temp <- function(i) file.path(adaptive_dir, sprintf("job_%02d_%s.csv", i, scenario_name))
  temps <- vapply(seq_len(num_parallel_jobs), make_temp, FUN.VALUE = "")
  for (i in seq_along(temps)) file.copy(initial_samples_file, temps[i], overwrite = TRUE)

  # ---------------- Local parallel execution ----------------
  # Launch parallel runs
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(min(num_parallel_jobs, max_cores))
    parallel::clusterExport(cl, c("adaptive_MC_sampling", "temps", "distance_matrix",
                                "mapping_max_iter", "relative_epsilon", "folds",
                                "output_dir", "scenario_name", "iterations"), envir = environment())
    parallel::clusterEvalQ(cl, library(topolow))
    parallel::parLapply(cl, seq_along(temps), function(i) {
      adaptive_MC_sampling(
        samples_file    = temps[i],
                          distance_matrix = distance_matrix,
                          iterations      = iterations,
                          batch_size      = 1,
                          mapping_max_iter= mapping_max_iter,
                          relative_epsilon= relative_epsilon,
                          folds           = folds,
                          num_cores       = 1,
                          scenario_name   = scenario_name,
                          #output_dir      = output_dir,
                          verbose         = TRUE)
    })
    parallel::stopCluster(cl)
  } else {
    parallel::mclapply(seq_along(temps), function(i) {
        adaptive_MC_sampling(
            samples_file    = temps[i],
            distance_matrix = distance_matrix,
            iterations      = iterations,
            batch_size      = 1,
            mapping_max_iter= mapping_max_iter,
            relative_epsilon= relative_epsilon,
            folds           = folds,
            num_cores       = 1,
            scenario_name   = scenario_name,
            #output_dir      = output_dir,
            verbose         = FALSE # Keep this TRUE for debugging
        )
    }, mc.cores = min(num_parallel_jobs, max_cores))
  }

  # Gather local results
  init2 <- read.csv(initial_samples_file, stringsAsFactors = FALSE)
  n0   <- nrow(init2)
  new_list <- lapply(temps, function(f) {
    if (file.exists(f)) {
        df <- tryCatch(read.csv(f, stringsAsFactors=FALSE), error=function(e) NULL)
        if (!is.null(df) && nrow(df)>n0) df[(n0+1):nrow(df), , drop=FALSE] else NULL
    } else {
        NULL
    }
  })
  all <- do.call(rbind, c(list(init2), new_list))
  # --- REORDER TO ORIGINAL INPUT HEADER ---
  all <- all[, orig_cols, drop=FALSE]
  write.csv(all, results_file, row.names=FALSE)
  file.remove(temps)
  if (verbose) cat("Local parallel jobs complete; results in", results_file, "\n")
  return(invisible(NULL)) # Exit after local run is complete
}



#' Sample from Weighted Distribution 
#'
#' @description
#' Generates new samples from a multivariate normal distribution with mean and 
#' covariance determined by weighting previous samples based on their likelihoods.
#'
#' @param dist List with mean and covariance matrix from previous samples
#' @param n Integer number of samples to generate 
#' @param epsilon Numeric probability of sampling from wider distribution
#' 
#' @importFrom stats runif
#' @importFrom MASS mvrnorm
#' @return Data frame of n new samples with columns for each parameter
#' @keywords internal
sample_from_distribution <- function(dist, n, epsilon) {
  # Validate inputs
  if (!is.list(dist) || !all(c("mean", "cov") %in% names(dist))) {
    stop("dist must be a list with mean and cov components")
  }
  
  if (!is.numeric(dist$mean)) {
    stop("dist$mean must be numeric")
  }
  
  if (!is.matrix(dist$cov)) {
    stop("dist$cov must be a matrix")
  }
  
  if (length(dist$mean) != nrow(dist$cov)) {
    stop("mean vector and covariance matrix dimensions must match")
  }
  
  if (!isSymmetric(dist$cov)) {
    stop("covariance matrix must be symmetric")
  }
  
  if (!is.numeric(n) || n < 1 || n != round(n)) {
    stop("n must be a positive integer")
  }
  
  if (!is.numeric(epsilon) || epsilon < 0 || epsilon > 1) {
    stop("epsilon must be between 0 and 1")
  }

  # epsilon % of time generate a random sample from a wider distribution
  if (runif(1) <= epsilon) {
    # Create diagonal matrix from diagonal of original covariance
    diag_elements <- diag(dist$cov)
    # Multiply diagonal elements by scalar 
    new_diag <- diag_elements * 2
    # Replace diagonal in original covariance matrix
    new_cov_matrix <- dist$cov
    diag(new_cov_matrix) <- new_diag
    samples <- mvrnorm(n, mu = dist$mean, Sigma = new_cov_matrix)
  } else {
    samples <- mvrnorm(n, mu = dist$mean, Sigma = dist$cov)
  }
  
  if(n==1){
    samples <- as.data.frame(t(samples))
  } else {
    samples <- as.data.frame(samples)
  }
  
  names(samples) <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  return(samples)
}


#' Generate New Parameter Samples Using KDE
#'
#' @description
#' Generates new parameter samples using weighted kernel density estimation
#' for each parameter independently.
#'
#' @param samples Data frame of previous samples with parameters and NLL
#' @param n Integer number of samples to generate
#' @param epsilon Numeric probability of wider bandwidth sampling
#' 
#' @importFrom stats na.omit runif bw.nrd0 approx
#' @return Data frame of n new samples
#' @keywords internal
generate_kde_samples <- function(samples, n, epsilon = 0) {
  # First, Remove outliers
  samples <- as.data.frame(lapply(samples, clean_data, k = 3))
  samples <- na.omit(samples)
  
  # Calculate weights from likelihoods
  log_likes <- -samples$NLL  
  std_log_likes <- log_likes - min(log_likes) + 0.05
  weights <- std_log_likes / sum(std_log_likes)
  
  # Get parameter names
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  
  # Initialize results
  new_samples <- data.frame(matrix(nrow = n, ncol = length(par_names)))
  names(new_samples) <- par_names
  
  # Generate samples for each parameter independently
  for (param in par_names) {
    # Generate wider bandwidth samples with probability epsilon
    if (runif(1) < epsilon) {
      bandwidth <- bw.nrd0(samples[[param]]) * 2
    } else {
      bandwidth <- NULL  # Use default
    }
    
    # Get kernel density estimate
    kde <- weighted_kde(samples[[param]], weights)
    
    # Sample from KDE using simple inverse transform sampling
    u <- runif(n)
    cdf <- cumsum(kde$y) / sum(kde$y)
    
    # Get samples through interpolation
    new_samples[[param]] <- approx(cdf, kde$x, u, rule = 2)$y
  }
  
  return(new_samples)
}


#' Weighted Kernel Density Estimation
#'
#' @description
#' Performs weighted kernel density estimation for univariate data. Uses parallel processing 
#' for efficiency. Useful for analyzing parameter distributions with importance weights.
#'
#' @param x Numeric vector of samples
#' @param weights Numeric vector of weights
#' @param n Integer number of evaluation points
#' @param from,to Numeric range for evaluation points
#' @return List containing:
#'   \item{x}{Vector of evaluation points}
#'   \item{y}{Vector of density estimates}
#' 
#' @importFrom stats sd dnorm
#' @importFrom parallel detectCores mclapply
#' @export
weighted_kde <- function(x, weights, n = 512, from = min(x), to = max(x)) {
  # Normalize weights
  weights <- weights / sum(weights)

  # Calculate bandwidth (Silverman's rule)
  bw <- 1.06 * sd(x) * length(x)^(-1/5)
  eval_points <- seq(from, to, length.out = n)
  
  compute_density <- function(z) {
    sum(weights * dnorm(z, mean = x, sd = bw))
  }

  # --- SAFE PARALLELISM CONTROL (CRAN Compliant) ---
  # In non-interactive sessions (like R CMD check), limit cores to 2.
  # Otherwise, respect the user's 'mc.cores' option or default to a safe number.
  num_cores <- if (!interactive()) {
    2L
  } else {
    getOption("mc.cores", max(1L, parallel::detectCores(logical = FALSE) - 1L))
  }

  # Use mclapply on non-Windows; fall back to sequential sapply on Windows
  if (.Platform$OS.type != "windows" && num_cores > 1) {
    density_est <- unlist(parallel::mclapply(eval_points, compute_density, mc.cores = num_cores))
  } else {
    density_est <- sapply(eval_points, compute_density)
  }
  # --- END OF SAFE PARALLELISM ---

  list(x = eval_points, y = density_est)
}


#' Unweighted Kernel Density Estimation 
#'
#' @description
#' Standard kernel density estimation for univariate data with various bandwidth
#' selection rules.
#'
#' @param x Numeric vector of samples
#' @param n Integer number of evaluation points
#' @param from,to Numeric range for evaluation points
#' @param bw Bandwidth selection ("nrd0", "nrd", "ucv", "bcv", "sj" or numeric)
#' @return List containing:
#'   \item{x}{Vector of evaluation points}
#'   \item{y}{Vector of density estimates}
#'   \item{bw}{Selected bandwidth}
#'   
#' @importFrom stats bw.nrd0 bw.nrd bw.ucv bw.bcv bw.SJ dnorm
#' @export
unweighted_kde <- function(x, n = 512, from = min(x), to = max(x), bw = "nrd0") {
  # Determine bandwidth
  if (is.character(bw)) {
    bw <- switch(bw,
                 nrd0 = bw.nrd0(x),
                 nrd = bw.nrd(x),
                 ucv = bw.ucv(x),
                 bcv = bw.bcv(x),
                 sj = bw.SJ(x),
                 stop("Unknown bandwidth rule")
    )
  }
  
  # Create evaluation points
  eval_points <- seq(from, to, length.out = n)
  
  # Compute KDE
  density_est <- sapply(eval_points, function(z) {
    mean(dnorm(z, mean = x, sd = bw))
  })
  
  list(x = eval_points, y = density_est, bw = bw)
}


#' Evaluate Likelihood with Cross-Validation 
#'
#' @description
#' Calculates cross-validated likelihood for a set of parameters by:
#' 1. Splitting data into training/validation sets
#' 2. Fitting model on training data
#' 3. Evaluating likelihood on validation set
#' 4. Repeating across folds
#' To calculate one NLL per set of parameters, the function uses a pooled errors approach which combines
#' all validation errors into one set, then calculate a single NLL. This approach has two main advantages:
#' 1- It treats all validation errors equally, respecting the underlying error distribution assumption
#' 2- It properly accounts for the total number of validation points
#' 
#' @param distance_matrix Distance matrix to fit
#' @param mapping_max_iter Maximum map optimization iterations
#' @param relative_epsilon Convergence threshold
#' @param N Number of dimensions
#' @param k0 Initial spring constant
#' @param cooling_rate Spring constant decay rate
#' @param c_repulsion Repulsion constant
#' @param folds Number of CV folds
#' @param num_cores Number of cores for parallel processing
#' @return List with:
#'   \item{Holdout_MAE}{Mean absolute error on validation data}
#'   \item{NLL}{Negative log likelihood}
#'   
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply
#' @keywords internal
likelihood_function <- function(distance_matrix, mapping_max_iter,
                              relative_epsilon, N, k0, cooling_rate,
                              c_repulsion, folds = 20, num_cores = 1) {
  # Create n folds in the data
  truth_matrix <- distance_matrix
  
  matrix_list <- replicate(folds, 
                          list(truth_matrix, NULL),
                          simplify = FALSE)
  
  num_elements <- sum(!is.na(distance_matrix))
  
  holdout_size <- floor(num_elements/(folds*2)) # 2 is because for each [i,j] we also null the [j,i] element
  
  # To cover n folds randomly, create a copy to remove fractions from it until finished:
  D_train <- distance_matrix
  
  for(i in 1:folds) {
    random_indices <- sample(which(!is.na(D_train)), size=holdout_size)
    input_matrix <- distance_matrix
    for(index in random_indices) {
      row <- (index - 1) %/% nrow(distance_matrix) + 1
      col <- (index - 1) %% ncol(distance_matrix) + 1
      input_matrix[row, col] <- NA
      input_matrix[col, row] <- NA
    }
    
    matrix_list[[i]][[2]] <- input_matrix
    
    for(index in random_indices) {
      row <- (index - 1) %/% nrow(D_train) + 1
      col <- (index - 1) %% ncol(D_train) + 1
      D_train[row, col] <- NA
      D_train[col, row] <- NA
    }
  }
  
  # Define the function for processing each fold - collecting error information for pooling
  process_sample <- function(i) {
    truth_matrix <- matrix_list[[i]][[1]]
    input_matrix <- matrix_list[[i]][[2]]
    
    # Call create_topolow_map with current parameters
    tryCatch({
      res_train <- create_topolow_map(
        distance_matrix = input_matrix, 
        ndim = N, 
        mapping_max_iter = mapping_max_iter, 
        k0 = k0, 
        cooling_rate = cooling_rate, 
        c_repulsion = c_repulsion,
        relative_epsilon = relative_epsilon,
        convergence_counter = 5,
        initial_positions = NULL,
        write_positions_to_csv = FALSE,
        verbose = FALSE
      )
      
      # Calculate errors
      p_dist_mat <- res_train$est_distances
      p_dist_mat <- as.matrix(p_dist_mat)
      
      errors <- error_calculator_comparison(p_dist_mat, truth_matrix, input_matrix)
      df <- errors$report_df
      
      # Extract actual errors for pooling
      out_sample_errors <- df$OutSampleError[!is.na(df$OutSampleError)]
      n_samples <- length(out_sample_errors)
      sum_abs_errors <- sum(abs(out_sample_errors))
      
      # Calculate fold-specific MAE for reference
      mae_holdout <- if(n_samples > 0) sum_abs_errors / n_samples else NA
      
      # Return fold results with information needed for pooling
      data.frame(
        Holdout_MAE = mae_holdout,
        n_samples = n_samples,
        sum_abs_errors = sum_abs_errors
      )
    }, error = function(e) {
      # Return NA result on error
      data.frame(
        Holdout_MAE = NA,
        n_samples = 0,
        sum_abs_errors = 0
      )
    })
  }
  
  # Process each fold - choose appropriate parallelization approach
  if(num_cores > 1) {
    if(.Platform$OS.type == "windows") {
      # For Windows, use a temporary cluster
      cl <- parallel::makeCluster(min(num_cores, folds))
      on.exit(parallel::stopCluster(cl))
      
      # Export required variables
      parallel::clusterExport(cl, 
                             c("matrix_list", "N", "k0", "cooling_rate", "c_repulsion", 
                               "mapping_max_iter", "relative_epsilon"), 
                             envir = environment())
      
      # Load packages
      parallel::clusterEvalQ(cl, {
        library(topolow)
      })
      
      # Process folds in parallel
      res_list <- parallel::parLapply(cl, 1:folds, process_sample)
    } else {
      # For Unix-like systems, use mclapply
      res_list <- parallel::mclapply(1:folds, process_sample, 
                                   mc.cores = min(num_cores, folds))
    }
  } else {
    # Sequential processing
    res_list <- lapply(1:folds, process_sample)
  }
  
  # Combine results with error handling
  valid_results <- !sapply(res_list, is.null) & 
                   !sapply(res_list, function(x) all(is.na(x$Holdout_MAE)))
                   
  if(sum(valid_results) == 0) {
    return(list(Holdout_MAE = NA, NLL = NA))
  }
  
  # Keep only valid results
  res_list <- res_list[valid_results]
  
  # Combine results
  res_df <- do.call(rbind, res_list)
  
  # Calculate pooled statistics
  total_samples <- sum(res_df$n_samples)
  total_abs_errors <- sum(res_df$sum_abs_errors)  
  # Calculate pooled MAE
  pooled_mae <- if(total_samples > 0) total_abs_errors / total_samples else NA
  
  # Calculate NLL using the correct formula and pooled MAE
  pooled_nll <- if(!is.na(pooled_mae)) total_samples*(1+log(2*pooled_mae)) else NA
  
  return(list(Holdout_MAE = pooled_mae, NLL = pooled_nll))
}


#' Perform Adaptive Monte Carlo Sampling
#'
#' @description
#' Main function implementing adaptive Monte Carlo sampling to explore parameter space.
#' Updates sampling distribution based on evaluated likelihoods. This is an internal
#' function called by `run_adaptive_sampling`.
#' 
#' @param samples_file Path to CSV with initial samples for this job.
#' @param distance_matrix Distance matrix to fit
#' @param iterations Number of sampling iterations per job
#' @param batch_size Samples per iteration (fixed to 1)
#' @param mapping_max_iter Maximum map optimization iterations 
#' @param relative_epsilon Convergence threshold
#' @param folds Number of CV folds
#' @param num_cores Number of cores for parallel processing
#' @param scenario_name Name for output files
#' @param verbose Logical. Whether to print progress messages. Default: FALSE
#'
#' @return A `data.frame` containing all samples (initial and newly generated) 
#' with their parameters and evaluated performance metrics. The data frame includes 
#' columns for the log-transformed parameters, `Holdout_MAE`, and `NLL`. 
#' Returns `NULL` if the results file was not created.
#' 
#' @importFrom filelock lock unlock
#' @importFrom utils read.csv write.table
#' @importFrom stats na.omit
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply
#' @export
adaptive_MC_sampling <- function(samples_file, 
                                 distance_matrix,
                                 iterations = 1, 
                                 batch_size = 1, # Now fixed to 1 by design
                                 mapping_max_iter, 
                                 relative_epsilon,
                                 folds = 20, 
                                 num_cores = 1,
                                 scenario_name, 
                                 verbose = FALSE) {
  # Require filelock for safe concurrent writes
  if (!requireNamespace("filelock", quietly = TRUE)) {
    stop("Package 'filelock' is required for safe writes. Please install it.")
  }

  # Handle parallel processing setup
  use_parallelism <- num_cores > 1
  if (use_parallelism) {
    if (verbose) cat("Setting up parallel processing\n")
    if (.Platform$OS.type == "windows") {
      if (verbose) cat("Using parallel cluster for Windows\n")
      cl <- parallel::makeCluster(num_cores)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl, c("distance_matrix", "mapping_max_iter", 
                                   "relative_epsilon", "folds"),
                             envir = environment())
      parallel::clusterEvalQ(cl, { library(topolow) })
    }
  }

  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  key_cols <- c(par_names, "Holdout_MAE", "NLL")
  
  for (iter in seq_len(iterations)) {
    if (verbose) cat(sprintf("\nStarting iteration %d of %d\n", iter, iterations))

    # Read current samples (always from samples_file)
    current_samples <- read.csv(samples_file, stringsAsFactors = FALSE)
    for (col in key_cols[key_cols %in% names(current_samples)]) {
      current_samples[[col]] <- as.numeric(as.character(current_samples[[col]]))
    }
    current_samples <- current_samples[apply(current_samples, 1, 
                                           function(row) all(is.finite(row))), ]
    current_samples <- na.omit(current_samples)

    if (nrow(current_samples) == 0) {
      warning("No valid samples remaining after filtering")
      break
    }

    # Burn-in and convergence checks
    if (nrow(current_samples) > 10) {
      burn_in <- min(round(nrow(current_samples) * 0.3), nrow(current_samples) - 5)
      current_samples <- current_samples[-seq_len(burn_in), ]
    }
    if (nrow(current_samples) > 500) {
      conv_check <- check_gaussian_convergence(
        data = current_samples[, par_names],
        window_size = 500,
        tolerance = 0.002
      )
      if (conv_check$converged) {
        if (verbose) cat("Convergence achieved at iteration", iter, "\n")
        break
      }
    }

    # Generate new samples via KDE
    new_samples <- generate_kde_samples(samples = current_samples, n = batch_size)
    # Ensure numeric columns
    for (col in par_names) {
      new_samples[[col]] <- as.numeric(new_samples[[col]])
    }
    
    # Define likelihood evaluation
    evaluate_sample <- function(i) {
      N <- round(exp(new_samples[i, "log_N"]))
      k0 <- exp(new_samples[i, "log_k0"])
      cooling_rate <- exp(new_samples[i, "log_cooling_rate"])
      c_repulsion <- exp(new_samples[i, "log_c_repulsion"])
      inner_cores <- if (use_parallelism) 1 else min(folds, num_cores)
      tryCatch({
        likelihood_function(
          distance_matrix = distance_matrix,
          mapping_max_iter = mapping_max_iter,
          relative_epsilon = relative_epsilon,
          N = N,
          k0 = k0,
          cooling_rate = cooling_rate,
          c_repulsion = c_repulsion,
          folds = folds,
          num_cores = inner_cores
        )
      }, error = function(e) {
        if (verbose) cat("Error in likelihood calculation:", e$message, "\n")
        NA
      })
    }
    
    # Evaluate samples with appropriate parallel method
    if (use_parallelism) {
      if (.Platform$OS.type == "windows") {
        parallel::clusterExport(cl, c("evaluate_sample", "new_samples"), envir = environment())
        new_likelihoods <- parallel::parLapply(cl, seq_len(nrow(new_samples)), evaluate_sample)
      } else {
        new_likelihoods <- parallel::mclapply(seq_len(nrow(new_samples)), evaluate_sample, mc.cores = num_cores)
      }
    } else {
      # Sequential processing
      new_likelihoods <- lapply(seq_len(nrow(new_samples)), evaluate_sample)
    }
    
    # Filter valid results
    valid_results <- !sapply(new_likelihoods, is.null) &
                     !sapply(new_likelihoods, function(x) all(is.na(unlist(x))))

    if (any(valid_results)) {
      likelihoods_mat <- do.call(rbind, lapply(new_likelihoods[valid_results], unlist))
      new_likelihoods_df <- as.data.frame(likelihoods_mat)
      colnames(new_likelihoods_df) <- c("Holdout_MAE", "NLL")

      valid_samples <- new_samples[valid_results, , drop = FALSE]
      valid_samples$Holdout_MAE <- as.numeric(new_likelihoods_df$Holdout_MAE)
      valid_samples$NLL <- as.numeric(new_likelihoods_df$NLL)

      for (col in key_cols[key_cols %in% names(valid_samples)]) {
        valid_samples[[col]] <- as.numeric(as.character(valid_samples[[col]]))
      }
      valid_samples <- valid_samples[apply(valid_samples, 1, function(row) all(is.finite(row))), ]

      if (nrow(valid_samples) > 0) {
        result_file <- samples_file
        lock_path <- paste0(result_file, ".lock")
        lock <- filelock::lock(lock_path)

        if (!file.exists(result_file)) {
          write.table(valid_samples, result_file,
                      sep = ",", row.names = FALSE, col.names = TRUE,
                      append = FALSE, quote = TRUE,
                      eol = "\n")
        } else {
          write.table(valid_samples, result_file,
                      sep = ",", row.names = FALSE, col.names = FALSE,
                      append = TRUE, quote = TRUE,
                      eol = "\n")
        }
        filelock::unlock(lock)

        if (verbose) cat(sprintf("Safely appended %d new valid samples to %s\n", 
                                 nrow(valid_samples), result_file))
      } else if (verbose) {
        cat("No valid samples in this iteration\n")
      }
    } else if (verbose) {
      cat("All likelihood evaluations failed in this iteration\n")
    }
  }

  # Return final samples
  result_file <- samples_file
  if (file.exists(result_file)) {
    final_samples <- read.csv(result_file, stringsAsFactors = FALSE)
    return(final_samples)
  } else {
    warning("No results file created")
    return(NULL)
  }
}



#' Calculate Weighted Marginal Distributions
#'
#' @description
#' Calculates marginal distributions for each parameter with weights derived from 
#' log-likelihoods.
#'
#' @param samples Data frame containing:
#'        - log_N, log_k0, log_cooling_rate, log_c_repulsion: Parameter columns
#'        - NLL: Negative log-likelihood column
#' @return Named list of marginal distributions, each containing:
#'   \item{x}{Vector of parameter values}
#'   \item{y}{Vector of density estimates}
#' @details 
#' Uses kernel density estimation weighted by normalized likelihoods.
#'
#' @importFrom stats na.omit
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply detectCores
#' @export
calculate_weighted_marginals <- function(samples) {
  # Input validation
  required_cols <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion", "NLL")
  if (!all(required_cols %in% names(samples))) {
    stop("Missing required columns: ", 
         paste(setdiff(required_cols, names(samples)), collapse = ", "))
  }
  
  if (!all(sapply(samples[required_cols], is.numeric))) {
    stop("All parameters must be numeric")
  }
  
  # Validate NLL values
  if (all(is.infinite(samples$NLL))) {
    stop("All NLL values are infinite")
  }
  
  if (any(is.na(samples$NLL))) {
    warning("NA values in NLL column will be removed")
    samples <- samples[!is.na(samples$NLL), ]
  }
  
  if (nrow(samples) == 0) {
    stop("No valid samples remain after filtering")
  }
  
  samples <- as.data.frame(lapply(samples, clean_data, k = 3))
  samples <- na.omit(samples)

  # Calculate weights from log-likelihoods
  log_likelihoods <- -samples$NLL
  std_log_likelihoods <- log_likelihoods - min(log_likelihoods) + 0.05
  weights <- std_log_likelihoods / sum(std_log_likelihoods)

  # Define parameter columns to process
  vars <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")

  # Use a simple, sequential loop. The parallelism is now handled safely
  # inside weighted_kde, which this loop calls four times.
  # This completely avoids the nested parallelism problem.
  marginals <- lapply(vars, function(var) {
    weighted_kde(samples[[var]], weights = weights)
  })

  # Set names and return
  names(marginals) <- vars
  return(marginals)
}


#' Find Mode of Density Distribution
#'
#' @description
#' Calculates the mode (maximum point) of a kernel density estimate.
#'
#' @param density List containing density estimate with components:
#'   \describe{
#'     \item{x}{Vector of values}
#'     \item{y}{Vector of density estimates}
#'   }
#' @return Numeric value of the mode
#' @export
find_mode <- function(density) {
  if (!is.list(density) || !all(c("x", "y") %in% names(density))) {
    stop("density must be a list with x and y components")
  }
  
  if (length(density$x) != length(density$y)) {
    stop("x and y components must have same length")
  }
  
  if (!all(is.finite(c(density$x, density$y)))) {
    stop("density values must be finite")
  }
  
  density$x[which.max(density$y)]
}


#' Create Grid Around Maximum Likelihood Estimate
#'
#' @description
#' Generates a sequence of values centered on the maximum likelihood estimate (MLE)
#' of a parameter, extending by specified factors in each direction.
#'
#' @param samples Data frame of MCMC samples with NLL column
#' @param param Character name of parameter column
#' @param num_points Integer number of points in grid
#' @param start_factor Numeric factor for lower bound relative to MLE
#' @param end_factor Numeric factor for upper bound relative to MLE
#' @return Numeric vector of grid points
#' @keywords internal
get_grid <- function(samples, param, num_points, start_factor, end_factor) {
  if (!is.data.frame(samples) || !all(c(param, "NLL") %in% names(samples))) {
    stop("samples must be a data frame containing param and NLL columns")
  }
  
  if (!is.numeric(samples[[param]])) {
    stop("Parameter column must be numeric")
  }
  
  if (start_factor >= end_factor) {
    stop("start_factor must be less than end_factor")
  }
  
  # Convert NLL to LL and calculate weighted marginals
  samples$LL <- -samples$NLL
  weighted_marginals <- calculate_weighted_marginals(samples)
  
  # Find MLE and create grid
  MLE <- find_mode(weighted_marginals[[param]])
  
  seq(start_factor * MLE, end_factor * MLE, length.out = num_points)
}



#' Profile Likelihood Analysis Results Class
#'
#' @description
#' S3 class for storing and manipulating profile likelihood analysis results.
#' Includes parameter values, log-likelihoods, and metadata for visualization.
#'
#' @keywords internal
profile_likelihood_result <- function(param_values, ll_values, param_name, 
                                    bandwidth, sample_counts) {
  structure(
    list(
      param = param_values,
      ll = ll_values,
      param_name = param_name,
      bandwidth = bandwidth,
      sample_counts = sample_counts
    ),
    class = "profile_likelihood"
  )
}

#' Profile Likelihood Analysis
#'
#' @description
#' Calculates profile likelihood for a parameter by evaluating conditional maximum 
#' likelihood across a grid of parameter values. Uses local sample windowing to
#' estimate conditional likelihoods. This implementation is not a classical profile likelihood 
#' calculation, but rather an "empirical profile likelihood" which estimates the profile 
#' likelihood at each point based on the many observations previously sampled in Monte Carlo simulations.
#'
#' @details
#' For each value in the parameter grid, the function:
#' 1. Identifies nearby samples using bandwidth window
#' 2. Calculates conditional maximum likelihood from these samples
#' 3. Tracks sample counts to assess estimate reliability
#' 4. Handles boundary conditions and sparse regions
#'
#' @param param Character name of parameter to analyze
#' @param samples Data frame containing parameter samples and log-likelihoods
#' @param grid_size Integer number of grid points (default: 48)
#' @param bandwidth_factor Numeric factor for local sample window (default: 0.03)
#' @param start_factor,end_factor Numeric range multipliers for parameter grid (default: 0.5, 1.2)
#' @param min_samples Integer minimum samples required for reliable estimate (default: 10)
#' @return Object of class "profile_likelihood" containing:
#'   \item{param}{Vector of parameter values}
#'   \item{ll}{Vector of log-likelihood values}
#'   \item{param_name}{Name of analyzed parameter}
#'   \item{bandwidth}{Bandwidth used for local windows}
#'   \item{sample_counts}{Number of samples per estimate}
#' @examples
#' # Create a sample data frame of MCMC samples
#' mcmc_samples <- data.frame(
#'   log_N = log(runif(50, 2, 10)),
#'   log_k0 = log(runif(50, 1, 5)),
#'   log_cooling_rate = log(runif(50, 0.01, 0.1)),
#'   log_c_repulsion = log(runif(50, 0.1, 1)),
#'   NLL = runif(50, 20, 100)
#' )
#' 
#' # Calculate profile likelihood for parameter "log_N"
#' pl <- profile_likelihood("log_N", mcmc_samples, 
#'                         grid_size = 10, # Smaller grid for a quick example
#'                         bandwidth_factor = 0.05)
#'                         
#' # Print the results
#' print(pl)
#'
#' @seealso 
#' \code{\link{plot.profile_likelihood}} for visualization
#' 
#' @importFrom stats sd approx
#' @export
profile_likelihood <- function(param, samples, grid_size = 40, 
                             bandwidth_factor = 0.05,
                             start_factor = 0.5, end_factor = 1.5,
                             min_samples = 5) {
  
  # Input validation
  if (!is.character(param) || length(param) != 1) {
    stop("param must be a single character string")
  }
  if (!param %in% names(samples)) {
    stop(sprintf("Parameter '%s' not found in samples", param))
  }
  if (!all(c("NLL") %in% names(samples))) {
    stop("Samples must contain 'NLL' column")
  }
  if (grid_size < 2) {
    stop("grid_size must be at least 2")
  }
  if (bandwidth_factor <= 0) {
    stop("bandwidth_factor must be positive")
  }
  if (start_factor >= end_factor) {
    stop("start_factor must be less than end_factor")
  }
  
  # Get parameter grid
  grid_values <- get_grid(samples, param, grid_size, start_factor, end_factor)
  
  # Initialize results
  ll_values <- numeric(length(grid_values))
  sample_counts <- numeric(length(grid_values))
  
  # Calculate profile likelihood
  for (i in seq_along(grid_values)) {
    val <- grid_values[i]
    
    # Get conditional samples
    param_range <- sd(samples[[param]]) * bandwidth_factor
    conditional_samples <- samples[abs(samples[[param]] - val) <= param_range, ]
    n_samples <- nrow(conditional_samples)
    
    if (n_samples < min_samples) {
      warning(sprintf("Too few samples (%d) near %s = %.3f", 
                     n_samples, param, val))
      ll_values[i] <- NA
    } else {
      # Calculate maximum likelihood from top n=3 samples
      x <- -conditional_samples$NLL[order(conditional_samples$NLL)[1:min(3, n_samples)]]
      max_x <- max(x)
      ll_values[i] <- max_x + log(sum(exp(x - max_x))) - log(length(x))
    }
    sample_counts[i] <- n_samples
  }
  
  # Handle boundary conditions - interpolate NAs if possible
  na_idx <- which(is.na(ll_values))
  if (length(na_idx) > 0 && length(na_idx) < length(ll_values)) {
    ll_values[na_idx] <- approx(grid_values[-na_idx], 
                               ll_values[-na_idx],
                               grid_values[na_idx])$y
  }
  
  # Create result object
  result <- profile_likelihood_result(
    param_values = grid_values,
    ll_values = ll_values,
    param_name = param,
    bandwidth = param_range,
    sample_counts = sample_counts
  )
  
  return(result)
}


#' Plot Method for Profile Likelihood Objects
#'
#' @description
#' Creates a visualization of profile likelihood for a parameter showing maximum
#' likelihood estimates and confidence intervals. Supports mathematical notation
#' for parameter names and configurable output settings.
#' 
#' Confidence interval is found using the likelihood ratio test: 
#' \eqn{LR(\theta_{ij}) = -2[log L_{max}(\theta_{ij}) - log L_{max}(\hat{\theta})]}
#' where \eqn{\hat{\theta}} is the maximum likelihood estimate for all parameters.
#' The 95% confidence interval is:
#' \eqn{\{\theta_{ij} : LR(\theta_{ij}) \leq \chi^2_{1,0.05} = 3.84\}}
#'
#' @param x A profile_likelihood object
#' @param LL_max Numeric maximum log-likelihood value
#' @param width Numeric width of output plot in inches (default: 3.5)
#' @param height Numeric height of output plot in inches (default: 3.5)
#' @param save_plot Logical. Whether to save plot to file. Default: FALSE
#' @param output_dir Character. Directory for output files. Required if `save_plot` is TRUE.
#' @param ... Additional arguments passed to plot
#' @return A ggplot object
#' @examples
#' # These examples take more than 5 seconds to run, so they are not run by default. 
#' # Use parallel processing (the default) to speed up.
#' \donttest{
#' # Create a sample data frame of MCMC samples
#' samples <- data.frame(
#'   log_N = log(runif(50, 2, 10)),
#'   log_k0 = log(runif(50, 1, 5)),
#'   log_cooling_rate = log(runif(50, 0.01, 0.1)),
#'   log_c_repulsion = log(runif(50, 0.1, 1)),
#'   NLL = runif(50, 20, 100)
#' )
#'
#' # Calculate profile likelihood
#' pl_result <- profile_likelihood("log_N", samples, grid_size = 10)
#' 
#' # Plot with maximum likelihood from samples
#' LL_max <- max(-samples$NLL)
#' # The plot function requires the ggplot2 package
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(pl_result, LL_max, width = 4, height = 3)
#' }
#' }
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_text labs theme_minimal theme element_text element_line element_blank element_rect margin scale_y_continuous ggsave
#' @importFrom scales label_number
#' @method plot profile_likelihood
#' @export 
plot.profile_likelihood <- function(x, LL_max, width = 3.5, height = 3.5,
                                    save_plot = FALSE, output_dir, ...) {
  # Convert profile likelihood object to data frame
  LL_list_param <- data.frame(
    param = x$param,
    LL = x$ll
  )
  
  # Convert parameter names to mathematical expressions
  param_expr <- switch(x$param_name,
                       "log_N" = expression(log(N)),
                       "log_c_repulsion" = expression(log(c)),
                       "log_cooling_rate" = expression(log(alpha)),
                       "log_k0" = expression(log(k)),
                       x$param_name)  # Default to original name if not matched
  
  # Create title expression
  title_expr <- switch(x$param_name,
                       "log_N" = expression(paste("Profile Likelihood: ", log(N))),
                       "log_c_repulsion" = expression(paste("Profile Likelihood: ", log(c))),
                       "log_cooling_rate" = expression(paste("Profile Likelihood: ", log(alpha))),
                       "log_k0" = expression(paste("Profile Likelihood: ", log(k))),
                       paste("Profile Likelihood:", x$param_name))
  
  CI_95_LL <- LL_max - 3.84 / 2
  
  p <- ggplot(LL_list_param, aes(x = .data$param, y = .data$LL)) +
    geom_line(color = "steelblue", size = 0.5) +
    geom_hline(yintercept = CI_95_LL, 
               linetype = "dashed", color = "black", size = 0.4) +
    geom_text(aes(x = min(.data$param), y = CI_95_LL + 0.02, 
                  label = "95% CI"),
              color = "black", vjust = -0.5, hjust = -0.05, size = 2.1) +
    labs(title = title_expr,
         x = param_expr,
         y = "Log Likelihood") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      plot.margin = margin(0.05, 0.05, 0.05, 0.05, "cm")
    ) +
    scale_y_continuous(labels = scales::label_number(big.mark = ""))
  
  if(save_plot) {
    if (missing(output_dir)) {
      stop("'output_dir' must be provided when save_plot is TRUE.", call. = FALSE)
    }
    filename <- file.path(output_dir, 
                          paste0("profile_likelihood_", x$param_name, ".pdf"))
    
    # Try to save using standard PDF device first (most compatible approach)
    tryCatch({
      ggsave_white_bg(filename, p, width = width, height = height, 
             device = "pdf", units = "in")
    }, error = function(e) {
      # If standard PDF fails, try without specifying device (system default)
      message("Standard PDF device failed. Trying system default device.")
      tryCatch({
        ggsave_white_bg(filename, p, width = width, height = height, units = "in")
      }, error = function(e) {
        message("Failed to save plot: ", e$message)
      })
    })
  }
  
  return(p)
}


#' Print Method for Profile Likelihood Objects
#'
#' @param x Profile likelihood object
#' @param ... Additional arguments passed to print
#' @method print profile_likelihood
#' @return The original `profile_likelihood` object (invisibly). Called for side effects (prints a summary to the console).

#' @export
print.profile_likelihood <- function(x, ...) {
  cat("Profile Likelihood Analysis\n")
  cat("Parameter:", x$param_name, "\n")
  cat("Grid points:", length(x$param), "\n")
  cat("Bandwidth:", round(x$bandwidth, 6), "\n")
  cat("Sample counts (min/median/max):",
      min(x$sample_counts), "/",
      median(x$sample_counts), "/", 
      max(x$sample_counts), "\n")
  invisible(x)
}


#' Parameter Sensitivity Analysis
#'
#' @description
#' Analyzes the sensitivity of model performance (MAE) to changes in a parameter.
#' Uses binning to identify the minimum MAE across parameter ranges and calculates
#' thresholds for acceptable parameter values.
#'
#' @details
#' The function performs these steps:
#' 1. Cleans the input data using MAD-based outlier detection
#' 2. Bins the parameter values into equal-width bins
#' 3. Calculates the minimum MAE within each bin. Analogous to "poorman's likelihood" approach, 
#' minimum MAE within each bin is an empirical estimate of the performance surface at this parameter
#' value when other parameters are at their optimal values.
#' 4. Identifies a threshold of acceptable performance (default: Topolow min. +5% MAE)
#' 5. Returns an object for visualization and further analysis
#'
#' @param param Character name of parameter to analyze
#' @param samples Data frame containing parameter samples and performance metrics
#' @param bins Integer number of bins for parameter range (default: 40)
#' @param mae_col Character name of column containing MAE values (default: "Holdout_MAE")
#' @param threshold_pct Numeric percentage above minimum for threshold calculation (default: 5)
#' @param min_samples Integer minimum number of samples required in a bin (default: 1)
#' @return Object of class "parameter_sensitivity" containing:
#'   \item{param_values}{Vector of parameter bin midpoints}
#'   \item{min_mae}{Vector of minimum MAE values per bin}
#'   \item{param_name}{Name of analyzed parameter}
#'   \item{threshold}{Threshold value (default: Topolow min. +5%)}
#'   \item{min_value}{Minimum MAE value across all bins}
#'   \item{sample_counts}{Number of samples per bin}
#'   
#' @importFrom stats na.omit
#' @importFrom graphics hist
#' @export
parameter_sensitivity_analysis <- function(param, samples, bins = 30, 
                                          mae_col = "Holdout_MAE",
                                          threshold_pct = 5,
                                          min_samples = 1) {
  # Validate inputs
  if (!is.character(param) || length(param) != 1) {
    stop("param must be a single character string")
  }
  if (!param %in% names(samples)) {
    stop(sprintf("Parameter '%s' not found in samples", param))
  }
  if (!mae_col %in% names(samples)) {
    stop(sprintf("MAE column '%s' not found in samples", mae_col))
  }
  if (bins < 2) {
    stop("bins must be at least 2")
  }
  
  # Clean the input data using MAD-based outlier detection
  # clean_samples <- as.data.frame(lapply(samples[, c(param, mae_col)], clean_data, k = 3))
  # clean_samples <- na.omit(clean_samples)
  clean_samples <- na.omit(samples)
  
  # Check if we have enough data
  if (nrow(clean_samples) < min_samples * 2) {
    stop("Insufficient data after cleaning")
  }
  
  # Find the range of parameter values to create bins
  min_param <- min(clean_samples[[param]])
  max_param <- max(clean_samples[[param]])
  
  # Create bins of equal width
  bin_width <- (max_param - min_param) / bins
  breaks <- seq(min_param, max_param + 1e-10, by = bin_width)  # Add a tiny amount to include max value
  
  # Use hist to get bin assignments (but don't plot)
  h <- hist(clean_samples[[param]], breaks = breaks, plot = FALSE)
  
  # For each bin, find the minimum MAE if there are data points
  bin_min_mae <- data.frame(
    bin_midpoint = h$mids,
    min_mae = NA,
    sample_count = h$counts
  )
  
  for (i in 1:bins) {
    # Find data points in this bin
    bin_data <- clean_samples[clean_samples[[param]] >= breaks[i] & 
                              clean_samples[[param]] < breaks[i+1], ]
    if (nrow(bin_data) >= min_samples) {
      bin_min_mae$min_mae[i] <- min(bin_data[[mae_col]])
    }
  }
  
  # Remove bins with insufficient data for analysis
  plot_data <- bin_min_mae[!is.na(bin_min_mae$min_mae), ]
  
  # Calculate min value and threshold
  min_value <- min(plot_data$min_mae)
  threshold_value <- min_value * (1 + threshold_pct/100)
  
  # Create result object
  result <- structure(
    list(
      param_values = plot_data$bin_midpoint,
      min_mae = plot_data$min_mae,
      param_name = param,
      threshold = threshold_value,
      min_value = min_value,
      sample_counts = plot_data$sample_count[!is.na(bin_min_mae$min_mae)]
    ),
    class = "parameter_sensitivity"
  )
  
  return(result)
}


#' Plot Method for Parameter Sensitivity Analysis
#'
#' @description
#' Creates a visualization of parameter sensitivity showing minimum MAE values
#' across parameter ranges with threshold indicators.
#'
#' @param x A parameter_sensitivity object
#' @param width Numeric width of output plot in inches (default: 3.5)
#' @param height Numeric height of output plot in inches (default: 3.5)
#' @param save_plot Logical. Whether to save plot to file. Default: FALSE
#' @param output_dir Character. Directory for output files. Required if `save_plot` is TRUE.
#' @param y_limit_factor Numeric. Factor to set the upper y-axis limit as a percentage above 
#'        the threshold value (e.g., 1.10 for 10% above). Default: NULL (automatic scaling)
#' @param ... Additional arguments passed to plot
#' @return A ggplot object
#' @method plot parameter_sensitivity
#' @importFrom ggplot2 ggplot aes geom_line geom_hline annotate labs theme_minimal theme element_text element_line element_blank element_rect margin scale_y_continuous ggsave
#' @importFrom scales comma label_number
#' @export
plot.parameter_sensitivity <- function(x, width = 3.5, height = 3.5,
                                     save_plot = FALSE, output_dir,
                                     y_limit_factor = NULL, ...) {
  # Convert to data frame for ggplot
  plot_data <- data.frame(
    param = x$param_values,
    min_mae = x$min_mae
  )
  
  # Convert parameter names to mathematical expressions
  param_expr <- switch(x$param_name,
                     "log_N" = expression(log(N)),
                     "log_c_repulsion" = expression(log(c)),
                     "log_cooling_rate" = expression(log(alpha)),
                     "log_k0" = expression(log(k)),
                     x$param_name)  # Default to original name if not matched
  
  # Create title expression
  title_expr <- switch(x$param_name,
                      "log_N" = expression(paste("Parameter Sensitivity: ", log(N))),
                      "log_c_repulsion" = expression(paste("Parameter Sensitivity: ", log(c))),
                      "log_cooling_rate" = expression(paste("Parameter Sensitivity: ", log(alpha))),
                      "log_k0" = expression(paste("Parameter Sensitivity: ", log(k))),
                      paste("Parameter Sensitivity:", x$param_name))
  
  # Calculate position for the threshold label
  label_x_pos <- min(plot_data$param) + 0.03 * (max(plot_data$param) - min(plot_data$param))
  
  # Create the plot
  p <- ggplot() +
    # Main data line
    geom_line(data = plot_data, 
              aes(x = .data$param, y = .data$min_mae),
              color = "steelblue",
              size = 0.8) +
    # Threshold line with direct label
    geom_hline(yintercept = x$threshold,
               linetype = "dashed", 
               color = "black", 
               size = 0.6) +
    # Add label for threshold line
    annotate("text", 
             x = label_x_pos, 
             y = x$threshold + 0.004, 
             label = "Threshold (min +5%)",
             hjust = 0, 
             size = 3.5,
             color = "black") +
    # Labels and theme
    labs(#title = title_expr,
         x = param_expr,
         y = "Validation MAE") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
    )
  
  # Set y-axis limits if y_limit_factor is provided
  if (!is.null(y_limit_factor)) {
    if (!is.numeric(y_limit_factor) || y_limit_factor <= 1) {
      warning("y_limit_factor must be a number greater than 1, using automatic scaling")
    } else {
      # Calculate upper limit based on threshold
      upper_limit <- x$threshold * y_limit_factor
      # Set lower limit to the minimum value or 0, whichever is smaller
      lower_limit <- min(plot_data$min_mae, na.rm = TRUE)*0.99
      
      p <- p + scale_y_continuous(
        limits = c(lower_limit, upper_limit),
        labels = scales::label_number(big.mark = "")
      )
    }
  } else {
    # Default continuous scale with comma formatting
    p <- p + scale_y_continuous(labels = scales::comma)
  }
  
  if(save_plot) {
    if (missing(output_dir)) {
      stop("'output_dir' must be provided when save_plot is TRUE.", call. = FALSE)
    }
    filename <- file.path(output_dir, 
                         paste0("parameter_sensitivity_", x$param_name, ".pdf"))
    
    # Use safe saving approach to handle device compatibility issues
    tryCatch({
      ggsave_white_bg(filename, p, width = width, height = height, 
             device = "pdf", units = "in")
    }, error = function(e) {
      message("Standard PDF device failed. Trying system default device.")
      tryCatch({
        ggsave_white_bg(filename, p, width = width, height = height, units = "in")
      }, error = function(e) {
        message("Failed to save plot: ", e$message)
      })
    })
  }
  
  return(p)
}



#' Print Method for Parameter Sensitivity Objects
#'
#' @param x A parameter_sensitivity object
#' @param ... Additional arguments passed to print
#' @method print parameter_sensitivity
#' @return No return value, called for side effects (prints a summary to the console).
#' @export
print.parameter_sensitivity <- function(x, ...) {
  cat("Parameter Sensitivity Analysis\n")
  cat("Parameter:", x$param_name, "\n")
  cat("Number of bins:", length(x$param_values), "\n")
  cat("Minimum MAE:", round(x$min_value, 6), "\n")
  cat("Threshold value (5%):", round(x$threshold, 6), "\n")
  cat("Sample counts (min/median/max):",
      min(x$sample_counts), "/",
      median(x$sample_counts), "/", 
      max(x$sample_counts), "\n")
}