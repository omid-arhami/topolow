# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# Licensed under MIT License
# R/adaptive_sampling.R

#' Latin Hypercube and Adaptive Monte Carlo Sampling Functions
#' 
#' @description
#' This file contains functions for performing Latin Hypercube and adaptive Monte Carlo 
#' sampling in parameter space. The AMC sampling adapts based on previous evaluations to focus
#' sampling in high-likelihood regions. The functions run either locally or via SLURM.
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
#' using k-fold cross-validation to assess prediction accuracy.
#'
#' @details
#' The function performs these steps:
#' 1. Generates LHS samples in parameter space 
#' 2. Creates k-fold splits of input data
#' 3. For each parameter set and fold:
#'    - Trains model on training set
#'    - Evaluates on validation set
#'    - Calculates MAE and negative log likelihood
#' 4. Can run computation locally or distribute via SLURM
#'
#' Parameters ranges are transformed to log scale where appropriate to handle
#' different scales effectively.
#'
#' @param distance_matrix Matrix or data frame. Input distance matrix. Must be square
#'        and symmetric. Can contain NA values for missing measurements.
#' @param max_iter Integer. Maximum number of optimization iterations.
#' @param relative_epsilon Numeric. Convergence threshold for relative change in error.
#' @param convergence_counter Integer. Number of iterations below threshold before 
#'        declaring convergence.
#' @param scenario_name Character. Name for output files and job identification.
#' @param N_min,N_max Integer. Range for number of dimensions parameter.
#' @param k0_min,k0_max Numeric. Range for initial spring constant parameter.
#' @param cqq_min,cqq_max Numeric. Range for repulsion constant parameter.
#' @param k_decay_min,k_decay_max Numeric. Range for spring decay parameter.
#' @param num_samples Integer. Number of LHS parameter samples to evaluate.
#' @param folds Integer. Number of cross-validation folds. Default: 10.
#' @param verbose Logical. Whether to print progress messages. Default: FALSE.
#' @param write_files Logical. Whether to save results to CSV. Default: TRUE.
#' @param output_dir Character. Directory where output and temporary files will be saved. If NULL,
#'        uses current working directory. Directory will be created if it doesn't exist.
#' @param num_cores Integer. Number of CPU cores to use for parallel processing.
#'        Default: 1.
#' @param use_slurm Logical. Whether to submit jobs via SLURM. Default: FALSE.
#' @param cider Logical. Whether to use cider queue in SLURM. Default: FALSE.
#'
#' @return If write_files=FALSE, returns a data frame with columns:
#'   \item{N}{Number of dimensions used}
#'   \item{k0}{Initial spring constant}
#'   \item{k_decay}{Spring decay rate}
#'   \item{cqq}{Repulsion constant}
#'   \item{Holdout_MAE}{Mean absolute error on validation sets}
#'   \item{NLL}{Negative log likelihood}
#'
#' If write_files=TRUE, results are saved to CSV files in the format:
#' "{scenario_name}_hyperparam_results.csv"
#'
#' @examples
#' \dontrun{
#' # Generate sample distance matrix
#' dist_mat <- matrix(runif(100), 10, 10)
#' dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
#' diag(dist_mat) <- 0
#'
#' # Run local optimization
#' results <- run_parameter_optimization(
#'   distance_matrix = dist_mat,
#'   max_iter = 1000,
#'   relative_epsilon = 1e-4,
#'   convergence_counter = 10,
#'   scenario_name = "test_opt",
#'   N_min = 2, N_max = 10,
#'   k0_min = 1, k0_max = 30,
#'   cqq_min = 0.00001, cqq_max = 0.2,
#'   k_decay_min = 0.00001, k_decay_max = 0.2,
#'   num_samples = 20,
#'   num_cores = 4
#' )
#'
#' # Run with SLURM
#' run_parameter_optimization(
#'   distance_matrix = dist_mat,
#'   max_iter = 1000,
#'   scenario_name = "slurm_opt",
#'   N_min = 2, N_max = 10,
#'   num_samples = 50,
#'   use_slurm = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{topolow_full}} for the core optimization algorithm
#'
#' @importFrom lhs maximinLHS
#' @importFrom parallel mclapply detectCores
#' @importFrom stats runif qunif
#' @export
run_parameter_optimization <- function(distance_matrix,
                                     max_iter,
                                     relative_epsilon,
                                     convergence_counter,
                                     scenario_name,
                                     N_min,
                                     N_max,
                                     k0_min,
                                     k0_max,
                                     cqq_min,
                                     cqq_max,
                                     k_decay_min,
                                     k_decay_max,
                                     num_samples,
                                     folds = 20,
                                     verbose = FALSE,
                                     write_files = TRUE,
                                     output_dir = NULL,
                                     num_cores = 1,
                                     use_slurm = FALSE,
                                     cider = FALSE) {
  # Input validation
  if (!is.matrix(distance_matrix)) {
    stop("distance_matrix must be a matrix")
  }
  if (nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix must be square")
  }
  
  # Validate integer parameters
  integer_params <- list(
    max_iter = max_iter,
    convergence_counter = convergence_counter,
    N_min = N_min,
    N_max = N_max,
    num_samples = num_samples,
    folds = folds,
    num_cores = num_cores
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
    cqq_min = cqq_min,
    cqq_max = cqq_max,
    k_decay_min = k_decay_min,
    k_decay_max = k_decay_max
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
  if (cqq_min >= cqq_max) {
    stop("cqq_min must be less than cqq_max")
  }
  if (k_decay_min >= k_decay_max) {
    stop("k_decay_min must be less than k_decay_max")
  }
  
  # Validate logical parameters
  logical_params <- list(
    verbose = verbose,
    write_files = write_files,
    use_slurm = use_slurm,
    cider = cider
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
  
  # Validate SLURM availability if requested
  if (use_slurm && !has_slurm()) {
    stop("SLURM requested but not available on this system")
  }
  
  # Check number of cores
  max_cores <- detectCores()
  if (num_cores > max_cores) {
    warning(sprintf("Requested %d cores but only %d available. Using %d cores.", 
                    num_cores, max_cores, max_cores))
    num_cores <- max_cores
  }
  
  # Handle output directory
  if (write_files || use_slurm) {
    # If no output directory specified, use current working directory
    if (is.null(output_dir)) {
      output_dir <- getwd()
    }
  
    # Create required directories
    hyperparam_dir <- file.path(output_dir, "hyperparam_results")
    run_topolow_dir <- file.path(output_dir, "run_topolow_arguments")
    
    if (write_files) {
      for (dir in c(hyperparam_dir, run_topolow_dir)) {
        if (!dir.exists(dir)) {
          dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        }
      }
    }
    
    # Verify write permissions
    if (write_files) {
      test_file <- file.path(hyperparam_dir, "test_write.txt")
      tryCatch({
        write.table(data.frame(test=1), test_file)
        unlink(test_file)
      }, error = function(e) {
        stop("No write permission in directory: ", hyperparam_dir)
      })
    }
    
    if (use_slurm) {
      test_file <- file.path(run_topolow_dir, "test_write.txt")
      tryCatch({
        write.table(data.frame(test=1), test_file)
        unlink(test_file)
      }, error = function(e) {
        stop("No write permission in directory: ", run_topolow_dir)
      })
    }
  }

  # Generate LHS samples
  lhs_samples <- maximinLHS(n = num_samples, k = 4)
  lhs_params <- data.frame(
    N = floor(qunif(lhs_samples[,1], min = N_min, max = N_max + 1)),
    k0 = qunif(lhs_samples[,2], min = k0_min, max = k0_max),
    cqq = qunif(lhs_samples[,3], min = cqq_min, max = cqq_max),
    k_decay = qunif(lhs_samples[,4], min = k_decay_min, max = k_decay_max)
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

  # Save matrix_list if using SLURM
  if(use_slurm) {
    # Save matrix_list
    matrix_list_file <- file.path(run_topolow_dir,
                                  paste0(scenario_name, "_matrix_list.rds")
    )
    saveRDS(matrix_list, matrix_list_file)
    
    # Submit jobs to SLURM
    submit_parameter_jobs(
      matrix_list_file = matrix_list_file,
      lhs_params = lhs_params,
      max_iter = max_iter,
      relative_epsilon = relative_epsilon,
      convergence_counter = convergence_counter,
      scenario_name = scenario_name,
      folds = folds,
      num_cores = num_cores,
      output_dir = output_dir,
      cider = cider
    )
    return(invisible(NULL))

  } else {

    # Process samples locally with same file structure as SLURM
    if(verbose) cat("Processing samples locally\n")
    
    # Process each sample and fold
    process_sample <- function(i) {
      sample_idx <- ((i - 1) %% num_samples) + 1
      fold_idx <- floor((i - 1) / num_samples) + 1
      
      N <- lhs_params$N[sample_idx]
      k0 <- lhs_params$k0[sample_idx]
      cqq <- lhs_params$cqq[sample_idx]
      k_decay <- lhs_params$k_decay[sample_idx]
      
      truth_matrix <- matrix_list[[fold_idx]][[1]]
      input_matrix <- matrix_list[[fold_idx]][[2]]
      
      tryCatch({
        res_train <- topolow_full(
          distance_matrix = input_matrix,
          ndim = N,
          max_iter = max_iter,
          k0 = k0,
          k_decay = k_decay,
          cqq = cqq,
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
        mae_holdout <- mean(abs(df$OutSampleError), na.rm = TRUE)
        
        n <- sum(!is.na(df$OutSampleError))
        NLL <- n * (1 + log(2) + log(mae_holdout))
        
        # Return valid results
        if(is.finite(mae_holdout) && is.finite(NLL)) {
          result <- data.frame(
            N = N,
            k0 = k0,
            k_decay = k_decay,
            cqq = cqq,
            Holdout_MAE = mae_holdout,
            NLL = NLL
          )
          
          # Save individual result if requested
          if(write_files) {
            result_file <- file.path(run_topolow_dir,
                                     sprintf("%d_hyperparams_%s.csv", i, scenario_name))
            write.csv(result, result_file, row.names = FALSE)
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
    
    # Process all samples in parallel
    if(verbose) cat("Starting parallel processing\n")
    
    # Determine OS and set up appropriate parallel method
    if (.Platform$OS.type == "windows") {
      if(num_cores > 1) {
        if(verbose) cat("Using parallel cluster on Windows\n")
        # Create cluster
        cl <- parallel::makeCluster(num_cores)
        
        # Export required functions and data to cluster
        parallel::clusterExport(cl, c("matrix_list", "lhs_params", "max_iter",
                                      "relative_epsilon", "convergence_counter",
                                      "scenario_name", "write_files", "verbose",
                                      "run_topolow_dir", "num_samples"),
                                envir = environment())
        
        # Load required packages on each cluster node
        parallel::clusterEvalQ(cl, {
          library(topolow)
        })
        
        # Run processing
        res_list <- parallel::parLapply(cl, 1:(num_samples * folds), process_sample)
        
        # Clean up
        parallel::stopCluster(cl)
      } else {
        if(verbose) cat("Running sequentially on Windows with single core\n")
        res_list <- lapply(1:(num_samples * folds), process_sample)
      }
    } else {
      # Use mclapply on Unix-like systems
      if(verbose) cat("Using mclapply on Unix-like system\n")
      res_list <- parallel::mclapply(
        1:(num_samples * folds),
        process_sample,
        mc.cores = num_cores
      )
    }
    
    # Remove NULL results
    res_list <- Filter(Negate(is.null), res_list)
    
    # Check if we have any valid results
    if(length(res_list) == 0) {
      stop("No valid results obtained from parameter optimization")
    }
    
    # Combine results
    res_list <- do.call(rbind, res_list)
    colnames(res_list) <- c("N", "k0", "k_decay", "cqq", "Holdout_MAE", "NLL")
    
    # Remove any remaining invalid values
    res_list <- res_list[complete.cases(res_list) & 
                           apply(res_list, 1, function(x) all(is.finite(x))), ]
    
    if(nrow(res_list) == 0) {
      stop("All results were invalid after filtering infinities and NAs")
    }
    
    # Calculate median results across folds with error checking
    tryCatch({
      res_list_median <- aggregate(
        cbind(Holdout_MAE, NLL) ~ N + k0 + k_decay + cqq,
        data = res_list,
        FUN = median
      )
    }, error = function(e) {
      stop("Failed to aggregate results: ", e$message, 
           "\nNumber of valid results: ", nrow(res_list))
    })
    
    # Write aggregated results
    if(write_files) {
      file_name <- file.path(hyperparam_dir,
                             paste0(scenario_name, "_hyperparam_results.csv"))
      
      if(file.exists(file_name)) {
        existing_data <- read.csv(file_name)
        colnames(existing_data) <- colnames(res_list_median)  
        updated_data <- rbind(existing_data, res_list_median)
        write.csv(updated_data, file_name, row.names = FALSE)
      } else {
        write.csv(res_list_median, file_name, row.names = FALSE)
      }
    }
    
    return(res_list_median)
  }
}



#' Submit Parameter Optimization Jobs to SLURM
#'
#' @description
#' Helper function to submit parameter optimization jobs to SLURM cluster.
#'
#' @param matrix_list_file Path to saved matrix list RDS file
#' @param lhs_params Data frame of LHS parameter samples
#' @param max_iter Maximum iterations
#' @param relative_epsilon Convergence threshold
#' @param convergence_counter Convergence counter
#' @param scenario_name Scenario name
#' @param folds Number of CV folds
#' @param num_cores Cores per job
#' @param output_dir Directory for output files. If NULL, uses current directory
#' @param cider Whether to use cider queue
#' @return Invisible NULL
#' @keywords internal
submit_parameter_jobs <- function(matrix_list_file,
                                lhs_params,
                                max_iter,
                                relative_epsilon, 
                                convergence_counter,
                                scenario_name,
                                folds,
                                num_cores,
                                output_dir = NULL,
                                cider=FALSE) {
  
  # Handle output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  # Create directories
  slurm_dir <- file.path(output_dir, "run_topolow_arguments")
  if (!dir.exists(slurm_dir)) {
    dir.create(slurm_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Get path to run script
  script_path <- system.file(
    "scripts/run_topolow_arguments.R",
    package = "topolow"
  )
  
  if (!file.exists(script_path)) {
    stop("Could not find run script in package: run_topolow_arguments.R")
  }
  
  # Verify matrix_list before submitting jobs
  matrix_list <- readRDS(matrix_list_file)
  if (!is.list(matrix_list) || length(matrix_list) != folds) {
    stop(sprintf("Invalid matrix_list: expected length %d, got %d", 
                 folds, length(matrix_list)))
  }
  
  num_samples <- nrow(lhs_params)
  total_jobs <- num_samples * folds
  
  # Submit jobs - one per parameter set/fold combination
  for(i in seq_len(total_jobs)) {
    # Calculate indices
    sample_idx <- ((i - 1) %% num_samples) + 1
    
    # Prepare arguments
    args <- c(
      matrix_list_file,
      as.character(lhs_params$N[sample_idx]),
      as.character(max_iter),
      as.character(lhs_params$k0[sample_idx]),
      as.character(lhs_params$k_decay[sample_idx]),
      as.character(lhs_params$cqq[sample_idx]),
      as.character(relative_epsilon),
      as.character(convergence_counter),
      "NULL", # initial_positions
      "FALSE", # write_positions
      "FALSE", # verbose
      scenario_name,
      as.character(i),
      output_dir,
      as.character(num_samples)
    )
    
    # Create job script
    job_name <- paste0(i, "_run_topolow_arguments_", scenario_name)
    output_file <- file.path("run_topolow_arguments",
                            paste0(i, "_", scenario_name, ".out"))
    error_file <- file.path("run_topolow_arguments",
                           paste0(i, "_", scenario_name, ".err"))
    
    script_file <- create_slurm_script(
      job_name = job_name,
      script_path = script_path,
      args = args,
      num_cores = num_cores,
      output_file = output_file,
      error_file = error_file
    )
    
    # Submit job
    submit_job(script_file, use_slurm = TRUE, cider = cider)
  }
  
  return(invisible(NULL))
}


#' Aggregate Results from Parameter Optimization Jobs
#'
#' @description
#' Combines results from multiple parameter optimization jobs executed via SLURM
#' into a single dataset. This function processes results from jobs submitted by
#' \code{\link{submit_parameter_jobs}}.
#'
#' @details
#' The function looks for CSV files in the run_topolow_arguments directory that match 
#' the pattern "hyperparams_{scenario_name}.csv". It combines all results into a 
#' single dataset, computes median values across folds, and optionally writes the 
#' aggregated results to a file.
#'
#' The output file is saved as:
#' "hyperparam_results/{scenario_name}_hyperparam_results.csv"
#'
#' @param scenario_name Character. Name used in parameter optimization jobs.
#' @param write_files Logical. Whether to save combined results (default: TRUE).
#' @return Data frame of aggregated results containing median values across folds:
#'   \item{N}{Number of dimensions}
#'   \item{k0}{Initial spring constant}
#'   \item{k_decay}{Spring decay rate}
#'   \item{cqq}{Repulsion constant}
#'   \item{Holdout_MAE}{Median holdout mean absolute error}
#'   \item{NLL}{Median negative log likelihood}
#'
#' @examples
#' \dontrun{
#' # After running parameter optimization jobs:
#' results <- aggregate_parameter_optimization_results("optimization_run1")
#' }
#'
#' @seealso
#' \code{\link{run_parameter_optimization}} for running the optimization
#' \code{\link{submit_parameter_jobs}} for job submission
#'
#' @export
aggregate_parameter_optimization_results <- function(scenario_name, write_files = TRUE,
                                                     output_dir = NULL) {
  # Handle output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  # Find result files
  pattern <- paste0("_hyperparams_", scenario_name, "\\.csv$")
  directory_path <- file.path(output_dir, "run_topolow_arguments")
  
  csv_files <- list.files(
    path = directory_path,
    pattern = pattern,
    full.names = TRUE
  )
  
  if(length(csv_files) == 0) {
    stop("No result files found for scenario: ", scenario_name)
  }
  
  # Read and combine results
  results <- do.call(rbind, lapply(csv_files, function(file) {
    tryCatch({
      df <- read.csv(file)
      # Ensure numeric columns
      df[] <- lapply(df, function(x) {
        if(is.character(x)) {
          as.numeric(x)
        } else {
          x
        }
      })
      df
    }, error = function(e) {
      warning("Error reading file ", file, ": ", e$message)
      NULL
    })
  }))
  
  # Remove any rows with NAs
  results <- results[complete.cases(results), ]
  
  # Calculate median results across folds
  results_median <- aggregate(
    cbind(Holdout_MAE, NLL) ~ N + k0 + k_decay + cqq,
    data = results,
    FUN = median
  )
  
  if(write_files) {
    output_file <- file.path(
      "hyperparam_results",
      paste0(scenario_name, "_hyperparam_results.csv")
    )
    
    # Create directory if it doesn't exist
    dir.create(
      "hyperparam_results", 
      showWarnings = FALSE, 
      recursive = TRUE
    )
    
    # Save results
    write.csv(results_median, output_file, row.names = FALSE)

    # Delete the files
    files_to_delete <- list.files(path = directory_path, pattern = pattern, full.names = TRUE)
    file.remove(files_to_delete)

  }
  
  return(invisible(NULL))
}


#' Submit Adaptive Monte Carlo Sampling Jobs
#'
#' @description 
#' Performs adaptive Monte Carlo sampling to explore parameter space, either locally
#' or distributed via SLURM. Samples are drawn adaptively based on previous evaluations
#' to focus sampling in high-likelihood regions. Results from all jobs accumulate in
#' a single output file.
#'
#' @details
#' The function:
#' 1. Takes initial parameter samples as starting points
#' 2. Creates n_iter batches of batch_size samples each
#' 3. Updates sampling distribution based on likelihoods
#' 4. Can distribute computation via SLURM for large-scale sampling
#'
#' Both local and SLURM executions append results to the same output file:
#' "hyperparam_results/{scenario_name}_hyperparam_results.csv"
#'
#' @param initial_samples_file Character. Path to CSV file containing initial samples.
#'        Must contain columns: log_N, log_k0, log_k_decay, log_cqq, NLL
#' @param distance_matrix Matrix. Distance matrix to optimize.
#' @param max_iter Integer. Maximum iterations per sample evaluation.
#' @param relative_epsilon Numeric. Convergence threshold.
#' @param folds Integer. Number of CV folds (default: 10).
#' @param num_cores Integer. Cores per job (default: 1).
#' @param scenario_name Character. Name for output files.
#' @param num_samples Integer. Total number of jobs to submit.
#' @param n_iter Integer. Number of sampling iterations per job.
#' @param batch_size Integer. Samples per iteration (default: 1).
#' @param use_slurm Logical. Whether to use SLURM (default: FALSE).
#' @param cider Logical. Whether to use cider queue (default: FALSE).
#'
#' @return Invisible NULL. Results are appended to:
#'         "hyperparam_results/{scenario_name}_hyperparam_results.csv"
#'
#' @examples
#' \dontrun{
#' # Read initial samples
#' init_file <- "initial_samples.csv"
#' 
#' # Create distance matrix
#' dist_mat <- matrix(runif(100), 10, 10)
#' dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
#' diag(dist_mat) <- 0
#'
#' # Run local sampling
#' run_adaptive_sampling(
#'   initial_samples_file = init_file,
#'   distance_matrix = dist_mat,
#'   max_iter = 1000,
#'   scenario_name = "test_sampling",
#'   num_samples = 10,
#'   n_iter = 5
#' )
#'
#' # Run with SLURM
#' run_adaptive_sampling(
#'   initial_samples_file = init_file, 
#'   distance_matrix = dist_mat,
#'   scenario_name = "slurm_sampling",
#'   num_samples = 50,
#'   use_slurm = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{adaptive_MC_sampling}} for the core sampling algorithm
#'
#' @export
run_adaptive_sampling <- function(initial_samples_file,
                                  distance_matrix,
                                  num_samples = 5,
                                  n_iter = 1, 
                                  batch_size = 1,
                                  max_iter, 
                                  relative_epsilon = 1e-4,
                                  folds = 20, 
                                  num_cores = 1,
                                  scenario_name,
                                  output_dir = NULL,
                                  use_slurm = FALSE,
                                  cider = FALSE,
                                  verbose = FALSE) {
  
  if(!is.matrix(distance_matrix)) {
    stop("distance_matrix must be a matrix")
  }
  
  # Validate numeric parameters
  integer_params <- list(
    max_iter = max_iter,
    folds = folds, 
    num_cores = num_cores,
    num_samples = num_samples,
    n_iter = n_iter,
    batch_size = batch_size
  )
  
  for(param_name in names(integer_params)) {
    param_value <- integer_params[[param_name]]
    if(!is.numeric(param_value) || param_value != round(param_value) || param_value < 1) {
      stop(sprintf("%s must be a positive integer", param_name))
    }
  }
  
  if(!is.numeric(relative_epsilon) || relative_epsilon <= 0) {
    stop("relative_epsilon must be a positive number")
  }
  
  if(!is.character(scenario_name) || length(scenario_name) != 1) {
    stop("scenario_name must be a single character string")
  }
  
  # Check SLURM availability if requested
  if(use_slurm && !has_slurm()) {
    stop("SLURM requested but not available")
  }
  
  # Handle output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  # Create required directories
  adaptive_sampling_dir <- file.path(output_dir, "adaptive_sampling_jobs")
  hyperparam_dir <- file.path(output_dir, "hyperparam_results")
  
  for (dir in c(adaptive_sampling_dir, hyperparam_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Verify input file exists and has valid data
  if(!file.exists(initial_samples_file)) {
    stop("initial_samples_file not found: ", initial_samples_file)
  }
  
  # Read and validate initial samples
  init_samples <- read.csv(initial_samples_file)
  required_cols <- c("log_N", "log_k0", "log_k_decay", "log_cqq", "NLL")
  if(!all(required_cols %in% names(init_samples))) {
    stop("Missing required columns: ", 
         paste(setdiff(required_cols, names(init_samples)), collapse=", "))
  }
  
  # Verify we have valid initial samples
  init_samples <- init_samples[complete.cases(init_samples) & 
                                 apply(init_samples[required_cols], 1, 
                                       function(x) all(is.finite(x))), ]
  if(nrow(init_samples) == 0) {
    stop("No valid samples in initial samples file after filtering NAs and infinities")
  }
  
  if(use_slurm) {
    # SLURM handling remains the same
    distance_matrix_file <- file.path(adaptive_sampling_dir,
                                      paste0(scenario_name, "_distance_matrix.rds"))
    saveRDS(distance_matrix, distance_matrix_file)
    
    for(i in 1:num_samples) {
      # Prepare arguments
      args <- c(
        initial_samples_file,
        distance_matrix_file,
        as.character(max_iter),
        as.character(relative_epsilon),
        as.character(num_cores),
        scenario_name,
        as.character(i),
        as.character(n_iter),
        as.character(batch_size),
        output_dir,
        as.character(folds)
      )
      
      job_name <- paste0(i, "_adaptive_sampling_", scenario_name)
      output_file <- file.path(adaptive_sampling_dir,
                               paste0(i, "_", scenario_name, ".out"))
      error_file <- file.path(adaptive_sampling_dir,
                              paste0(i, "_", scenario_name, ".err"))
      
      script_file <- create_slurm_script(
        job_name = job_name,
        script_path = system.file("scripts/run_adaptive_sampling_jobs.R",
                                  package = "topolow"),
        args = args,
        num_cores = num_cores,
        output_file = output_file,
        error_file = error_file
      )
      
      submit_job(script_file, use_slurm = TRUE, cider = cider)
      
      if(verbose) cat(sprintf("Submitted sample %d with %d iterations\n", i, n_iter))
    }
    
    return(invisible(NULL))
    
  } else {
    # Run sampling locally, all num_samples will be processed by AMC as one batch
    tryCatch({
      results <- adaptive_MC_sampling(
        samples_file = initial_samples_file,
        distance_matrix = distance_matrix,
        n_iter = n_iter,
        batch_size = num_samples,
        max_iter = max_iter,
        relative_epsilon = relative_epsilon,
        folds = folds,
        num_cores = num_cores,
        scenario_name = scenario_name,
        output_dir = output_dir,
        verbose = verbose
      )
      
      if(is.null(results) || nrow(results) == 0) {
        warning("No valid results obtained from adaptive sampling")
        return(NULL)
      }
      
      return(results)
      
    }, error = function(e) {
      warning("Error in adaptive sampling: ", e$message)
      return(NULL)
    })
  }
}



#' Check Status of Adaptive Sampling Jobs
#'
#' @description
#' Monitors the status of submitted adaptive sampling jobs by checking
#' output files and SLURM queue.
#'
#' @param scenario_name Character. Name of sampling scenario.
#' @param expected_jobs Integer. Number of expected jobs.
#' @return Named list with components:
#'   \item{completed}{Number of completed jobs}
#'   \item{running}{Number of running jobs}
#'   \item{failed}{Number of failed jobs}
#'   \item{missing}{Number of missing jobs}
#'
#' @keywords internal
check_sampling_jobs <- function(scenario_name, expected_jobs) {
  # Initialize counters
  status <- list(
    completed = 0,
    running = 0,
    failed = 0,
    missing = 0
  )
  
  # Check for result files
  pattern <- paste0("_hyperparam_results_", scenario_name, "\\.csv$")
  result_files <- list.files(
    path = "hyperparam_results",
    pattern = pattern
  )
  status$completed <- length(result_files)
  
  # If SLURM available, check queue
  if(has_slurm()) {
    # Get running jobs
    running_jobs <- system(
      sprintf("squeue -n '*%s*' -h | wc -l", scenario_name),
      intern = TRUE
    )
    status$running <- as.numeric(running_jobs)
    
    # Check for failed jobs in error files
    error_files <- list.files(
      path = "adaptive_sampling_jobs",
      pattern = paste0("_", scenario_name, "\\.err$")
    )
    status$failed <- sum(file.size(error_files) > 0)
  }
  
  # Calculate missing jobs
  status$missing <- expected_jobs - sum(status$completed, status$running)
  
  return(status)
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
  
  names(samples) <- c("log_N", "log_k0", "log_k_decay", "log_cqq")
  return(samples)
}


#' Calculate Sampling Distribution Parameters
#'
#' @description
#' Calculates mean vector and covariance matrix for next sampling distribution
#' by weighting previous samples according to their likelihoods.
#'
#' @param samples Data frame of previous samples with parameter columns and NLL
#' @return List with:
#'   \item{mean}{Vector of weighted means for each parameter}
#'   \item{cov}{Weighted covariance matrix}
#' @keywords internal
get_sampling_prob_legacy <- function(samples) {
  # First, Remove outliers
  samples <- as.data.frame(lapply(samples, clean_data, k = 3))
  samples <- na.omit(samples)
  
  # Calculate weights based on likelihoods
  log_likelihoods <- -samples$NLL  
  
  std_log_likelihoods <- log_likelihoods - min(log_likelihoods) + 0.05
  weights <- std_log_likelihoods / sum(std_log_likelihoods)

  # Calculate weighted mean
  weighted_mean_vector <- sapply(samples[, 1:4], weighted_mean, weights = weights)

  # Calculate weighted covariance
  weighted_cov <- cov.wt(samples[, 1:4], wt = weights)$cov
  
  return(list(mean = weighted_mean_vector, cov = weighted_cov))
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
  par_names <- c("log_N", "log_k0", "log_k_decay", "log_cqq")
  
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


#' Calculate Weighted Mean
#'
#' @description Helper function to calculate weighted mean of a vector.
#' @param x Numeric vector
#' @param weights Vector of weights
#' @return Weighted mean value
#' @keywords internal
weighted_mean <- function(x, weights) {
  sum(x * weights) / sum(weights)
}


#' Weighted Kernel Density Estimation
#'
#' @description
#' Performs weighted kernel density estimation for univariate data. Useful for
#' analyzing parameter distributions with importance weights.
#'
#' @param x Numeric vector of samples
#' @param weights Numeric vector of weights
#' @param n Integer number of evaluation points
#' @param from,to Numeric range for evaluation points
#' @return List containing:
#'   \item{x}{Vector of evaluation points}
#'   \item{y}{Vector of density estimates}
#' @export
weighted_kde <- function(x, weights, n = 512, from = min(x), to = max(x)) {
  # Normalize weights
  weights <- weights / sum(weights)
  
  # Calculate bandwidth (Silverman's rule)
  bw <- 1.06 * sd(x) * length(x)^(-1/5)
  
  # Generate evaluation points
  eval_points <- seq(from, to, length.out = n)
  
  # Compute density
  compute_density <- function(z) {
    sum(weights * dnorm(z, mean = x, sd = bw))
  }
  # Limit cores to prevent excessive process spawning
  # Use minimum of: available cores - 1, 2, or what system allows
  num_cores <- min(parallel::detectCores() - 1, 2)
  
  if (.Platform$OS.type == "windows") {
    # Run sequentially on Windows
    density_est <- sapply(eval_points, compute_density)
  } else {
    # Run in parallel on Unix-like systems
    density_est <- parallel::mclapply(eval_points, compute_density, 
                                    mc.cores = num_cores)
    density_est <- unlist(density_est)
  }
  
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
#' @export
unweighted_kde <- function(x, n = 512, from = min(x), to = max(x), 
                          bw = "nrd0") {
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
#'
#' @param distance_matrix Distance matrix to fit
#' @param max_iter Maximum optimization iterations
#' @param relative_epsilon Convergence threshold
#' @param N Number of dimensions
#' @param k0 Initial spring constant
#' @param k_decay Spring constant decay rate
#' @param cqq Repulsion constant
#' @param folds Number of CV folds
#' @param num_cores Number of cores for parallel processing
#' @return List with:
#'   \item{Holdout_MAE}{Mean absolute error on validation data}
#'   \item{NLL}{Negative log likelihood}
#' @keywords internal
likelihood_function <- function(distance_matrix, max_iter,
                              relative_epsilon, N, k0, k_decay,
                              cqq, folds = 20, num_cores) {
  # Create n folds in the data
  truth_matrix <- distance_matrix
  
  matrix_list <- replicate(folds, 
                          list(truth_matrix, NULL),
                          simplify = FALSE)
  
  num_elements <- sum(!is.na(distance_matrix))
  
  holdout_size <- floor(num_elements/(folds*2)) # 2 is because for each [i,j] we also null the [j,i] element
  # To cover n folds randomly,create a copy to remove fractions from it until finished:
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
  
  # Process each fold in parallel
  process_sample <- function(i) {
    truth_matrix <- matrix_list[[i]][[1]]
    input_matrix <- matrix_list[[i]][[2]]
    
    res_train <- topolow_full(distance_matrix = input_matrix, 
                             ndim = N, 
                             max_iter = max_iter, 
                             k0 = k0, 
                             k_decay = k_decay, 
                             cqq = cqq,
                             relative_epsilon = relative_epsilon,
                             convergence_counter = 5,
                             initial_positions = NULL,
                             write_positions_to_csv = FALSE,
                             verbose = FALSE)
    
    p_dist_mat <- res_train$est_distances
    p_dist_mat <- as.matrix(p_dist_mat)
    
    errors <- error_calculator_comparison(p_dist_mat, truth_matrix, input_matrix)
    df <- errors$report_df
    mae_holdout <- mean(abs(df$OutSampleError), na.rm = TRUE)
    
    data.frame(N = N, k0 = k0, k_decay = k_decay, cqq = cqq, 
               Holdout_MAE = mae_holdout)
  }
  
  # Run in parallel
  res_list <- mclapply(1:folds, process_sample, mc.cores = num_cores)
  
  # Combine results
  res_list <- do.call(rbind, res_list)
  colnames(res_list) <- c("N", "k0", "k_decay", "cqq", "Holdout_MAE")
  
  res_list_median <- aggregate(Holdout_MAE ~ N + k0 + k_decay + cqq, 
                             data = res_list, FUN = median)
  
  n <- holdout_size
  Holdout_MAE <- res_list_median$Holdout_MAE
  NLL <- n*(1 + log(2) + log(Holdout_MAE))
  
  return(list(Holdout_MAE=Holdout_MAE, NLL=NLL))
}


#' Safe Wrapper for Likelihood Evaluation
#'
#' @description
#' Wraps likelihood calculation in error handler to return NA on failure.
#' @param ... Arguments passed to likelihood_function
#' @return Same as likelihood_function or NA if error
#' @keywords internal
safe_likelihood_function <- function(...) {
  tryCatch({
    likelihood_function(...)
  }, error = function(e) {
    cat("Error in safe_likelihood_function:", e$message, "\n")
    # If an error occurs, return NA or a placeholder value
    NA  # You can return NA or a small likelihood if preferred, e.g. 0.01
  })
}


#' Perform Adaptive Monte Carlo Sampling
#'
#' @description
#' Main function implementing adaptive Monte Carlo sampling (https://www.sciencedirect.com/science/article/pii/0167473088900203) 
#' to explore parameter space. Updates sampling distribution based on evaluated likelihoods.
#'
#' @param samples_file Path to CSV with initial samples
#' @param distance_matrix Distance matrix to fit
#' @param n_iter Number of sampling iterations
#' @param batch_size Samples per iteration
#' @param max_iter Maximum optimization iterations 
#' @param relative_epsilon Convergence threshold
#' @param folds Number of CV folds
#' @param num_cores Number of cores for parallel processing
#' @param scenario_name Name for output files
#' @param replace_csv Whether to replace existing CSV
#' @return Data frame of samples with evaluated likelihoods
#' @export
adaptive_MC_sampling_legacy <- function(samples_file, distance_matrix,
                               n_iter = 1, batch_size = 1,
                               max_iter, relative_epsilon,
                               folds = 20, num_cores,
                               scenario_name, replace_csv) {

  par_names <- c("log_N", "log_k0", "log_k_decay", "log_cqq")
  
  for (iter in 1:n_iter) {
    # Get updated samples
    current_samples <- read.csv(file.path("~/scripts/Model1", samples_file), 
                              stringsAsFactors = FALSE)
    current_samples <- current_samples[apply(current_samples, 1, 
                                           function(row) all(is.finite(row))), ]
    current_samples <- na.omit(current_samples)
    
    # Check convergence
    if(nrow(current_samples) > 1500){
      conv_check <- check_gaussian_convergence(data=current_samples[,par_names], 
                                             window_size = 500, 
                                             tolerance = 0.002)
      if (conv_check$converged) {
        cat("Convergence achieved at iteration", iter, "\n")
        break
      }
    }

    # Get adjusted proposal distribution
    dist <- get_sampling_prob(samples= current_samples[, c(par_names, "NLL")])
    
    # Sample from updated proposal
    new_samples <- sample_from_distribution(dist, batch_size, epsilon = 0.05)
    
    new_samples$log_N <- as.numeric(new_samples$log_N)
    new_samples$log_k0 <- as.numeric(new_samples$log_k0)
    new_samples$log_k_decay <- as.numeric(new_samples$log_k_decay)
    new_samples$log_cqq <- as.numeric(new_samples$log_cqq)

    # Transform parameters
    new_samples$N <- round(exp(new_samples$log_N))
    new_samples$k0 <- exp(new_samples$log_k0)
    new_samples$k_decay = exp(new_samples$log_k_decay)
    new_samples$cqq = exp(new_samples$log_cqq)

    # Evaluate likelihood
    new_likelihoods <- sapply(1:nrow(new_samples), function(i) 
      safe_likelihood_function(
        distance_matrix=distance_matrix,
        max_iter=max_iter,
        relative_epsilon=relative_epsilon, 
        N= new_samples[i, "N"],
        k0= new_samples[i, "k0"],
        k_decay= new_samples[i, "k_decay"],
        cqq= new_samples[i, "cqq"],
        folds = folds,
        num_cores=num_cores))
    
    if(!all(is.na(new_likelihoods))){
      # Process results
      new_likelihoods_df <- as.data.frame(t(new_likelihoods))
      colnames(new_likelihoods_df) <- c("Holdout_MAE", "NLL")
      
      new_samples <- cbind(new_samples, new_likelihoods_df)
      new_samples$Holdout_MAE <- as.numeric(new_samples$Holdout_MAE)
      new_samples$NLL <- as.numeric(new_samples$NLL)

      new_samples <- new_samples[, names(current_samples)]
      
      # Save results
      for (i in 1:nrow(new_samples)) {
        file_name <- paste0("hyperparam_results/", scenario_name, 
                           "_hyperparam_results.csv")
        
        if (file.exists(file_name)) {
          existing_data <- read.csv(file_name)
          colnames(existing_data) <- colnames(new_samples)
          updated_data <- rbind(existing_data, new_samples)
        } else {
          updated_data <- new_samples
        }
        write.csv(updated_data, 
                 file.path("~/scripts/Model1", file_name), 
                 row.names = FALSE)
      }
    }
  }
  
  return(new_samples)
}


#' Perform Adaptive Monte Carlo Sampling
#'
#' @description
#' Main function implementing adaptive Monte Carlo sampling (https://www.sciencedirect.com/science/article/pii/0167473088900203) 
#' to explore parameter space. Updates sampling distribution based on evaluated likelihoods.
#'
#' @param samples_file Path to CSV with initial samples
#' @param distance_matrix Distance matrix to fit
#' @param n_iter Number of sampling iterations
#' @param batch_size Samples per iteration
#' @param max_iter Maximum optimization iterations 
#' @param relative_epsilon Convergence threshold
#' @param folds Number of CV folds
#' @param num_cores Number of cores for parallel processing
#' @param scenario_name Name for output files
#' @param output_dir Directory for output files. If NULL, uses current directory
#' @param replace_csv Whether to replace existing CSV
#' @return Data frame of samples with evaluated likelihoods
#' @export
adaptive_MC_sampling <- function(samples_file, 
                                 distance_matrix,
                                 n_iter = 1, 
                                 batch_size = 1,
                                 max_iter, 
                                 relative_epsilon,
                                 folds = 20, 
                                 num_cores = 1,
                                 scenario_name, 
                                 output_dir = NULL,
                                 verbose = FALSE) {
  
  # Handle parallel processing setup
  if(.Platform$OS.type == "windows" && num_cores > 1) {
    if(verbose) cat("Setting up parallel cluster for Windows\n")
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl))
    
    parallel::clusterExport(cl, c("distance_matrix", "max_iter", 
                                  "relative_epsilon", "folds"),
                            envir = environment())
    
    parallel::clusterEvalQ(cl, {
      library(topolow)
    })
  }
  
  # Handle output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  # Create results directory
  hyperparam_dir <- file.path(output_dir, "hyperparam_results")
  if (!dir.exists(hyperparam_dir)) {
    dir.create(hyperparam_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  par_names <- c("log_N", "log_k0", "log_k_decay", "log_cqq")
  
  for (iter in 1:n_iter) {
    if(verbose) cat(sprintf("\nStarting iteration %d of %d\n", iter, n_iter))
    
    # Read current samples
    current_samples <- read.csv(samples_file)
    current_samples <- current_samples[apply(current_samples, 1, 
                                             function(row) all(is.finite(row))), ]
    current_samples <- na.omit(current_samples)
    
    if(nrow(current_samples) == 0) {
      warning("No valid samples remaining after filtering")
      break
    }
    
    # Check convergence
    if(nrow(current_samples) > 1000) {
      conv_check <- check_gaussian_convergence(data=current_samples[,par_names], 
                                               window_size = 500, 
                                               tolerance = 0.002)
      if (conv_check$converged) {
        if(verbose) cat("Convergence achieved at iteration", iter, "\n")
        break
      }
    }
    
    # Generate new samples using KDE
    new_samples <- generate_kde_samples(
      samples = current_samples,
      n = batch_size
    )
    
    # Ensure numeric columns
    for (col in par_names) {
      new_samples[[col]] <- as.numeric(new_samples[[col]])
    }
    
    # Evaluate samples with appropriate parallel method
    evaluate_sample <- function(i) {
      safe_likelihood_function(
        distance_matrix = distance_matrix,
        max_iter = max_iter,
        relative_epsilon = relative_epsilon, 
        N = round(exp(new_samples[i, "log_N"])),
        k0 = exp(new_samples[i, "log_k0"]),
        k_decay = exp(new_samples[i, "log_k_decay"]),
        cqq = exp(new_samples[i, "log_cqq"]),
        folds = folds,
        num_cores = 1  # Use single core within parallel processes
      )
    }
    
    if(.Platform$OS.type == "windows" && num_cores > 1) {
      new_likelihoods <- parallel::parLapply(cl, 1:nrow(new_samples), 
                                             evaluate_sample)
    } else if(.Platform$OS.type != "windows" && num_cores > 1) {
      new_likelihoods <- parallel::mclapply(1:nrow(new_samples), 
                                            evaluate_sample,
                                            mc.cores = num_cores)
    } else {
      new_likelihoods <- lapply(1:nrow(new_samples), evaluate_sample)
    }
    
    # Convert list to matrix
    new_likelihoods <- do.call(rbind, new_likelihoods)
    
    if (!all(is.na(new_likelihoods))) {
      # Process results
      new_likelihoods_df <- as.data.frame(new_likelihoods)
      colnames(new_likelihoods_df) <- c("Holdout_MAE", "NLL")
      
      new_samples <- cbind(new_samples, new_likelihoods_df)
      new_samples$Holdout_MAE <- as.numeric(new_samples$Holdout_MAE)
      new_samples$NLL <- as.numeric(new_samples$NLL)
      
      # Match columns with current samples
      new_samples <- new_samples[, names(current_samples)]
      
      # Remove invalid results
      new_samples <- new_samples[complete.cases(new_samples) & 
                                   apply(new_samples, 1, function(x) all(is.finite(x))), ]
      
      if (nrow(new_samples) > 0) {
        # Save results
        result_file <- file.path(hyperparam_dir,
                                 paste0(scenario_name, "_hyperparam_results.csv"))
        
        if (file.exists(result_file)) {
          existing_data <- read.csv(result_file)
          colnames(existing_data) <- colnames(new_samples)
          updated_data <- rbind(existing_data, new_samples)
          write.csv(updated_data, result_file, row.names = FALSE)
        } else {
          write.csv(new_samples, result_file, row.names = FALSE)
        }
        
        if(verbose) {
          cat(sprintf("Added %d new valid samples\n", nrow(new_samples)))
        }
      } else {
        if(verbose) cat("No valid samples in this iteration\n")
      }
    } else {
      if(verbose) cat("All likelihood evaluations failed in this iteration\n")
    }
  }
  
  # Return final samples if any
  result_file <- file.path(hyperparam_dir,
                           paste0(scenario_name, "_hyperparam_results.csv"))
  if(file.exists(result_file)) {
    final_samples <- read.csv(result_file)
    return(final_samples)
  } else {
    warning("No results file created")
    return(NULL)
  }
}


#' Calculate Weighted Marginal Distributions in Parallel
#'
#' @description
#' Calculates marginal distributions for each parameter with weights derived from 
#' log-likelihoods. Uses parallel processing for efficiency.
#'
#' @param samples Data frame containing:
#'        - log_N, log_k0, log_k_decay, log_cqq: Parameter columns
#'        - NLL: Negative log-likelihood column
#' @return Named list of marginal distributions, each containing:
#'   \item{x}{Vector of parameter values}
#'   \item{y}{Vector of density estimates}
#' @details 
#' Uses kernel density estimation weighted by normalized likelihoods.
#' Parallelizes computation across parameter dimensions using mclapply.
#'
#' @importFrom parallel mclapply detectCores
#' @importFrom stats dnorm sd
#' @export
calculate_weighted_marginals <- function(samples) {
  # Input validation
  required_cols <- c("log_N", "log_k0", "log_k_decay", "log_cqq", "NLL")
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
  vars <- c("log_N", "log_k0", "log_k_decay", "log_cqq")
  
  # Define KDE calculation function
  calculate_kde <- function(var, weights) {
    weighted_kde(samples[[var]], weights = weights)
  }
  
  # Calculate marginals in parallel
  marginals <- mclapply(vars, calculate_kde, 
                       weights = weights, 
                       mc.cores = detectCores())
  
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
#' estimate conditional likelihoods.
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
#' \dontrun{
#' # Calculate profile likelihood for parameter "log_N"
#' pl <- profile_likelihood("log_N", mcmc_samples, 
#'                         grid_size = 60,
#'                         bandwidth_factor = 0.02)
#'                         
#' # Plot results
#' plot(pl)
#' }
#' @seealso 
#' \code{\link{plot.profile_likelihood}} for visualization
#' @export
profile_likelihood <- function(param, samples, grid_size = 48, 
                             bandwidth_factor = 0.03,
                             start_factor = 0.5, end_factor = 1.5,
                             min_samples = 10) {
  
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

#' Print Method for Profile Likelihood Objects
#'
#' @param x Profile likelihood object
#' @param ... Additional arguments passed to print
#' @method print profile_likelihood
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
}