# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/adaptive_sampling.R

#' Parameter Space Sampling and Optimization Functions for topolow
#'


#' Initial Parameter Optimization using Latin Hypercube Sampling
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
#' **Note**: As of version 2.0.0, this function returns log-transformed parameters directly,
#' eliminating the need to call `log_transform_parameters()` separately.
#'
#' @details
#' The function performs these steps:
#' 1. Generates LHS samples in the parameter space (original scale for sampling).
#' 2. If \code{opt_subsample} is specified, each parameter evaluation gets a different
#'    random subsample for robustness testing.
#' 3. For each parameter set, \code{likelihood_function} is called which internally
#'    creates k-fold splits and evaluates via cross-validation.
#' 4. Parameter evaluations are processed in batches for memory efficiency.
#' 5. Computations are run in parallel within each batch.
#' 6. Automatically log-transforms the final results for direct use with adaptive sampling.
#'
#' **Note on Cross-Validation with Subsampling:**
#' When \code{opt_subsample} is used, each parameter evaluation receives a different
#' random subsample. Since the underlying data differs between evaluations, each evaluation
#' also gets its own CV fold structure (created internally by \code{likelihood_function}).
#' This approach tests parameter robustness across different data samples and fold structures,
#' which is appropriate for the exploration phase of parameter optimization.
#'
#' @param dissimilarity_matrix Matrix. Input dissimilarity matrix. Must be square and symmetric.
#' @param mapping_max_iter Integer. Maximum number of optimization iterations for each map.
#' @param relative_epsilon Numeric. Convergence threshold for relative change in error.
#' @param convergence_counter Integer. Number of iterations below threshold before declaring convergence.
#' @param scenario_name Character. Name for output files and job identification.
#' @param N_min,N_max Integer. Range for the number of dimensions parameter.
#' @param k0_min,k0_max Numeric. Range for the initial spring constant parameter.
#' @param c_repulsion_min,c_repulsion_max Numeric. Range for the repulsion constant parameter.
#' @param cooling_rate_min,cooling_rate_max Numeric. Range for the cooling rate parameter.
#' @param num_samples Integer. Number of LHS samples to generate. Default: 20.
#' @param max_cores Integer. Maximum number of cores for parallel processing. Default: NULL (uses all but one).
#' @param folds Integer. Number of cross-validation folds. Default: 20.
#' @param opt_subsample Integer or NULL. If specified, randomly samples this many
#'   points from the dissimilarity matrix for each parameter evaluation to reduce
#'   computational cost. The function automatically validates that subsampled data
#'   forms a connected graph. Each parameter evaluation uses a different random
#'   subsample for robustness. Default: NULL (use full data).
#'
#'   **Notes:**
#'   \itemize{
#'     \item If connectivity cannot be achieved, sample size is adaptively increased
#'     \item Minimum recommended: \code{opt_subsample >= max(100, folds)}
#'     \item Each parameter set gets a different subsample, but all CV folds within
#'       that parameter set use the same subsample
#'     \item The actual subsample size used is reported in output columns
#'   }
#' @param verbose Logical. Whether to print progress messages. Default: FALSE.
#' @param write_files Logical. Whether to save results to a CSV file. Default: FALSE.
#' @param output_dir Character. Directory for output files. Required if `write_files` is TRUE.
#'
#' @return A `data.frame` containing the log-transformed parameter sets and their performance metrics.
#'   Columns include: `log_N`, `log_k0`, `log_cooling_rate`, `log_c_repulsion`, `Holdout_MAE`, `NLL`,
#'   and if \code{opt_subsample} was used: `opt_subsample`,`original_n_points`.
#'
#' @note
#' \strong{Breaking Change in v2.0.0:} This function now returns log-transformed parameters directly.
#' The returned data frame has columns `log_N`, `log_k0`, `log_cooling_rate`, `log_c_repulsion`
#' instead of the original scale parameters. This eliminates the need to call `log_transform_parameters()`
#' separately before using `run_adaptive_sampling()`.
#'
#' \strong{Breaking Change in v2.0.0:} The parameter \code{distance_matrix} has been renamed to
#' \code{dissimilarity_matrix}. Please update your code accordingly.
#'
#' @examples
#' \donttest{
#' # This example can exceed 5 seconds on some systems.
#' # 1. Create a simple synthetic dataset for the example
#' synth_coords <- matrix(rnorm(60), nrow = 20, ncol = 3)
#' dist_mat <- coordinates_to_matrix(synth_coords)
#'
#' # 2. Run the optimization on the synthetic data (full data)
#' results <- initial_parameter_optimization(
#'   dissimilarity_matrix = dist_mat,
#'   mapping_max_iter = 100,
#'   relative_epsilon = 1e-3,
#'   convergence_counter = 2,
#'   scenario_name = "test_opt_synthetic",
#'   N_min = 2, N_max = 5,
#'   k0_min = 1, k0_max = 10,
#'   c_repulsion_min = 0.001, c_repulsion_max = 0.05,
#'   cooling_rate_min = 0.001, cooling_rate_max = 0.02,
#'   num_samples = 4,
#'   max_cores = 1,
#'   verbose = FALSE
#' )
#'
#' # 3. With subsampling for faster computation
#' results_sub <- initial_parameter_optimization(
#'   dissimilarity_matrix = dist_mat,
#'   mapping_max_iter = 100,
#'   relative_epsilon = 1e-3,
#'   convergence_counter = 2,
#'   scenario_name = "test_opt_subsampled",
#'   N_min = 2, N_max = 5,
#'   k0_min = 1, k0_max = 10,
#'   c_repulsion_min = 0.001, c_repulsion_max = 0.05,
#'   cooling_rate_min = 0.001, cooling_rate_max = 0.02,
#'   num_samples = 4,
#'   max_cores = 1,
#'   folds = 10,
#'   opt_subsample = 15,  # Use only 15 points
#'   verbose = TRUE
#' )
#' }
#'
#' @seealso \code{\link{euclidean_embedding}} for the core optimization algorithm,
#' \code{\link{subsample_dissimilarity_matrix}} for subsampling details.
#'
#' @importFrom lhs maximinLHS
#' @importFrom stats qunif complete.cases aggregate
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply
#' @importFrom future availableCores
#' @export
initial_parameter_optimization <- function(dissimilarity_matrix,
                                           mapping_max_iter = 500,
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
                                           num_samples = 20,
                                           max_cores = NULL,
                                           folds = 20,
                                           opt_subsample = NULL,
                                           verbose = FALSE,
                                           write_files = FALSE,
                                           output_dir) {
  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================
  
  # Matrix validation
  if (!is.matrix(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be a matrix")
  }
  if (nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be square")
  }
  matrix_size <- nrow(dissimilarity_matrix)
  
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
    if (!is.numeric(param_value) || length(param_value) != 1 ||
        is.na(param_value) || param_value <= 0) {
      stop(sprintf("%s must be a positive number", param_name))
    }
  }
  
  # Validate ranges
  if (N_min > N_max) stop("N_min must be less than N_max")
  if (k0_min >= k0_max) stop("k0_min must be less than k0_max")
  if (c_repulsion_min >= c_repulsion_max) stop("c_repulsion_min must be less than c_repulsion_max")
  if (cooling_rate_min >= cooling_rate_max) stop("cooling_rate_min must be less than cooling_rate_max")
  
  # Validate logical parameters
  logical_params <- list(verbose = verbose, write_files = write_files)
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
  
  # Validate opt_subsample
  use_subsampling <- FALSE
  if (!is.null(opt_subsample)) {
    if (!is.numeric(opt_subsample) || length(opt_subsample) != 1 ||
        opt_subsample < 5) {
      stop("opt_subsample must be NULL or a numeric value >= 5")
    }
    opt_subsample <- as.integer(floor(opt_subsample))
    
    if (opt_subsample >= matrix_size) {
      if (verbose) {
        cat(sprintf("opt_subsample (%d) >= matrix size (%d). Using full matrix.\n",
                    opt_subsample, matrix_size))
      }
      opt_subsample <- NULL
    } else {
      use_subsampling <- TRUE
      
      # Warn if subsample seems too small
      recommended_min <- max(50, folds)
      if (opt_subsample < recommended_min) {
        warning(sprintf(
          paste0("opt_subsample (%d) is smaller than recommended minimum (%d). ",
                 "This may lead to connectivity issues or unreliable results."),
          opt_subsample, recommended_min
        ), call. = FALSE)
      }
    }
  }
  

  # ==========================================================================
  # PRE-FLIGHT DATA QUALITY CHECK
  # ==========================================================================
  if (verbose) {
    cat("\nPerforming pre-flight data quality checks...\n")
    
    # Check sparsity
    n_missing <- sum(is.na(dissimilarity_matrix))
    sparsity <- n_missing / (matrix_size^2)
    cat(sprintf("  Sparsity: %.1f%% missing values\n", sparsity * 100))
    
    # Check connectivity if using subsampling
    if (use_subsampling && requireNamespace("igraph", quietly = TRUE)) {
      cat("  Checking full matrix connectivity...\n")
      full_connectivity <- tryCatch({
        check_matrix_connectivity(dissimilarity_matrix, min_completeness = 0.01)
      }, error = function(e) {
        list(is_connected = NA, completeness = 0)
      })
      
      if (!is.na(full_connectivity$is_connected)) {
        cat(sprintf("    Full matrix: %s (%d components, %.1f%% complete)\n",
                    ifelse(full_connectivity$is_connected, "CONNECTED", "DISCONNECTED"),
                    full_connectivity$n_components,
                    full_connectivity$completeness * 100))
        
        # Warn if very sparse
        if (sparsity > 0.95) {
          warning(
            sprintf(
              paste0("\n  DATA SPARSITY WARNING:\n",
                    "  Your matrix is %.1f%% sparse (missing values).\n",
                    "  Subsampling may create disconnected graphs.\n",
                    "  Consider pruning with prune_sparse_matrix() first.\n"),
              sparsity * 100
            ),
            call. = FALSE, immediate. = TRUE
          )
        }
        
        # Test subsample feasibility
        if (use_subsampling && sparsity > 0.90) {
          cat("  Testing subsample feasibility...\n")
          test_subsample <- tryCatch({
            subsample_dissimilarity_matrix(
              dissimilarity_matrix = dissimilarity_matrix,
              sample_size = opt_subsample,
              max_attempts = 2,
              verbose = FALSE
            )
          }, error = function(e) NULL)
          
          if (is.null(test_subsample) || !test_subsample$is_connected) {
            warning(
              sprintf(
                paste0("\n  SUBSAMPLE CONNECTIVITY WARNING:\n",
                      "  Test subsample (n=%d) failed connectivity check.\n",
                      "  All %d parameter evaluations may fail.\n",
                      "  STRONGLY RECOMMEND: Prune matrix before optimization.\n"),
                opt_subsample, num_samples
              ),
              call. = FALSE, immediate. = TRUE
            )
          } else {
            cat(sprintf("    Test subsample: CONNECTED (%.1f%% complete)\n",
                        test_subsample$completeness * 100))
          }
        }
      }
    }
    
    cat("Pre-flight checks complete.\n\n")
  }

  # ==========================================================================
  # PARALLEL PROCESSING SETUP
  # ==========================================================================
  available_cores <- future::availableCores()
  if (is.null(max_cores)) {
    max_cores <- max(1, available_cores - 1, na.rm = TRUE)
  } else {
    max_cores <- min(max_cores, available_cores)
  }
  
  if (verbose) {
    cat(sprintf("Processing %d parameter samples with up to %d cores\n",
                num_samples, max_cores))
    if (use_subsampling) {
      cat(sprintf("Using subsampling: %d points (from %d total)\n",
                  opt_subsample, matrix_size))
    }
  }
  
  # ==========================================================================
  # DIRECTORY SETUP (if writing files)
  # ==========================================================================
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
  }
  
  # ==========================================================================
  # LHS PARAMETER SAMPLING
  # ==========================================================================
  lhs_samples <- lhs::maximinLHS(n = num_samples, k = 4)
  lhs_params <- data.frame(
    N = floor(stats::qunif(lhs_samples[,1], min = N_min, max = N_max + 1)),
    k0 = stats::qunif(lhs_samples[,2], min = k0_min, max = k0_max),
    c_repulsion = stats::qunif(lhs_samples[,3], min = c_repulsion_min, max = c_repulsion_max),
    cooling_rate = stats::qunif(lhs_samples[,4], min = cooling_rate_min, max = cooling_rate_max)
  )
  
  # ==========================================================================
  # PARAMETER EVALUATION FUNCTION (PROCESSES ONE PARAMETER SET)
  # ==========================================================================
  # Store reference to full matrix for subsampling
  full_dissimilarity_matrix <- dissimilarity_matrix
  
  # Initialize failure tracking
  failure_diagnostics <- list(
    subsample_failures = 0,
    cv_fold_failures = 0,
    embedding_failures = 0,
    other_failures = 0,
    failure_messages = c()
  )

  process_param_set <- function(i) {
    tryCatch({
      # Get parameters for this evaluation
      N <- lhs_params$N[i]
      k0 <- lhs_params$k0[i]
      c_repulsion <- lhs_params$c_repulsion[i]
      cooling_rate <- lhs_params$cooling_rate[i]
      
      if (verbose) {
        cat(sprintf("\nEvaluating sample %d/%d: N=%d, k0=%.3f, c_rep=%.4f, cool=%.4f\n",
                    i, num_samples, N, k0, c_repulsion, cooling_rate))
      }
      
      # ======================================================================
      # SUBSAMPLING (if enabled)
      # ======================================================================
      local_matrix <- full_dissimilarity_matrix
      
      if (use_subsampling) {
        if (verbose) {
          cat(sprintf("  Subsampling to %d points...\n", opt_subsample))
        }
        
        subsample_result <- tryCatch({
          subsample_dissimilarity_matrix(
            dissimilarity_matrix = full_dissimilarity_matrix,
            sample_size = opt_subsample,
            max_attempts = 5,
            verbose = verbose
          )
        }, error = function(e) {
          # CHANGE: Track subsample failure with details
          failure_diagnostics$subsample_failures <<- failure_diagnostics$subsample_failures + 1
          failure_diagnostics$failure_messages <<- c(
            failure_diagnostics$failure_messages,
            sprintf("Sample %d: Subsample failed - %s", i, e$message)
          )
          if (verbose) {
            cat("  X Subsampling failed:", e$message, "\n")
          }
          return(NULL)
        })
        
        if (is.null(subsample_result)) {
          return(NULL)
        }

        if (!subsample_result$is_connected) {
          # Track connectivity failure with details
          failure_diagnostics$subsample_failures <<- failure_diagnostics$subsample_failures + 1
          failure_diagnostics$failure_messages <<- c(
            failure_diagnostics$failure_messages,
            sprintf("Sample %d: Disconnected subsample (%d components, %.1f%% complete)",
                    i, subsample_result$n_components, 
                    subsample_result$completeness * 100)
          )
          if (verbose) {
            cat("  X Failed to obtain connected subsample. Skipping this parameter set.\n")
          }
          return(NULL)
        }
        
        local_matrix <- subsample_result$subsampled_matrix
        
        # Sanity check the subsampled data
        sanity_result <- sanity_check_subsample(
          local_matrix,
          folds = folds,
          verbose = verbose
        )
        
        if (!sanity_result$all_checks_passed && verbose) {
          cat("  [!] Subsample sanity checks raised warnings (proceeding anyway)\n")
        }
      }
      
      # ======================================================================
      # CROSS-VALIDATION EVALUATION
      # ======================================================================
      result <- tryCatch({
        likelihood_function(
          dissimilarity_matrix = local_matrix,
          mapping_max_iter = mapping_max_iter,
          relative_epsilon = relative_epsilon,
          N = N,
          k0 = k0,
          cooling_rate = cooling_rate,
          c_repulsion = c_repulsion,
          folds = folds,
          num_cores = 1
        )
      }, error = function(e) {
        # Track CV/embedding failure with details
        if (grepl("fold|sparse", e$message, ignore.case = TRUE)) {
          failure_diagnostics$cv_fold_failures <<- failure_diagnostics$cv_fold_failures + 1
        } else {
          failure_diagnostics$embedding_failures <<- failure_diagnostics$embedding_failures + 1
        }
        failure_diagnostics$failure_messages <<- c(
          failure_diagnostics$failure_messages,
          sprintf("Sample %d: Evaluation failed - %s", i, 
                  substr(e$message, 1, 100))
        )
        if (verbose) {
          cat("  X Evaluation error:", e$message, "\n")
        }
        return(list(Holdout_MAE = NA, NLL = NA))
      })
      
      # ======================================================================
      # CHECK MAE REASONABLENESS
      # ======================================================================
      # Check if result is valid
      if (is.null(result) || is.na(result$Holdout_MAE) || is.na(result$NLL)) {
        if (is.null(result)) {
          failure_diagnostics$other_failures <<- failure_diagnostics$other_failures + 1
          failure_diagnostics$failure_messages <<- c(
            failure_diagnostics$failure_messages,
            sprintf("Sample %d: Returned NULL result", i)
          )
        }
        return(NULL)
      }

      if (!is.na(result$Holdout_MAE) && !is.null(result$Holdout_MAE)) {
        # Get median of observed dissimilarities for comparison
        # (handles threshold indicators):
        observed_dissim <- extract_numeric_values(local_matrix)
        median_dissim <- median(observed_dissim[!is.na(observed_dissim)])
        
        if (!is.na(median_dissim) && result$Holdout_MAE > (2 * median_dissim)) {
          if (verbose) {
            cat(sprintf("  [!] High MAE (%.3f) vs median dissimilarity (%.3f). Poor fit suspected.\n",
                        result$Holdout_MAE, median_dissim))
          }
        }
      }
      
      # ======================================================================
      # COMPILE RESULTS
      # ======================================================================
      result_row <- data.frame(
        N = N,
        k0 = k0,
        cooling_rate = cooling_rate,
        c_repulsion = c_repulsion,
        Holdout_MAE = result$Holdout_MAE,
        NLL = result$NLL
      )
      
      # Add subsampling metadata if used
      if (use_subsampling) {
        result_row$opt_subsample <- opt_subsample
        result_row$original_n_points <- matrix_size
      }
      
      if (verbose) {
        cat(sprintf("  [OK] MAE: %.4f, NLL: %.2f\n",
                    result$Holdout_MAE, result$NLL))
      }
      
      return(result_row)
      
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("  X Error evaluating sample %d: %s\n", i, e$message))
      }
      return(NULL)
    })
  }
  
  # ==========================================================================
  # PARALLEL EXECUTION WITH BATCHING
  # ==========================================================================
  if (verbose) cat("\nStarting parameter evaluations with batching...\n")
  
  # Determine batch size for memory management
  cores_to_use <- min(max_cores, num_samples)
  batch_size <- min(cores_to_use * 4, num_samples)  # Process 4 per core per batch
  num_batches <- ceiling(num_samples / batch_size)
  
  if (verbose && num_batches > 1) {
    cat(sprintf("Processing in %d batches (batch size: %d)\n", 
                num_batches, batch_size))
  }
  
  all_results <- list()
  
  for (batch in 1:num_batches) {
    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, num_samples)
    batch_indices <- batch_start:batch_end
    
    if (verbose && num_batches > 1) {
      cat(sprintf("\nProcessing batch %d/%d (samples %d-%d)...\n",
                  batch, num_batches, batch_start, batch_end))
    }
    
    # Process batch in parallel
    if (cores_to_use > 1 && length(batch_indices) > 1) {
      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(min(cores_to_use, length(batch_indices)))
        on.exit(parallel::stopCluster(cl), add = TRUE)
        
        # Load topolow on cluster
        parallel::clusterEvalQ(cl, library(topolow))
        
        # Export required objects
        parallel::clusterExport(cl,
                                c("lhs_params", "full_dissimilarity_matrix", "matrix_size",
                                  "mapping_max_iter", "relative_epsilon", "folds",
                                  "use_subsampling", "opt_subsample", "verbose",
                                  "likelihood_function", "euclidean_embedding",
                                  "subsample_dissimilarity_matrix",
                                  "check_matrix_connectivity",
                                  "analyze_network_structure",
                                  "sanity_check_subsample"),
                                envir = environment())
        
        batch_results <- parallel::parLapply(cl, batch_indices, process_param_set)
      } else {
        # Unix-like systems
        batch_results <- parallel::mclapply(batch_indices, process_param_set,
                                            mc.cores = min(cores_to_use, length(batch_indices)))
      }
    } else {
      # Sequential processing for small batches
      batch_results <- lapply(batch_indices, process_param_set)
    }
    
    # Collect batch results
    all_results <- c(all_results, batch_results)
    
    # Garbage collection between batches for memory management
    gc()
    
    if (verbose && num_batches > 1) {
      valid_in_batch <- sum(!sapply(batch_results, is.null))
      cat(sprintf("Batch %d complete: %d/%d successful\n", 
                  batch, valid_in_batch, length(batch_indices)))
    }
  }
  
  # ==========================================================================
  # AGGREGATE RESULTS
  # ==========================================================================
  # Remove NULL results (failed evaluations)
  valid_results <- Filter(Negate(is.null), all_results)
  
  
  if (length(valid_results) == 0) {
    # Provide detailed diagnostic error message
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("ERROR: ALL PARAMETER EVALUATIONS FAILED\n")
    cat(rep("=", 80), "\n\n", sep = "")
    
    cat("Failure Summary:\n")
    cat(sprintf("  - Subsample/connectivity failures: %d/%d (%.1f%%)\n",
                failure_diagnostics$subsample_failures, num_samples,
                100 * failure_diagnostics$subsample_failures / num_samples))
    cat(sprintf("  - CV fold creation failures:       %d/%d (%.1f%%)\n",
                failure_diagnostics$cv_fold_failures, num_samples,
                100 * failure_diagnostics$cv_fold_failures / num_samples))
    cat(sprintf("  - Embedding failures:              %d/%d (%.1f%%)\n",
                failure_diagnostics$embedding_failures, num_samples,
                100 * failure_diagnostics$embedding_failures / num_samples))
    cat(sprintf("  - Other failures:                  %d/%d (%.1f%%)\n",
                failure_diagnostics$other_failures, num_samples,
                100 * failure_diagnostics$other_failures / num_samples))
    
    cat("\nData Characteristics:\n")
    cat(sprintf("  - Matrix size: %d x %d\n", matrix_size, matrix_size))
    cat(sprintf("  - Missing values: %d/%d (%.1f%%)\n",
                sum(is.na(dissimilarity_matrix)),
                matrix_size^2,
                100 * sum(is.na(dissimilarity_matrix)) / matrix_size^2))
    
    if (use_subsampling) {
      cat(sprintf("  - Subsample size: %d (%.1f%% of full matrix)\n",
                  opt_subsample, 100 * opt_subsample / matrix_size))
      cat(sprintf("  - CV folds: %d (each removes ~%d points)\n",
                  folds, floor(opt_subsample / folds)))
    } else {
      cat(sprintf("  - CV folds: %d (each removes ~%d points)\n",
                  folds, floor(matrix_size / folds)))
    }
    
    # Determine primary failure mode
    primary_issue <- "unknown"
    if (failure_diagnostics$subsample_failures > num_samples * 0.5) {
      primary_issue <- "subsampling"
    } else if (failure_diagnostics$cv_fold_failures > num_samples * 0.5) {
      primary_issue <- "cv_folds"
    } else if (failure_diagnostics$embedding_failures > num_samples * 0.5) {
      primary_issue <- "embedding"
    }
    
    cat("\nLikely Root Cause:\n")
    if (primary_issue == "subsampling") {
      cat("  >> DATA TOO SPARSE FOR SUBSAMPLING <<\n\n")
      cat("  Your data has disconnected components when subsampled.\n")
      cat("  This happens when observations are too sparse to maintain\n")
      cat("  connectivity after random sampling.\n\n")
      cat("Recommended Solutions (in order of preference):\n")
      cat("  1. Prune sparse matrix BEFORE optimization:\n")
      cat("       pruned <- prune_sparse_matrix(your_matrix, \n")
      cat("                                      target_completeness = 0.02,\n")
      cat("                                      min_points = 1000)\n")
      cat("       # Then use pruned$pruned_matrix for optimization\n\n")
      cat("  2. Increase opt_subsample (current: ", opt_subsample, ")\n", sep = "")
      cat("       Try: opt_subsample = ", ceiling(opt_subsample * 1.5), "\n\n", sep = "")
      cat("  3. Reduce CV folds (current: ", folds, ")\n", sep = "")
      cat("       Try: folds = ", max(5, floor(folds / 2)), "\n\n", sep = "")
      cat("  4. Disable subsampling (use full matrix):\n")
      cat("       opt_subsample = NULL\n")
      cat("       Warning: This will be slow with large matrices\n\n")
      
    } else if (primary_issue == "cv_folds") {
      cat("  >> INSUFFICIENT DATA FOR CROSS-VALIDATION <<\n\n")
      cat("  The matrix becomes too sparse when splitting into CV folds.\n\n")
      cat("Recommended Solutions:\n")
      cat("  1. Reduce number of CV folds (current: ", folds, ")\n", sep = "")
      cat("       Try: folds = ", max(5, floor(folds / 2)), "\n\n", sep = "")
      cat("  2. Prune sparse matrix to increase density:\n")
      cat("       pruned <- prune_sparse_matrix(your_matrix,\n")
      cat("                                      min_connections = 15)\n\n")
      
    } else if (primary_issue == "embedding") {
      cat("  >> EMBEDDING OPTIMIZATION FAILURES <<\n\n")
      cat("  The optimization process is failing, possibly due to:\n")
      cat("    - Extreme parameter values\n")
      cat("    - Non-metric data structure\n")
      cat("    - Insufficient iterations\n\n")
      cat("Recommended Solutions:\n")
      cat("  1. Check parameter ranges are reasonable\n")
      cat("  2. Increase mapping_max_iter (current: ", mapping_max_iter, ")\n", sep = "")
      cat("  3. Analyze data structure with analyze_network_structure()\n\n")
    }
    
    # Show first few failure messages for debugging
    if (length(failure_diagnostics$failure_messages) > 0) {
      cat("\nFirst Few Failure Messages:\n")
      cat(rep("-", 80), "\n", sep = "")
      n_show <- min(5, length(failure_diagnostics$failure_messages))
      for (msg in failure_diagnostics$failure_messages[1:n_show]) {
        cat("  ", msg, "\n", sep = "")
      }
      if (length(failure_diagnostics$failure_messages) > n_show) {
        cat(sprintf("  ... and %d more failures\n", 
                    length(failure_diagnostics$failure_messages) - n_show))
      }
    }
    
    cat("\n", rep("=", 80), "\n\n", sep = "")
    
    stop("No valid parameter evaluations completed successfully. See diagnostic output above.")
  }

  if (verbose && length(valid_results) < num_samples) {
    cat(sprintf("\nWarning: Only %d/%d parameter evaluations completed successfully.\n",
                length(valid_results), num_samples))
  }
  
  # Combine into data frame
  final_df <- do.call(rbind, valid_results)
  
  # ==========================================================================
  # LOG-TRANSFORM PARAMETERS
  # ==========================================================================
  params_to_transform <- c("N", "k0", "cooling_rate", "c_repulsion")
  
  for (param in params_to_transform) {
    if (!is.numeric(final_df[[param]])) {
      final_df[[param]] <- suppressWarnings(as.numeric(final_df[[param]]))
      if (any(is.na(final_df[[param]]))) {
        stop("Non-numeric values found in column: ", param)
      }
    }
    if (any(final_df[[param]] <= 0, na.rm = TRUE)) {
      stop("Non-positive values found in column: ", param,
           ". Log transform requires positive values.")
    }
  }
  
  # Create log-transformed columns
  final_df$log_N <- log(final_df$N)
  final_df$log_k0 <- log(final_df$k0)
  final_df$log_cooling_rate <- log(final_df$cooling_rate)
  final_df$log_c_repulsion <- log(final_df$c_repulsion)
  
  # Reorder columns: log-transformed params first, then MAE, NLL, then metadata
  base_cols <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion",
                 "Holdout_MAE", "NLL")
  metadata_cols <- setdiff(names(final_df), c(base_cols, params_to_transform))
  final_df <- final_df[, c(base_cols, metadata_cols)]
  
  # ==========================================================================
  # SAVE RESULTS (if requested)
  # ==========================================================================
  if (write_files) {
    output_file <- file.path(param_dir,
                             paste0(scenario_name, "_model_parameters.csv"))
    write.csv(final_df, output_file, row.names = FALSE)
    if (verbose) {
      cat(sprintf("\n[OK] Results saved to: %s\n", output_file))
    }
  }
  
  # ==========================================================================
  # SUMMARY OUTPUT
  # ==========================================================================
  if (verbose) {
    cat("\n" , "=== Optimization Summary ===", "\n")
    cat(sprintf("Successful evaluations: %d/%d\n",
                length(valid_results), num_samples))
    cat(sprintf("Best MAE: %.4f\n", min(final_df$Holdout_MAE, na.rm = TRUE)))
    cat(sprintf("Best NLL: %.2f\n", min(final_df$NLL, na.rm = TRUE)))
    if (use_subsampling) {
      cat(sprintf("Subsampling used: %d points (from %d)\n",
                  opt_subsample, matrix_size))
      
    }
    cat("\n")
  }
  
  return(final_df)
}


#' Run Adaptive Monte Carlo Sampling for Parameter Refinement
#'
#' @description
#' Refines parameter estimates using adaptive Monte Carlo sampling. This function
#' uses initial parameter samples and iteratively generates new samples concentrated
#' in high-likelihood regions of parameter space. Multiple parallel jobs explore
#' different regions simultaneously. Results from all parallel jobs accumulate in a 
#' single output file.
#'
#'
#' @param initial_samples_file Character. Path to a CSV file containing initial samples
#'   (typically from \code{\link{initial_parameter_optimization}}).
#' @param scenario_name Character. Name for the output files.
#' @param dissimilarity_matrix Matrix. The input dissimilarity matrix.
#' @param max_cores Integer. Number of cores to use for parallel execution. If NULL,
#'   uses all available cores minus 1.
#' @param num_samples Integer. Number of new samples to generate via adaptive sampling
#'   per parallel job. Default: 10.
#' @param mapping_max_iter Integer. Maximum number of map optimization iterations.
#' @param relative_epsilon Numeric. Convergence threshold for relative change in error.
#'   Default: 1e-4.
#' @param folds Integer. Number of cross-validation folds. Default: 20.
#' @param opt_subsample Integer or NULL. If specified, randomly samples this many
#'   points from the dissimilarity matrix for each adaptive sampling iteration to
#'   reduce computational cost. The function automatically validates that subsampled
#'   data forms a connected graph. Default: NULL (use full data).
#'
#'   **Important Notes:**
#'   \itemize{
#'     \item If connectivity cannot be achieved, sample size is adaptively increased
#'     \item Each parallel job uses a different subsample
#'     \item Recommended: use same value as in \code{initial_parameter_optimization}
#'     \item The actual subsample size used is reported in output columns
#'   }
#' @param output_dir Character. Required directory for output files.
#' @param verbose Logical. Whether to print progress messages. Default: FALSE.
#'
#' @return No return value. Called for its side effect of writing results to a CSV file
#'   in \code{output_dir}.
#'
#' @details
#' The function runs multiple parallel jobs, each performing adaptive sampling
#' starting from the initial samples. The jobs are independent and can be run
#' simultaneously. Results from all jobs are combined with the initial samples
#' and written to a single output file.
#'
#' If \code{opt_subsample} is used, each parallel job will independently subsample
#' the data, ensuring diversity in the parameter evaluations while maintaining
#' computational efficiency.
#'
#' @examples
#' \donttest{
#' # 1. Locate the example initial samples file included with the package
#' # In a real scenario, this file would be from an 'initial_parameter_optimization' run.
#' initial_file <- system.file(
#'   "extdata", "initial_samples_example.csv",
#'   package = "topolow"
#' )
#'
#' # 2. Create a temporary directory for the function's output
#' # This function requires a writable directory for its results.
#' temp_out_dir <- tempdir()
#'
#' # 3. Create a sample dissimilarity matrix
#' dissim_mat <- matrix(runif(900, 1, 10), 30, 30)
#' diag(dissim_mat) <- 0
#'
#' # 4. Run adaptive sampling (full data)
#' if (nzchar(initial_file)) {
#'   run_adaptive_sampling(
#'     initial_samples_file = initial_file,
#'     scenario_name = "adaptive_test_example",
#'     dissimilarity_matrix = dissim_mat,
#'     output_dir = temp_out_dir,
#'     max_cores = 1,
#'     num_samples = 1,
#'     verbose = FALSE
#'   )
#' }
#'
#' # 5. Run adaptive sampling with subsampling
#' if (nzchar(initial_file)) {
#'   run_adaptive_sampling(
#'     initial_samples_file = initial_file,
#'     scenario_name = "adaptive_test_subsample",
#'     dissimilarity_matrix = dissim_mat,
#'     output_dir = temp_out_dir,
#'     max_cores = 1,
#'     num_samples = 1,
#'     opt_subsample = 15,  # Use only 8 points
#'     verbose = TRUE
#'   )
#' }
#'
#' # Clean up
#' unlink(temp_out_dir, recursive = TRUE)
#' }
#'
#' @seealso \code{\link{adaptive_MC_sampling}} for the underlying sampling algorithm,
#' \code{\link{initial_parameter_optimization}} for generating initial samples.
#'
#' @importFrom utils read.csv write.csv write.table
#' @importFrom future availableCores
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply
#' @export
run_adaptive_sampling <- function(initial_samples_file,
                                  scenario_name,
                                  dissimilarity_matrix,
                                  max_cores = NULL,
                                  num_samples = 10,
                                  mapping_max_iter = 500,
                                  relative_epsilon = 1e-4,
                                  folds = 20,
                                  opt_subsample = NULL,
                                  output_dir,
                                  verbose = FALSE) {
  # ==========================================================================
  # CAPTURE ORIGINAL COLUMN ORDER (CRITICAL FOR CRAN)
  # ==========================================================================
  orig_cols <- names(read.csv(initial_samples_file,
                              stringsAsFactors = FALSE,
                              nrows = 1))
  
  # Parameter names
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  
  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================
  
  # Validate output_dir (CRAN requirement: must be user-specified)
  if (missing(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("'output_dir' must be a character string specifying the directory.", call. = FALSE)
  }
  
  # Validate dissimilarity matrix
  if (!is.matrix(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be a matrix")
  }
  if (nrow(dissimilarity_matrix) != ncol(dissimilarity_matrix)) {
    stop("dissimilarity_matrix must be square")
  }
  matrix_size <- nrow(dissimilarity_matrix)
  
  # Validate initial samples file
  if (!file.exists(initial_samples_file)) {
    stop("initial_samples_file not found: ", initial_samples_file)
  }
  
  # Read and validate initial samples
  init <- read.csv(initial_samples_file, stringsAsFactors = FALSE)
  req <- c(par_names, "NLL")
  if (!all(req %in% names(init))) {
    stop("Missing columns in initial samples: ",
         paste(setdiff(req, names(init)), collapse = ", "))
  }
  
  # Filter valid rows
  init <- init[complete.cases(init[, req]) & 
                 apply(init[, req], 1, function(x) all(is.finite(x))), ]
  if (nrow(init) == 0) {
    stop("No valid initial samples after filtering")
  }
  
  # Validate opt_subsample
  use_subsampling <- FALSE
  if (!is.null(opt_subsample)) {
    if (!is.numeric(opt_subsample) || length(opt_subsample) != 1 ||
        opt_subsample < 10) {
      stop("opt_subsample must be NULL or a numeric value >= 10")
    }
    opt_subsample <- as.integer(floor(opt_subsample))
    
    if (opt_subsample >= matrix_size) {
      if (verbose) {
        cat(sprintf("opt_subsample (%d) >= matrix size (%d). Using full matrix.\n",
                    opt_subsample, matrix_size))
      }
      opt_subsample <- NULL
    } else {
      use_subsampling <- TRUE
      
      # Warn if subsample seems too small
      recommended_min <- max(50, folds)
      if (opt_subsample < recommended_min) {
        warning(sprintf(
          paste0("opt_subsample (%d) is smaller than recommended minimum (%d). ",
                 "This may lead to connectivity issues or unreliable results."),
          opt_subsample, recommended_min
        ), call. = FALSE)
      }
    }
  }
  
  # ==========================================================================
  # PARALLEL PROCESSING SETUP (using future package)
  # ==========================================================================
  
  available_cores <- future::availableCores()
  if (is.null(max_cores)) {
    max_cores <- max(1, available_cores - 1)
  } else {
    max_cores <- min(max_cores, available_cores)
  }
  num_parallel_jobs <- max_cores
  
  # iterations = num_samples per job (following original naming)
  iterations <- num_samples
  
  if (verbose) {
    cat(sprintf("Running %d parallel jobs on %d cores\n",
                num_parallel_jobs, max_cores))
    cat(sprintf("Each job will generate %d new samples\n", iterations))
    cat(sprintf("Total new samples: %d\n", iterations * num_parallel_jobs))
    if (use_subsampling) {
      cat(sprintf("Using subsampling: %d points per job (from %d total)\n",
                  opt_subsample, matrix_size))
    }
  }
  
  # ==========================================================================
  # DIRECTORY SETUP (FOLLOWING ORIGINAL PATTERN)
  # ==========================================================================
  
  # Use original directory naming convention
  adaptive_dir <- file.path(output_dir, "adaptive_sampling_jobs")
  param_dir <- file.path(output_dir, "model_parameters")
  
  # Create directories if they don't exist
  for (dir in c(adaptive_dir, param_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # ==========================================================================
  # COMPREHENSIVE CLEANUP (CRITICAL FOR CRAN)
  # ==========================================================================
  
  if (verbose) cat("Cleaning adaptive sampling jobs directory...\n")
  
  # List ALL files in the directory and remove them
  files_to_remove <- list.files(adaptive_dir, full.names = TRUE, recursive = TRUE)
  if (length(files_to_remove) > 0) {
    file.remove(files_to_remove)
  }
  
  # ==========================================================================
  # SETUP RESULT FILE
  # ==========================================================================
  
  results_file <- file.path(param_dir,
                            paste0(scenario_name, "_model_parameters.csv"))
  
  # Copy initial samples to final destination if it doesn't exist
  if (!file.exists(results_file)) {
    file.copy(initial_samples_file, results_file)
  }
  
  # ==========================================================================
  # CREATE TEMPORARY FILES FOR PARALLEL JOBS (ORIGINAL PATTERN)
  # ==========================================================================
  
  # Use original naming pattern: job_XX_scenario.csv
  make_temp <- function(i) {
    file.path(adaptive_dir, sprintf("job_%02d_%s.csv", i, scenario_name))
  }
  temps <- vapply(seq_len(num_parallel_jobs), make_temp, FUN.VALUE = "")
  
  # Copy initial samples to each temp file
  for (i in seq_along(temps)) {
    file.copy(initial_samples_file, temps[i], overwrite = TRUE)
  }
  
  # ==========================================================================
  # RUN PARALLEL ADAPTIVE SAMPLING JOBS
  # ==========================================================================
  
  if (verbose) cat("\nStarting parallel adaptive sampling jobs...\n")
  
  # ---------------- Local parallel execution ----------------
  # Following original pattern for Windows vs Unix
  
  if (.Platform$OS.type == "windows") {
    # Windows: use PSOCK cluster
    cl <- parallel::makeCluster(num_parallel_jobs)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Load topolow package in each worker
    parallel::clusterEvalQ(cl, library(topolow))
    
    # Export necessary variables (following original pattern)
    parallel::clusterExport(cl,
                            c("adaptive_MC_sampling", "temps", "dissimilarity_matrix",
                              "mapping_max_iter", "relative_epsilon", "folds",
                              "output_dir", "scenario_name", "iterations",
                              "opt_subsample","subsample_dissimilarity_matrix",
                          "check_matrix_connectivity",
                          "analyze_network_structure",
                          "sanity_check_subsample"), 
                            envir = environment())
    
    # Run jobs
    parallel::parLapply(cl, seq_along(temps), function(i) {
      adaptive_MC_sampling(
        samples_file = temps[i],
        dissimilarity_matrix = dissimilarity_matrix,
        iterations = iterations,
        mapping_max_iter = mapping_max_iter,
        relative_epsilon = relative_epsilon,
        folds = folds,
        num_cores = 1,
        scenario_name = scenario_name,
        opt_subsample = opt_subsample,  # Pass subsampling parameter
        verbose = FALSE  # Keep individual jobs quiet
      )
    })
    
    parallel::stopCluster(cl)
    
  } else {
    # For Unix-like systems (Linux, macOS)
    parallel::mclapply(seq_along(temps), function(i) {
      adaptive_MC_sampling(
        samples_file = temps[i],
        dissimilarity_matrix = dissimilarity_matrix,
        iterations = iterations,
        mapping_max_iter = mapping_max_iter,
        relative_epsilon = relative_epsilon,
        folds = folds,
        num_cores = 1,
        scenario_name = scenario_name,
        opt_subsample = opt_subsample,  # Pass subsampling parameter
        verbose = FALSE  # Keep individual jobs quiet
      )
    }, mc.cores = num_parallel_jobs)
  }
  
  # ==========================================================================
  # GATHER AND COMBINE RESULTS (FOLLOWING ORIGINAL PATTERN)
  # ==========================================================================
  
  if (verbose) cat("\nLocal parallel jobs complete; gathering results...\n")
  
  # Read initial samples to get original row count
  init2 <- read.csv(initial_samples_file, stringsAsFactors = FALSE)
  n0 <- nrow(init2)
  
  # Read new samples from each temp file
  new_list <- lapply(temps, function(f) {
    if (file.exists(f)) {
      df <- tryCatch(
        read.csv(f, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      # Extract only new rows (beyond initial samples)
      if (!is.null(df) && nrow(df) > n0) {
        df[(n0 + 1):nrow(df), , drop = FALSE]
      } else {
        NULL
      }
    } else {
      NULL
    }
  })
  
  # Combine all results
  all <- do.call(rbind, c(list(init2), new_list))
  
  # ==========================================================================
  # REORDER TO ORIGINAL INPUT HEADER (CRITICAL FOR CRAN)
  # ==========================================================================
  
  # Ensure all original columns are present before reordering
  if (all(orig_cols %in% names(all))) {
    all <- all[, orig_cols, drop = FALSE]
  }
  
  # Write final results
  write.csv(all, results_file, row.names = FALSE)
  
  # ==========================================================================
  # CLEANUP TEMPORARY FILES
  # ==========================================================================
  
  file.remove(temps)
  
  if (verbose) {
    cat("Results written to:", results_file, "\n")
    cat(sprintf("Total samples: %d (initial: %d, new: %d)\n",
                nrow(all), n0, nrow(all) - n0))
    
    if (use_subsampling) {
      cat(sprintf("Subsampling was used: %d points per job\n", opt_subsample))

    }
    
    # Show best results
    if ("Holdout_MAE" %in% names(all)) {
      best_mae <- min(all$Holdout_MAE, na.rm = TRUE)
      cat(sprintf("Best MAE achieved: %.4f\n", best_mae))
    }
    if ("NLL" %in% names(all)) {
      best_nll <- min(all$NLL, na.rm = TRUE)
      cat(sprintf("Best NLL achieved: %.2f\n", best_nll))
    }
  }
  
  return(invisible(NULL))
}


#' Perform Adaptive Monte Carlo Sampling (Internal)
#'
#' @description
#' Core implementation of the adaptive Monte Carlo sampling algorithm. This internal
#' function explores the parameter space by updating the sampling distribution
#' based on evaluated likelihoods. It is called by the main \code{run_adaptive_sampling}
#' function.
#'
#' @param samples_file Path to the CSV file with samples for the current job.
#' @param dissimilarity_matrix The dissimilarity matrix to be fitted.
#' @param iterations Number of sampling iterations per job.
#' @param mapping_max_iter Maximum optimization iterations for the embedding.
#' @param relative_epsilon Convergence threshold for the optimization.
#' @param folds Number of cross-validation folds.
#' @param num_cores Number of cores for parallel processing within this job.
#' @param scenario_name Name for output files, used for context.
#' @param opt_subsample Integer or NULL. If specified, subsamples the data for
#'   this adaptive sampling job. A single subsample is created at the start and
#'   reused for all iterations within this job for consistency. Default: NULL.
#' @param verbose Logical. If TRUE, prints progress messages. Default: FALSE.
#'
#' @return A \code{data.frame} containing all samples (initial and newly generated)
#' with their parameters and evaluated performance metrics. The data frame includes
#' columns for the log-transformed parameters, \code{Holdout_MAE}, and \code{NLL}.
#' Returns \code{NULL} if the results file was not created.
#'
#' @details
#' The function performs these steps:
#' 1. Reads existing samples from the samples file
#' 2. If \code{opt_subsample} is specified, creates a single subsample for this job
#' 3. Calculates weighted marginal distributions for each parameter
#' 4. For each iteration:
#'    - Samples new parameter values from the marginal distributions
#'    - Evaluates the parameter set via cross-validation
#'    - Appends results to the samples file
#'    - Updates the marginal distributions
#'
#' The adaptive aspect comes from using likelihood-weighted KDE on previous samples
#' to concentrate new samples in high-likelihood regions.
#'
#' **Subsampling behavior**: When \code{opt_subsample} is used, a single subsample
#' is created at the start of this job and reused for all iterations. This ensures
#' consistency within the job while still allowing different jobs (in
#' \code{run_adaptive_sampling}) to explore different subsamples.
#'
#' @keywords internal
#' @importFrom filelock lock unlock
#' @importFrom utils read.csv write.table
#' @importFrom stats na.omit
#' @export
adaptive_MC_sampling <- function(samples_file,
                                 dissimilarity_matrix,
                                 iterations = 1,
                                 mapping_max_iter,
                                 relative_epsilon,
                                 folds = 20,
                                 num_cores = 1,
                                 scenario_name,
                                 opt_subsample = NULL,
                                 verbose = FALSE) {
  # ==========================================================================
  # INPUT VALIDATION AND SETUP
  # ==========================================================================
  
  # Require filelock for safe concurrent writes
  if (!requireNamespace("filelock", quietly = TRUE)) {
    stop("Package 'filelock' is required for safe concurrent writes. Please install it.")
  }
  
  if (!file.exists(samples_file)) {
    stop("samples_file does not exist: ", samples_file)
  }
  
  # Store reference to full matrix
  full_dissimilarity_matrix <- dissimilarity_matrix
  matrix_size <- nrow(dissimilarity_matrix)
  
  # Validate opt_subsample
  use_subsampling <- FALSE
  
  if (!is.null(opt_subsample)) {
    if (!is.numeric(opt_subsample) || length(opt_subsample) != 1 ||
        opt_subsample < 10) {
      stop("opt_subsample must be NULL or a numeric value >= 10")
    }
    opt_subsample <- as.integer(floor(opt_subsample))
    
    if (opt_subsample >= matrix_size) {
      if (verbose) {
        cat(sprintf("opt_subsample (%d) >= matrix size (%d). Using full matrix.\n",
                    opt_subsample, matrix_size))
      }
      opt_subsample <- NULL
    } else {
      use_subsampling <- TRUE
    }
  }
  
  # ==========================================================================
  # CREATE SUBSAMPLE (if enabled) - ONCE FOR THIS JOB
  # ==========================================================================
  
  if (use_subsampling) {
    if (verbose) {
      cat(sprintf("Creating subsample for this adaptive job: %d points from %d\n",
                  opt_subsample, matrix_size))
    }
    
    subsample_result <- subsample_dissimilarity_matrix(
      dissimilarity_matrix = full_dissimilarity_matrix,
      sample_size = opt_subsample,
      max_attempts = 5,
      verbose = verbose
    )
    
    if (!subsample_result$is_connected) {
      stop(sprintf(
        paste0("Failed to obtain connected subsample for adaptive sampling job '%s'.\n",
               "This job will be skipped."),
        scenario_name
      ))
    }
    
    # Use this subsample for all iterations in this job
    dissimilarity_matrix <- subsample_result$subsampled_matrix
    
    if (verbose) {
      cat(sprintf("  Using this subsample for all %d iterations\n", iterations))
    }
    
    # Sanity check
    sanity_result <- sanity_check_subsample(
      dissimilarity_matrix,
      folds = folds,
      verbose = verbose
    )
    
    if (!sanity_result$all_checks_passed && verbose) {
      cat("  [!] Subsample sanity checks raised warnings (proceeding anyway)\n")
    }
  }
  
  # ==========================================================================
  # ADAPTIVE SAMPLING LOOP
  # ==========================================================================
  
  if (verbose) {
    cat(sprintf("\nStarting %d adaptive sampling iterations for job '%s'\n",
                iterations, scenario_name))
  }
  
  for (iter in 1:iterations) {
    if (verbose) {
      cat(sprintf("\n--- Iteration %d/%d ---\n", iter, iterations))
    }
    
    # ========================================================================
    # READ CURRENT SAMPLES
    # ========================================================================
    samples <- tryCatch({
      read.csv(samples_file, stringsAsFactors = FALSE)
    }, error = function(e) {
      stop("Failed to read samples file: ", e$message)
    })
    
    # Validate required columns
    required_cols <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion", "NLL")
    if (!all(required_cols %in% names(samples))) {
      stop("Samples file missing required columns: ",
           paste(setdiff(required_cols, names(samples)), collapse = ", "))
    }
    
    # ========================================================================
    # CALCULATE WEIGHTED MARGINALS
    # ========================================================================
    if (verbose) cat("  Calculating weighted marginals...\n")
    
    marginals <- tryCatch({
      calculate_weighted_marginals(samples)
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("  Warning: Failed to calculate marginals: %s\n", e$message))
      }
      return(NULL)
    })
    
    if (is.null(marginals)) {
      if (verbose) cat("  Skipping iteration due to marginal calculation failure\n")
      next
    }
    
    # ========================================================================
    # SAMPLE NEW PARAMETERS
    # ========================================================================
    if (verbose) cat("  Sampling new parameter values...\n")
    
    # Sample from each marginal distribution
    log_N_new <- sample(marginals$log_N$x, size = 1,
                        prob = marginals$log_N$y)
    log_k0_new <- sample(marginals$log_k0$x, size = 1,
                         prob = marginals$log_k0$y)
    log_cooling_rate_new <- sample(marginals$log_cooling_rate$x, size = 1,
                                   prob = marginals$log_cooling_rate$y)
    log_c_repulsion_new <- sample(marginals$log_c_repulsion$x, size = 1,
                                  prob = marginals$log_c_repulsion$y)
    
    # Convert back to original scale for embedding
    N_new <- round(exp(log_N_new))
    k0_new <- exp(log_k0_new)
    cooling_rate_new <- exp(log_cooling_rate_new)
    c_repulsion_new <- exp(log_c_repulsion_new)
    
    if (verbose) {
      cat(sprintf("  New params: N=%d, k0=%.3f, cool=%.4f, c_rep=%.4f\n",
                  N_new, k0_new, cooling_rate_new, c_repulsion_new))
    }
    
    # ========================================================================
    # EVALUATE NEW PARAMETERS
    # ========================================================================
    if (verbose) cat("  Evaluating via cross-validation...\n")
    
    result <- tryCatch({
      likelihood_function(
        dissimilarity_matrix = dissimilarity_matrix,  # Uses subsample if set
        mapping_max_iter = mapping_max_iter,
        relative_epsilon = relative_epsilon,
        N = N_new,
        k0 = k0_new,
        cooling_rate = cooling_rate_new,
        c_repulsion = c_repulsion_new,
        folds = folds,
        num_cores = num_cores
      )
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("  X Evaluation failed: %s\n", e$message))
      }
      return(NULL)
    })
    
    if (is.null(result) || is.na(result$Holdout_MAE) || is.na(result$NLL)) {
      if (verbose) cat("  X Evaluation returned invalid results. Skipping.\n")
      next
    }
    
    # ========================================================================
    # CHECK MAE REASONABLENESS
    # ========================================================================
    #observed_dissim <- as.numeric(dissimilarity_matrix)
    observed_dissim <- extract_numeric_values(dissimilarity_matrix)
    median_dissim <- median(observed_dissim[!is.na(observed_dissim)])
    
    if (!is.na(median_dissim) && result$Holdout_MAE > (2 * median_dissim)) {
      if (verbose) {
        cat(sprintf("  [!] High MAE (%.3f) vs median dissimilarity (%.3f)\n",
                    result$Holdout_MAE, median_dissim))
      }
    }
    
    if (verbose) {
      cat(sprintf("  [OK] MAE: %.4f, NLL: %.2f\n",
                  result$Holdout_MAE, result$NLL))
    }
    
    # ========================================================================
    # APPEND NEW SAMPLE TO FILE
    # ========================================================================
    
    # Create new row with log-transformed parameters
    new_row <- data.frame(
      log_N = log_N_new,
      log_k0 = log_k0_new,
      log_cooling_rate = log_cooling_rate_new,
      log_c_repulsion = log_c_repulsion_new,
      Holdout_MAE = result$Holdout_MAE,
      NLL = result$NLL
    )
    
    # Add subsampling metadata if used
    if (use_subsampling) {
      new_row$opt_subsample <- opt_subsample
      new_row$original_n_points <- matrix_size
    }
    
    # Match columns to existing samples
    for (col in names(samples)) {
      if (!(col %in% names(new_row))) {
        new_row[[col]] <- NA
      }
    }
    new_row <- new_row[, names(samples), drop = FALSE]
    
    # Append to file with file locking for thread safety
    tryCatch({
      # Use filelock if available for thread-safe writing
      if (requireNamespace("filelock", quietly = TRUE)) {
        lock_file <- paste0(samples_file, ".lock")
        lock <- filelock::lock(lock_file, timeout = 30000)
        on.exit({
          if (!is.null(lock)) filelock::unlock(lock)
          if (file.exists(lock_file)) file.remove(lock_file)
        }, add = TRUE)
      }
      
      # Append to CSV
      write.table(new_row, file = samples_file,
                  sep = ",", append = TRUE,
                  row.names = FALSE, col.names = FALSE)
      
      if (verbose) cat("  [OK] Sample appended to file\n")
      
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("  X Failed to append sample: %s\n", e$message))
      }
    })
  }
  
  # ==========================================================================
  # FINALIZATION
  # ==========================================================================
  
  if (verbose) {
    cat(sprintf("\n[OK] Adaptive sampling complete for job '%s'\n", scenario_name))
    cat(sprintf("  %d iterations completed\n", iterations))

  }
  
  # Read and return final results
  final_samples <- tryCatch({
    read.csv(samples_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    if (verbose) cat("Warning: Could not read final samples file\n")
    return(NULL)
  })
  
  return(final_samples)
}



#' Generate New Parameter Samples Using KDE
#'
#' @description
#' Generates new parameter samples using weighted kernel density estimation
#' for each parameter independently. This is an internal helper function for the
#' adaptive sampling process.
#'
#' @param samples A data frame of previous samples containing parameter columns and an "NLL" column.
#' @param n The integer number of new samples to generate.
#' @param epsilon A numeric probability (0-1) of sampling with a wider bandwidth to
#'   encourage exploration. Default is 0.
#'
#' @return A data frame containing `n` new parameter samples.
#'
#' @keywords internal
#' @importFrom stats na.omit runif bw.nrd0 approx
generate_kde_samples <- function(samples, n, epsilon = 0) {
  # First, remove outliers from the input samples to stabilize the KDE
  samples <- as.data.frame(lapply(samples, clean_data, k = 3))
  samples <- na.omit(samples)

  # Check if we have enough samples after cleaning
  if (nrow(samples) < 2) {
    stop("Insufficient samples remaining after removing NA & outliers (need at least 2)")
  }

  # Calculate importance weights from the Negative Log-Likelihood (NLL)
  log_likes <- -samples$NLL
  # Shift log-likelihoods to be positive for stability
  std_log_likes <- log_likes - min(log_likes) + 0.05
  weights <- std_log_likes / sum(std_log_likes)

   # Check for NA weights
  if (any(is.na(weights))) {
    warning("NA values in weights, replacing with uniform weights")
    weights <- rep(1/nrow(samples), nrow(samples))
  }

  # Define the names of the parameters to be sampled
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")

  # Initialize a data frame to store the new samples
  new_samples <- data.frame(matrix(nrow = n, ncol = length(par_names)))
  names(new_samples) <- par_names

  # Generate samples for each parameter independently
  for (param in par_names) {
    # With a small probability (epsilon), use a wider bandwidth to escape local optima
    if (runif(1) < epsilon) {
      bandwidth <- bw.nrd0(samples[[param]]) * 2
    } else {
      bandwidth <- NULL  # Use the default bandwidth in the KDE function
    }

    # Get the weighted kernel density estimate for the current parameter
    kde <- weighted_kde(samples[[param]], weights)

     # Check for NA or negative values in density estimates
    if (any(is.na(kde$y)) || any(kde$y < 0)) {
      warning(paste("Invalid KDE values for parameter", param, "- using uniform sampling"))
      # Fallback to uniform sampling from the range
      new_samples[[param]] <- runif(n, min(samples[[param]]), max(samples[[param]]))
      next
    }

    # Sample from the KDE using inverse transform sampling
    u <- runif(n)
    cdf <- cumsum(kde$y) / sum(kde$y)

    # Final check for NA in CDF
    if (any(is.na(cdf))) {
      warning(paste("NA values in CDF for parameter", param, "- using uniform sampling"))
      new_samples[[param]] <- runif(n, min(samples[[param]]), max(samples[[param]]))
      next
    }

    # Use linear interpolation to find the sample values from the CDF
    new_samples[[param]] <- approx(cdf, kde$x, u, rule = 2)$y
  }

  return(new_samples)
}


#' Weighted Kernel Density Estimation
#'
#' @description
#' Performs weighted kernel density estimation for univariate data. This is useful for
#' analyzing parameter distributions where each sample has an associated importance
#' weight (e.g., a likelihood).
#'
#' @param x A numeric vector of samples.
#' @param weights A numeric vector of weights corresponding to each sample in x.
#' @param n The integer number of points at which to evaluate the density.
#' @param from,to The range over which to evaluate the density.
#'
#' @return A list containing the evaluation points (`x`) and the estimated density values (`y`).
#'
#' @importFrom stats sd dnorm
#' @importFrom parallel detectCores mclapply
#' @export
weighted_kde <- function(x, weights, n = 512, from = min(x), to = max(x)) {
  # Normalize weights to ensure they sum to 1
  weights <- weights / sum(weights)

  # Calculate bandwidth using Silverman's rule of thumb
  bw <- 1.06 * sd(x) * length(x)^(-1/5)
  eval_points <- seq(from, to, length.out = n)

  # Define the function to compute density at a single point z
  compute_density <- function(z) {
    sum(weights * dnorm(z, mean = x, sd = bw))
  }

  # --- SAFE PARALLELISM CONTROL ---
  # In non-interactive sessions (like R CMD check), limit cores to 2.
  num_cores <- if (!interactive()) {
    2L
  } else {
    getOption("mc.cores", max(1L, future::availableCores() - 1L))
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



#' Profile Likelihood Analysis
#'
#' @description
#' Calculates the profile likelihood for a given parameter by evaluating the conditional
#' maximum likelihood across a grid of parameter values. This "empirical profile
#' likelihood" estimates the likelihood surface based on samples from Monte Carlo
#' simulations.
#'
#' @details
#' For each value in the parameter grid, the function:
#' 1. Identifies nearby samples using a bandwidth window.
#' 2. Calculates the conditional maximum likelihood from these samples.
#' 3. Tracks sample counts to assess the reliability of the estimate.
#'
#' @param param The character name of the parameter to analyze (e.g., "log_N").
#' @param samples A data frame containing parameter samples and a log-likelihoods column named "NLL".
#' @param grid_size The integer number of grid points for the analysis.
#' @param bandwidth_factor A numeric factor for the local sample window size.
#' @param start_factor,end_factor Numeric range multipliers for parameter grid (default: 0.5, 1.2)
#' @param min_samples Integer minimum samples required for reliable estimate (default: 10)
#' @return Object of class "profile_likelihood" containing:
#'   \item{param}{Vector of parameter values}
#'   \item{ll}{Vector of log-likelihood values}
#'   \item{param_name}{Name of analyzed parameter}
#'   \item{bandwidth}{Bandwidth used for local windows}
#'   \item{sample_counts}{Number of samples per estimate}
#' @examples
#' # Create a sample data frame of parameter samples
#' mcmc_samples <- data.frame(
#'   log_N = log(runif(50, 2, 10)),
#'   log_k0 = log(runif(50, 1, 5)),
#'   log_cooling_rate = log(runif(50, 0.01, 0.1)),
#'   log_c_repulsion = log(runif(50, 0.1, 1)),
#'   NLL = runif(50, 20, 100)
#' )
#'
#' # Calculate profile likelihood for the "log_N" parameter
#' pl <- profile_likelihood("log_N", mcmc_samples,
#'                          grid_size = 10, # Smaller grid for a quick example
#'                          bandwidth_factor = 0.05)
#'
#' # Print the results
#' print(pl)
#'
#' @seealso
#' The S3 methods `print.profile_likelihood` and `summary.profile_likelihood` for viewing results.
#'
#' @importFrom stats sd approx
#' @export
profile_likelihood <- function(param, samples, grid_size = 40,
                               bandwidth_factor = 0.05,
                               start_factor = 0.5, end_factor = 1.5,
                               min_samples = 5) {

  # --- Input Validation ---
  if (!is.character(param) || length(param) != 1) {
    stop("param must be a single character string.")
  }
  if (!param %in% names(samples)) {
    stop(sprintf("Parameter '%s' not found in samples.", param))
  }
  if (!("NLL" %in% names(samples))) {
    stop("Samples data frame must contain an 'NLL' column.")
  }
  if (grid_size < 2) stop("grid_size must be at least 2.")
  if (bandwidth_factor <= 0) stop("bandwidth_factor must be positive.")
    if (start_factor >= end_factor) {
    stop("start_factor must be less than end_factor")
  }
  # --- Get Parameter Grid ---
  # Generate a grid of values centered around the parameter's MLE.
  grid_values <- get_grid(samples, param, grid_size, start_factor, end_factor)

  # Initialize result vectors
  ll_values <- numeric(length(grid_values))
  sample_counts <- numeric(length(grid_values))

  # --- Calculate Profile Likelihood ---
  for (i in seq_along(grid_values)) {
    val <- grid_values[i]

    # Define a window around the current grid value to get conditional samples
    param_range <- sd(samples[[param]], na.rm = TRUE) * bandwidth_factor
    conditional_samples <- samples[abs(samples[[param]] - val) <= param_range, ]
    n_samples <- nrow(conditional_samples)

    # Ensure there are enough samples in the window for a reliable estimate
    if (n_samples < min_samples) {
      warning(sprintf("Too few samples (%d) near %s = %.3f. Consider a larger bandwidth_factor.",
                      n_samples, param, val))
      ll_values[i] <- NA
    } else {
      # Calculate max likelihood from the top n=3 samples in the window using log-sum-exp trick for stability
      x <- -conditional_samples$NLL[order(conditional_samples$NLL)[1:min(3, n_samples)]]
      max_x <- max(x, na.rm = TRUE)
      ll_values[i] <- max_x + log(sum(exp(x - max_x), na.rm = TRUE)) - log(length(x))
    }
    sample_counts[i] <- n_samples
  }

  # Handle case where all likelihood values are NA
  if (all(is.na(ll_values))) {
    warning("All likelihood values are NA. Consider increasing bandwidth_factor or using more samples.")
    # Set all values to a default
    ll_values[] <- -Inf
  }

  # Interpolate NA values caused by sparse regions in the parameter space
  na_idx <- which(is.na(ll_values))
  if (length(na_idx) > 0 && length(na_idx) < length(ll_values)) {
    # Only interpolate if we have at least 2 non-NA values
    non_na_count <- sum(!is.na(ll_values))
    if (non_na_count >= 2) {
      ll_values[na_idx] <- approx(grid_values[-na_idx],
                                  ll_values[-na_idx],
                                  grid_values[na_idx])$y
    } else {
      # If insufficient non-NA values, set all to the mean of available values
      if (non_na_count == 1) {
        ll_values[na_idx] <- ll_values[!is.na(ll_values)][1]
      }
      # If no non-NA values, they remain NA
    }
  }

  # Create result object using the S3 constructor
  result <- profile_likelihood_result(
    param_values = grid_values,
    ll_values = ll_values,
    param_name = param,
    bandwidth = param_range,
    sample_counts = sample_counts
  )

  return(result)
}


#' Create Grid Around Maximum Likelihood Estimate (Internal)
#'
#' @description
#' Internal helper to generate a sequence of values for a parameter. The grid is
#' centered on the parameter's Maximum Likelihood Estimate (MLE), which is found
#' by calculating the mode of its weighted marginal distribution.
#'
#' @param samples Data frame of parameter samples with an NLL column.
#' @param param Character name of the parameter column.
#' @param num_points Integer number of points for the grid.
#' @param start_factor Numeric factor for grid's lower boundary relative to MLE.
#' @param end_factor Numeric factor for grid's upper boundary relative to MLE.
#' @return A numeric vector of grid points.
#' @keywords internal
get_grid <- function(samples, param, num_points, start_factor, end_factor) {
  # Validate inputs
  if (!is.data.frame(samples) || !all(c(param, "NLL") %in% names(samples))) {
    stop("`samples` must be a data frame with the specified parameter and NLL columns.")
  }
  if (start_factor >= end_factor) {
    stop("`start_factor` must be less than `end_factor`.")
  }

  # Validate NLL values
  if (all(is.infinite(samples$NLL))) {
    stop("All NLL values are infinite")
  }
  if (any(is.na(samples$NLL))) {
    warning("NA values in the NLL column will be removed.")
    samples <- samples[!is.na(samples$NLL), ]
  }

  if (nrow(samples) < 2) {
    stop("At least two valid samples are required after filtering.")
  }

  # --- Calculate weighted marginal for the specific parameter ---
  # Remove outliers to stabilize density estimation
  clean_samples <- as.data.frame(lapply(samples[, c(param, "NLL")], clean_data, k = 3))
  clean_samples <- na.omit(clean_samples)

  if (nrow(clean_samples) < 2) {
    stop("Insufficient samples remaining after cleaning for parameter: ", param)
  }

  # Calculate importance weights from log-likelihoods
  log_likelihoods <- -clean_samples$NLL
  # Shift to be positive and normalize
  std_log_likelihoods <- log_likelihoods - min(log_likelihoods, na.rm = TRUE) + 0.05
  weights <- std_log_likelihoods / sum(std_log_likelihoods, na.rm = TRUE)

  # Calculate weighted marginal for the specific parameter
  marginal <- weighted_kde(clean_samples[[param]], weights = weights)

  # Find the mode (MLE) for the specific parameter
  MLE <- marginal$x[which.max(marginal$y)]

  # Create the grid sequence centered around the MLE
  seq(from = start_factor * MLE, to = end_factor * MLE, length.out = num_points)
}


#' S3 Constructor for Profile Likelihood Results
#'
#' @description
#' Internal S3 constructor for storing results from `profile_likelihood`.
#'
#' @param param_values Vector of parameter values.
#' @param ll_values Vector of log-likelihood values.
#' @param param_name Name of the analyzed parameter.
#' @param bandwidth Bandwidth used for local windows.
#' @param sample_counts Number of samples per estimate.
#' @return An object of class `profile_likelihood`.
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


#' Parameter Sensitivity Analysis
#'
#' @description
#' Analyzes the sensitivity of the model performance (measured by MAE) to changes
#' in a single parameter. This function bins the parameter range to identify the
#' minimum MAE for each bin, helping to understand how robust the model is to
#' parameter choices.
#'
#' @details
#' The function performs these steps:
#' 1. Cleans the input data using Median Absolute Deviation (MAD) to remove outliers.
#' 2. Bins the parameter values into equal-width bins.
#' 3. Calculates the minimum MAE within each bin to create an empirical performance curve.
#' 4. Identifies a performance threshold based on a percentage above the global minimum MAE.
#' 5. Returns an S3 object for plotting and further analysis.
#'
#' @param param The character name of the parameter to analyze.
#' @param samples A data frame containing parameter samples and performance metrics.
#' @param bins The integer number of bins to divide the parameter range into.
#' @param mae_col The character name of the column containing the Mean Absolute Error (MAE) values.
#' @param threshold_pct A numeric percentage above the minimum MAE to define an acceptable performance threshold.
#' @param min_samples The integer minimum number of samples required in a bin for it to be included in the analysis.
#'
#' @return An object of class "parameter_sensitivity" containing:
#'   \item{param_values}{Vector of parameter bin midpoints}
#'   \item{min_mae}{Vector of minimum MAE values per bin}
#'   \item{param_name}{Name of analyzed parameter}
#'   \item{threshold}{Threshold value (default: min. +5%)}
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
  # --- Input Validation ---
  if (!is.character(param) || length(param) != 1) {
    stop("`param` must be a single character string.")
  }
  if (!param %in% names(samples)) {
    stop(sprintf("Parameter column '%s' not found in samples.", param))
  }
  if (!mae_col %in% names(samples)) {
    stop(sprintf("MAE column '%s' not found in samples.", mae_col))
  }
  if (bins < 2) {
    stop("`bins` must be at least 2.")
  }

  # --- Data Cleaning ---
  # Clean the relevant columns using MAD-based outlier detection for robustness
  clean_samples <- as.data.frame(lapply(samples[, c(param, mae_col)], clean_data, k = 3))
  clean_samples <- na.omit(clean_samples)

  if (nrow(clean_samples) < min_samples * 2) {
    stop("Insufficient data remaining for analysis after cleaning.")
  }

  # --- Binning and Analysis ---
  # Find the range of parameter values to create bins
  min_param <- min(clean_samples[[param]])
  max_param <- max(clean_samples[[param]])

  # Create breaks for the bins, adding a small epsilon to include the max value
  bin_width <- (max_param - min_param) / bins
  breaks <- seq(min_param, max_param + 1e-10, by = bin_width)

  # Use hist to get bin assignments and counts without plotting
  h <- hist(clean_samples[[param]], breaks = breaks, plot = FALSE)

  # Initialize a data frame to store the minimum MAE for each bin
  bin_min_mae <- data.frame(
    bin_midpoint = h$mids,
    min_mae = NA_real_,
    sample_count = h$counts
  )

  # Iterate through each bin to find the minimum MAE
  for (i in 1:bins) {
    # Isolate the data points that fall within the current bin
    bin_data <- clean_samples[clean_samples[[param]] >= breaks[i] &
                              clean_samples[[param]] < breaks[i+1], ]

    if (nrow(bin_data) >= min_samples) {
      bin_min_mae$min_mae[i] <- min(bin_data[[mae_col]], na.rm = TRUE)
    }
  }

  # Remove bins that had insufficient data for a reliable calculation
  plot_data <- bin_min_mae[!is.na(bin_min_mae$min_mae), ]

  # Calculate the global minimum MAE and the performance threshold
  min_value <- min(plot_data$min_mae)
  threshold_value <- min_value * (1 + threshold_pct / 100)

  # --- Create S3 Result Object ---
  result <- structure(
    list(
      param_values = plot_data$bin_midpoint,
      min_mae = plot_data$min_mae,
      param_name = param,
      threshold = threshold_value,
      min_value = min_value,
      sample_counts = plot_data$sample_count
    ),
    class = "parameter_sensitivity"
  )

  return(result)
}


#' Plot Parameter Sensitivity Analysis
#'
#' @description
#' The S3 plot method for `parameter_sensitivity` objects. It creates a visualization
#' showing how the model's performance (minimum MAE) changes across the range of a
#' single parameter. A threshold line is included to indicate the region of acceptable
#' performance.
#'
#' @param x A `parameter_sensitivity` object, typically from `parameter_sensitivity_analysis()`.
#' @param width The numeric width of the output plot in inches.
#' @param height The numeric height of the output plot in inches.
#' @param save_plot A logical indicating whether to save the plot to a file.
#' @param output_dir A character string specifying the directory for output files. Required if `save_plot` is TRUE.
#' @param y_limit_factor A numeric factor to set the upper y-axis limit as a percentage
#'   above the threshold value (e.g., 1.10 for 10% above). If NULL, scaling is automatic.
#' @param ... Additional arguments (not currently used).
#'
#' @return A `ggplot` object representing the sensitivity plot.
#'
#' @method plot parameter_sensitivity
#' @importFrom ggplot2 ggplot aes geom_line geom_hline annotate labs theme_minimal theme element_text element_line element_blank element_rect margin scale_y_continuous ggsave
#' @export
plot.parameter_sensitivity <- function(x, width = 3.5, height = 3.5,
                                       save_plot = FALSE, output_dir,
                                       y_limit_factor = NULL, ...) {
  # check if required packages are installed
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for formatting y-axis labels. Please install it.")
  }

  # Convert the S3 object to a data frame for ggplot2
  plot_data <- data.frame(
    param = x$param_values,
    min_mae = x$min_mae
  )

  # Convert internal parameter names to mathematical expressions for cleaner labels
  param_expr <- switch(x$param_name,
                       "log_N" = expression(log(N)),
                       "log_c_repulsion" = expression(log(c)),
                       "log_cooling_rate" = expression(log(alpha)),
                       "log_k0" = expression(log(k)),
                       x$param_name) # Default to the original name if not matched

  # Calculate position for the threshold label for consistent placement
  label_x_pos <- min(plot_data$param) + 0.03 * (max(plot_data$param) - min(plot_data$param))

  # Create the plot object
  p <- ggplot() +
    # Main data line showing minimum MAE across the parameter range
    geom_line(data = plot_data,
              aes(x = .data$param, y = .data$min_mae),
              color = "steelblue",
              linewidth = 0.8) +
    # Dashed line indicating the performance threshold
    geom_hline(yintercept = x$threshold,
               linetype = "dashed",
               color = "black",
               linewidth = 0.6) +
    # Add a text label for the threshold line
    annotate("text",
             x = label_x_pos,
             y = x$threshold + 0.004, # Position slightly above the line
             label = "Threshold (min +5%)",
             hjust = 0,
             size = 3.5,
             color = "black") +
    # Labels and theme
    labs(x = param_expr, y = "Validation MAE") +
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

  # Set y-axis limits if a factor is provided
  if (!is.null(y_limit_factor)) {
    if (!is.numeric(y_limit_factor) || y_limit_factor <= 1) {
      warning("`y_limit_factor` must be a number > 1. Using automatic scaling.")
    } else {
      upper_limit <- x$threshold * y_limit_factor
      lower_limit <- min(plot_data$min_mae, na.rm = TRUE) * 0.99
      p <- p + scale_y_continuous(
        limits = c(lower_limit, upper_limit),
        labels = scales::label_number(big.mark = "")
      )
    }
  } else {
    # Default y-axis formatting
    p <- p + scale_y_continuous(labels = scales::comma)
  }

  # Save the plot if requested
  if (save_plot) {
    if (missing(output_dir)) {
      stop("'output_dir' must be provided when save_plot is TRUE.", call. = FALSE)
    }
    filename <- file.path(output_dir,
                          paste0("parameter_sensitivity_", x$param_name, ".pdf"))

    # Save with a white background using the standard ggsave function
    ggplot2::ggsave(filename, p, width = width, height = height,
                    device = "pdf", units = "in", bg = "white")
  }

  return(p)
}

#' Print Method for Parameter Sensitivity Objects
#'
#' @description
#' The S3 print method for `parameter_sensitivity` objects. It displays a concise
#' summary of the analysis results, including the parameter analyzed, the minimum
#' error found, and the performance threshold.
#'
#' @param x A `parameter_sensitivity` object.
#' @param ... Additional arguments passed to the print function (not currently used).
#'
#' @return Invisibly returns the original object. Called for its side effect of
#'   printing a summary to the console.
#'
#' @method print parameter_sensitivity
#' @export
print.parameter_sensitivity <- function(x, ...) {
  cat("--- Parameter Sensitivity Analysis ---\n")
  cat("Parameter Analyzed:", x$param_name, "\n")
  cat("Number of Bins:", length(x$param_values), "\n")
  cat("------------------------------------\n")
  cat("Minimum MAE Found:", round(x$min_value, 6), "\n")
  cat("Performance Threshold (min +5%):", round(x$threshold, 6), "\n")
  cat("------------------------------------\n")
  cat("Sample Counts per Bin (Min/Median/Max):",
      min(x$sample_counts), "/",
      median(x$sample_counts), "/",
      max(x$sample_counts), "\n")
  invisible(x)
}


#' Calculate Weighted Marginal Distributions
#'
#' @description
#' Calculates the marginal probability distribution for each model parameter. The
#' distributions are weighted by the likelihood of each sample, making this useful for
#' identifying the most probable parameter values from a set of Monte Carlo samples.
#'
#' @param samples A data frame containing parameter samples (e.g., `log_N`, `log_k0`)
#'   and a negative log-likelihood column named `NLL`.
#'
#' @return A named list where each element is a density object (a list with `x`
#'   and `y` components) corresponding to a model parameter.
#'   \item{x}{Vector of parameter values}
#'   \item{y}{Vector of density estimates}
#'
#' @details
#' This function uses the `weighted_kde` helper to perform kernel density
#' estimation for each parameter, with weights derived from the normalized
#' likelihoods of the samples.
#'
#' @importFrom stats na.omit
#' @export
calculate_weighted_marginals <- function(samples) {
  # --- Input Validation ---
  required_cols <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion", "NLL")
  if (!all(required_cols %in% names(samples))) {
    stop("Missing required columns: ",
         paste(setdiff(required_cols, names(samples)), collapse = ", "))
  }

  # Ensure required columns are numeric
  if (!all(sapply(samples[required_cols], is.numeric))) {
    stop("All required parameter and NLL columns must be numeric.")
  }
  # Validate NLL values
  if (all(is.infinite(samples$NLL))) {
    stop("All NLL values are infinite")
  }
  if (any(is.na(samples$NLL))) {
    warning("NA values in the NLL column will be removed.")
    samples <- samples[!is.na(samples$NLL), ]
  }

  if (nrow(samples) < 2) {
    stop("At least two valid samples are required after filtering.")
  }

  # --- Data Cleaning and Weight Calculation ---
  # Remove outliers to stabilize density estimation
  samples <- as.data.frame(lapply(samples, clean_data, k = 3))
  samples <- na.omit(samples)

  # Calculate importance weights from log-likelihoods
  log_likelihoods <- -samples$NLL
  # Shift to be positive and normalize
  std_log_likelihoods <- log_likelihoods - min(log_likelihoods, na.rm = TRUE) + 0.05
  weights <- std_log_likelihoods / sum(std_log_likelihoods, na.rm = TRUE)

  # Define the parameter columns to process
  vars <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")

  # Calculate the marginal density for each parameter by applying weighted_kde.
  # A simple loop is used here because the `weighted_kde` function itself handles
  # parallelism internally, which avoids nested parallel calls.
  marginals <- lapply(vars, function(var) {
    weighted_kde(samples[[var]], weights = weights)
  })

  names(marginals) <- vars
  return(marginals)
}


#' Evaluate a Parameter Set with Cross-Validation
#'
#' @description
#' This internal function calculates the cross-validated likelihood for a given
#' set of parameters. It splits the data into training and validation sets across
#' multiple folds, fits the topolow model on each training set, and evaluates the
#' error on the corresponding validation set.
#'
#' @details
#' To calculate a single Negative Log-Likelihood (NLL) value per parameter set,
#' the function uses a "pooled errors" approach. It combines all out-of-sample
#' errors from every fold into a single set before calculating the NLL and the
#' overall Mean Absolute Error (MAE). This method respects the underlying error
#' distribution and correctly accounts for the total number of validation points.
#'
#' @param dissimilarity_matrix The input dissimilarity matrix to fit.
#' @param mapping_max_iter The maximum number of optimization iterations.
#' @param relative_epsilon The convergence threshold for optimization.
#' @param N The number of dimensions for the embedding.
#' @param k0 The initial spring constant.
#' @param cooling_rate The spring constant decay rate.
#' @param c_repulsion The repulsion constant.
#' @param folds The number of cross-validation folds.
#' @param num_cores The number of cores for parallel processing.
#'
#' @return A list containing the pooled `Holdout_MAE` and the `NLL`.
#'
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply mclapply
#' @keywords internal
likelihood_function <- function(dissimilarity_matrix, mapping_max_iter,
                                relative_epsilon, N, k0, cooling_rate,
                                c_repulsion, folds = 20, num_cores = 1) {
  # --- Manual k-Fold Creation ---
  # This approach creates k folds by iteratively removing a random subset of
  # non-NA values from the full matrix for each fold's training set.
  truth_matrix <- dissimilarity_matrix

  # Create a list to hold the [truth, training] matrix pairs for each fold
  matrix_list <- replicate(folds, list(truth_matrix, NULL), simplify = FALSE)

  num_elements <- sum(!is.na(dissimilarity_matrix))
  # Size of the holdout set for each fold. /2 because each entry (i,j) has a symmetric counterpart (j,i).
  holdout_size <- floor(num_elements / (folds * 2))

  # Create a temporary copy of the matrix to draw samples from without replacement
  D_train <- dissimilarity_matrix

  for(i in 1:folds) {
    # It's possible that D_train becomes too sparse to sample from
    if (sum(!is.na(D_train)) < holdout_size) {
        warning("Could not create all folds due to data sparsity. Using fewer folds.")
        folds <- i - 1
        matrix_list <- matrix_list[1:folds]
        break
    }

    random_indices <- sample(which(!is.na(D_train)), size = holdout_size)
    input_matrix <- dissimilarity_matrix # Start with the full matrix

    # Create the training matrix for the current fold by setting holdout values to NA
    for(index in random_indices) {
      row <- (index - 1) %/% nrow(dissimilarity_matrix) + 1
      col <- (index - 1) %% ncol(dissimilarity_matrix) + 1
      input_matrix[row, col] <- NA
      input_matrix[col, row] <- NA
      # Remove these indices from the pool for the next iteration
      D_train[row, col] <- NA
      D_train[col, row] <- NA
    }
    matrix_list[[i]][[2]] <- input_matrix
  }

  # --- Process Each Fold ---
  # Define the function to run topolow on a single fold
  process_sample <- function(i) {
    truth_matrix <- matrix_list[[i]][[1]]
    input_matrix <- matrix_list[[i]][[2]]

    tryCatch({
      # Call euclidean_embedding with the current parameters
      res_train <- euclidean_embedding(
        dissimilarity_matrix = input_matrix,
        ndim = N,
        mapping_max_iter = mapping_max_iter,
        k0 = k0,
        cooling_rate = cooling_rate,
        c_repulsion = c_repulsion,
        relative_epsilon = relative_epsilon,
        convergence_counter = 5,
        verbose = FALSE
      )

      # Use the topolow return object 'est_distances'
      p_dist_mat <- as.matrix(res_train$est_distances)

      # Calculate errors on the holdout set
      errors <- error_calculator_comparison(p_dist_mat, truth_matrix, input_matrix)
      df <- errors$report_df

      out_sample_errors <- df$OutSampleError[!is.na(df$OutSampleError)]
      n_samples <- length(out_sample_errors)
      sum_abs_errors <- sum(abs(out_sample_errors))

      # Calculate fold-specific MAE for reference
      mae_holdout <- if(n_samples > 0) sum_abs_errors / n_samples else NA

      # Return results needed for pooling
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

  # --- Parallel Execution ---
  if(num_cores > 1 && folds > 1) {
    if(.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(min(num_cores, folds))
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Export required variables to the cluster
      parallel::clusterExport(cl,
                              c("matrix_list", "N", "k0", "cooling_rate", "c_repulsion",
                                "mapping_max_iter", "relative_epsilon", "euclidean_embedding",
                                "error_calculator_comparison"),
                              envir = environment())

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

  pooled_mae <- if(total_samples > 0) total_abs_errors / total_samples else NA
  # NLL is derived from the properties of the Laplace distribution, where MAE is the MLE for the scale parameter b.
  pooled_nll <- if(!is.na(pooled_mae)) total_samples * (1 + log(2 * pooled_mae)) else NA

  return(list(Holdout_MAE = pooled_mae, NLL = pooled_nll))
}


#' Plot Method for profile_likelihood Objects
#'
#' @description
#' Creates a visualization of the profile likelihood for a parameter, showing the
#' maximum likelihood estimates and the 95% confidence interval. It supports
#' mathematical notation for parameter names for clearer plot labels.
#'
#' @details
#' The 95% confidence interval is determined using the likelihood ratio test, where the
#' cutoff is based on the chi-squared distribution:
#' \eqn{LR(\theta_{ij}) = -2[log L_{max}(\theta_{ij}) - log L_{max}(\hat{\theta})]}.
#' The interval includes all parameter values \eqn{\theta_{ij}} for which
#' \eqn{LR(\theta_{ij}) \leq \chi^2_{1,0.05} \approx 3.84}.
#'
#' @param x A `profile_likelihood` object returned by `profile_likelihood()`.
#' @param LL_max The global maximum log-likelihood value from the entire sample set,
#'   used as the reference for calculating the confidence interval.
#' @param width,height Numeric. The width and height of the output plot in inches.
#' @param save_plot Logical. If TRUE, the plot is saved to a file.
#' @param output_dir Character. The directory where the plot will be saved. Required if `save_plot` is TRUE.
#' @param ... Additional arguments passed to the plot function.
#'
#' @return A ggplot object representing the profile likelihood plot.
#'
#' @examples
#' \donttest{
#' # This example can take more than 5 seconds to run.
#' # Create a sample data frame of MCMC samples
#' samples <- data.frame(
#'   log_N = log(runif(50, 2, 10)),
#'   log_k0 = log(runif(50, 1, 5)),
#'   log_cooling_rate = log(runif(50, 0.01, 0.1)),
#'   log_c_repulsion = log(runif(50, 0.1, 1)),
#'   NLL = runif(50, 20, 100)
#' )
#'
#' # Calculate profile likelihood for the "log_N" parameter
#' pl_result <- profile_likelihood("log_N", samples, grid_size = 10)
#'
#' # Provide the global maximum log-likelihood from the samples
#' LL_max <- max(-samples$NLL)
#'
#' # The plot function requires the ggplot2 package
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(pl_result, LL_max, width = 4, height = 3)
#' }
#' }
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_text labs theme_minimal theme element_text element_line element_blank element_rect margin scale_y_continuous ggsave
#' @method plot profile_likelihood
#' @export
plot.profile_likelihood <- function(x, LL_max, width = 3.5, height = 3.5,
                                    save_plot = FALSE, output_dir, ...) {
                                      # check if required packages are loaded
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("The 'scales' package is required for formatting numbers. Please install it first.",
         call. = FALSE)
  }
  # Convert the profile_likelihood object to a data frame for plotting
  plot_df <- data.frame(
    param = x$param,
    LL = x$ll
  )

  # Use mathematical expressions for common topolow parameter names for better plot labels
  param_expr <- switch(x$param_name,
                       "log_N" = expression(log(N)),
                       "log_c_repulsion" = expression(log(c[repulsion])),
                       "log_cooling_rate" = expression(log(alpha)),
                       "log_k0" = expression(log(k[0])),
                       x$param_name) # Default to the original name if not matched

  # Calculate the log-likelihood value for the 95% confidence interval cutoff
  CI_95_LL <- LL_max - (qchisq(0.95, df=1) / 2) # More precisely ~1.92

  # Create the plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$param, y = .data$LL)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = CI_95_LL,
               linetype = "dashed", color = "black", linewidth = 0.5) +
    ggplot2::annotate("text", x = min(plot_df$param), y = CI_95_LL,
                      label = "95% CI", color = "black", vjust = -0.5, hjust = -0.05, size = 3) +
    ggplot2::labs(title = paste("Profile Likelihood:", x$param_name),
         x = param_expr,
         y = "Log Likelihood") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 7, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_text(size = 6),
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ""))

  if (save_plot) {
    if (missing(output_dir)) {
      stop("'output_dir' must be provided when save_plot is TRUE.", call. = FALSE)
    }
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    filename <- file.path(output_dir, paste0("profile_likelihood_", x$param_name, ".pdf"))

    # Use the ggsave_white_bg wrapper for saving
    ggsave_white_bg(filename, p, width = width, height = height, device = "pdf", units = "in")
  }

  return(p)
}


#' Print Method for profile_likelihood Objects
#'
#' @description
#' Provides a concise summary of a `profile_likelihood` object.
#'
#' @param x A `profile_likelihood` object.
#' @param ... Additional arguments passed to `print`.
#' @return The original `profile_likelihood` object (invisibly). Called for its
#'   side effect of printing a summary to the console.
#' @method print profile_likelihood
#' @export
print.profile_likelihood <- function(x, ...) {
  cat("--- Profile Likelihood Analysis ---\n")
  cat("Parameter:", x$param_name, "\n")
  cat("Grid Points Evaluated:", length(x$param), "\n")
  cat("Bandwidth Used:", format(x$bandwidth, digits = 6), "\n")
  cat("Sample Counts per Window (Min/Median/Max):",
      min(x$sample_counts), "/",
      median(x$sample_counts), "/",
      max(x$sample_counts), "\n")
  invisible(x)
}
