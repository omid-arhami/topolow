# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# R/algorithm_comparison_helpers.R

#' Algorithm Comparison Helper Functions
#' 
#' @description
#' Helper functions for running RACMACS and Topolow during algorithm comparisons.
#' Handles both local and SLURM execution.
#'
#' @keywords internal
"_PACKAGE"


#' Run Topolow via SLURM  (LEGACY)
#'
#' @description
#' Submits Topolow job to SLURM cluster and collects results.
#'
#' @param truth_matrix Matrix of true distances
#' @param input_matrix Matrix for optimization
#' @param best_params List of optimal parameters
#' @param opt_params List of optimization parameters
#' @param scenario_name Character scenario identifier
#' @param fold Integer fold number
#' @param cider Logical; whether to use cider queue
#' @return List of performance metrics
#' @keywords internal
run_topolow_slurm <- function(truth_matrix, input_matrix, best_params,
                              opt_params, scenario_name, fold, cider,
                              time = "8:00:00", memory = "18G",
                              matrix_dir, results_dir) {
  
  # Save matrices for this fold
  matrix_list <- list(list(truth_matrix, input_matrix))
  matrix_file <- file.path(matrix_dir,
                           sprintf("%s_fold%d_matrices.rds", 
                                   scenario_name, fold))
  saveRDS(matrix_list, matrix_file)
  
  # Create argument list with absolute paths
  args_list <- c(
    normalizePath(matrix_file),
    as.character(best_params$N),
    as.character(opt_params$max_topolow_iter),
    as.character(best_params$k0),
    as.character(best_params$cooling_rate),
    as.character(best_params$c_repulsion),
    scenario_name,
    as.character(fold),
    normalizePath(results_dir)  # Pass results directory to SLURM job
  )
  
  # Create and submit job
  job_name <- sprintf("%d_topolow_%s", fold, scenario_name)
  slurm_script <- create_slurm_script(
    job_name = job_name,
    script_path = system.file("scripts/run_topolow_comparison.R",
                              package = "topolow"),
    args = args_list,
    num_cores = 1,
    time = time,
    memory = memory,
    output_file = file.path(results_dir,
                            sprintf("%s_fold%d.out", scenario_name, fold)),
    error_file = file.path(results_dir,
                           sprintf("%s_fold%d.err", scenario_name, fold))
  )
  
  submit_job(slurm_script, cider = cider)
  
  # Wait for and collect results
  result_file <- file.path(results_dir,
                           sprintf("%s_fold%d_results.csv", scenario_name, fold))
  
  # Wait for result file with timeout
  wait_time <- 0
  max_wait <- 360000  # 100 hour timeout
  while (!file.exists(result_file) && wait_time < max_wait) {
    Sys.sleep(350) # Check every 5 minutes
    wait_time <- wait_time + 350
  }
  
  if (!file.exists(result_file)) {
    stop("Timeout waiting for Topolow SLURM job results")
  }
  
  # Read and return results
  results <- read.csv(result_file)
  list(
    mae = results$mae,
    coverage = results$coverage,
    correlation = results$correlation
  )
}


#' Run Topolow Locally
#'
#' @description
#' Runs Topolow optimization locally for algorithm comparison.
#'
#' @param truth_matrix Matrix of true distances
#' @param input_matrix Matrix for optimization
#' @param best_params List of optimal parameters
#' @param opt_params List of optimization parameters
#' @return List of performance metrics
#' @keywords internal
run_topolow_local <- function(truth_matrix, input_matrix, best_params, opt_params) {
  # Run Topolow optimization
  result <- topolow_full(
    distance_matrix = input_matrix,
    ndim = best_params$N,
    max_iter = opt_params$max_topolow_iter,
    k0 = best_params$k0,
    cooling_rate = best_params$cooling_rate,
    c_repulsion = best_params$c_repulsion,
    write_positions_to_csv = FALSE,
    verbose = FALSE
  )
  
  # Calculate performance metrics
  topolow_errors <- error_calculator_comparison(
    p_dist_mat = result$est_distances,
    truth_matrix = truth_matrix,
    input_matrix = input_matrix
  )
  
  list(
    errors = errors$report_df,
    coverage = errors$coverage,
    correlation = errors$OutSampleCor,
    positions = result$positions
  )
}


#' Run RACMACS Optimization
#'
#' @description
#' Runs RACMACS optimization for algorithm comparison.
#'
#' @param truth_matrix Matrix of true distances
#' @param input_matrix Matrix for optimization
#' @param opt_params List of optimization parameters
#' @return List of performance metrics
#' @importFrom Racmacs acmap optimizeMap keepBestOptimization
#' @keywords internal
run_racmacs <- function(truth_matrix, input_matrix, opt_params) {
  # Convert to titer format
  titer_table <- dist_to_titer_table(input_matrix, base = 2, tens = 10)
  
  # Run RACMACS optimization
  map <- acmap(titer_table = titer_table)
  
  map <- optimizeMap(
    map = map,
    number_of_dimensions = 2,  # RACMACS uses 2D
    number_of_optimizations = opt_params$racmacs_opt_rounds,
    minimum_column_basis = "none"
  )
  
  map <- keepBestOptimization(map)
  
  # Extract coordinates
  coords <- rbind(
    cbind(name = map$antigens$name,
          x = map$antigens$coords[,1],
          y = map$antigens$coords[,2],
          type = "antigen"),
    cbind(name = map$sera$name,
          x = map$sera$coords[,1],
          y = map$sera$coords[,2],
          type = "serum")
  )
  
  coords <- as.data.frame(coords)
  coords$x <- as.numeric(coords$x)
  coords$y <- as.numeric(coords$y)
  
  # Calculate distances
  p_dist_mat <- as.matrix(dist(coords[,c("x", "y")]))
  rownames(p_dist_mat) <- colnames(p_dist_mat) <- coords$name
  
  # Calculate performance metrics
  racmacs_errors <- error_calculator_comparison(
    p_dist_mat = p_dist_mat,
    truth_matrix = truth_matrix,
    input_matrix = input_matrix
  )
  
  list(
    mae = mean(abs(racmacs_errors$report_df$OutSampleError), na.rm = TRUE),
    coverage = racmacs_errors$coverage,
    correlation = racmacs_errors$OutSampleCor
  )
}
