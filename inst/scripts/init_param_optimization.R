# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# inst/scripts/init_param_optimization.R

# Load required packages
library(topolow)
library(data.table)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

matrix_list_file <- args[1]
N <- as.integer(args[2])
max_iter <- as.integer(args[3])
k0 <- as.numeric(args[4])
cooling_rate <- as.numeric(args[5])
c_repulsion <- as.numeric(args[6])
relative_epsilon <- as.numeric(args[7])
convergence_counter <- as.integer(args[8])
initial_positions <- NULL # args[9] is "NULL"
write_positions <- FALSE # args[10] is "FALSE"
verbose <- FALSE # args[11] is "FALSE"
scenario_name <- args[12]
i <- as.numeric(args[13])
output_dir <- args[14]
num_samples <- as.integer(args[15])

# Create output directories if they don't exist
param_dir <- file.path(output_dir, "init_param_optimization")
if (!dir.exists(param_dir)) {
  dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)
}

# Add debug output
cat("Loading matrix_list from:", matrix_list_file, "\n")
cat("i value:", i, "\n")

# Load matrix list and check its length
matrix_list <- readRDS(matrix_list_file)
cat("Length of matrix_list:", length(matrix_list), "\n")

# Calculate fold index
data_idx <- floor((i-1) / num_samples) + 1

# Validate data_idx before using it
if (data_idx > length(matrix_list)) {
  stop(sprintf("Data fold index i=%d. matrix_list has length %d", 
               data_idx, length(matrix_list)))
}

# Get matrices for this sample-fold combination
truth_matrix <- matrix_list[[data_idx]][[1]]
input_matrix <- matrix_list[[data_idx]][[2]]

# Run optimization
res_train <- topolow_full(
  distance_matrix = input_matrix,
  ndim = N,
  max_iter = max_iter,
  k0 = k0,
  cooling_rate = cooling_rate,
  c_repulsion = c_repulsion,
  relative_epsilon = relative_epsilon,
  convergence_counter = convergence_counter,
  initial_positions = initial_positions,
  write_positions_to_csv = write_positions,
  verbose = verbose
)

p_dist_mat <- res_train$est_distances

# Calculate errors
errors <- error_calculator_comparison(
  p_dist_mat = p_dist_mat,
  truth_matrix = truth_matrix,
  input_matrix = input_matrix
)

df <- errors$report_df
mae_holdout <- mean(abs(df$OutSampleError), na.rm = TRUE)

# Calculate negative log likelihood
n <- sum(!is.na(df$OutSampleError))
NLL <- n * (1 + log(2) + log(mae_holdout))

# Create results data frame
result <- data.frame(
  N = N,
  k0 = k0,
  cooling_rate = cooling_rate,
  c_repulsion = c_repulsion,
  Holdout_MAE = mae_holdout,
  NLL = NLL
)

# Save results
output_file <- file.path(param_dir,
                         paste0(i, "_params_", scenario_name, ".csv"))
write.csv(result, output_file, row.names = FALSE)