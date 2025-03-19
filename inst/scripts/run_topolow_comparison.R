# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# Licensed under MIT License

# inst/scripts/run_topolow_comparison.R

# Check and install required packages if needed
source(system.file("scripts", "check_dependencies.R", package = "topolow"))

# Load required packages
library(topolow)
library(tidyverse)
library(data.table)
library(reshape2)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

matrix_list_file_name = as.character(args[1])
N = as.integer(args[2])
max_iter = as.integer(args[3])
k0 = as.numeric(args[4])
cooling_rate = as.numeric(args[5])
c_repulsion =  as.numeric(args[6])
scenario_name = as.character(args[7])
i = as.numeric(args[8])
results_dir = as.character(args[9])

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}

matrix_list <- readRDS(matrix_list_file_name)
truth_matrix <- matrix_list[[i]][[1]]
input_matrix <- matrix_list[[i]][[2]]

result <- create_topolow_map(distance_matrix=input_matrix, 
                       ndim=N, 
                       max_iter=max_iter, 
                       k0=k0 , 
                       cooling_rate=cooling_rate, 
                       c_repulsion=c_repulsion, 
                       relative_epsilon=10^-9,
                       convergence_counter= 5,
                       initial_positions=NULL,
                       write_positions_to_csv=FALSE,
                       verbose = FALSE)

p_dist_mat <- result$est_distances
p_dist_mat <- as.matrix(p_dist_mat)

## calc errors:
topolow_errors <- error_calculator_comparison(p_dist_mat, truth_matrix, input_matrix)
topolow_df <- topolow_errors$report_df
topolow_df$Dimension <- N
topolow_df$Algorithm <- "TopoLow"
topolow_df$Scenario <- scenario_name
topolow_df$Fold <- i

write.csv(topolow_df, file.path(results_dir,
                                sprintf("%s_fold%d_errors.csv", scenario_name, i)), 
          row.names = FALSE)

# count the recovered:
topolow_coverage <- topolow_errors$coverage

mapped_objects_df= data.table(Dimension=N, Algorithm="TopoLow", 
                              Scenario = scenario_name, Fold = i,
                                Mapped= nrow(p_dist_mat),
                                Total= nrow(truth_matrix),
                                Validation_Coverage= topolow_coverage)

write.csv(mapped_objects_df, file.path(results_dir,
                                       sprintf("%s_fold%d_coverage.csv", scenario_name, i)), 
          row.names = FALSE)

# Cor
topolow_cor <- topolow_errors$OutSampleCor

cor_df <- data.table(Dimension=N, Algorithm="TopoLow",
                     Scenario = scenario_name, Fold = i, Validation_Cor= topolow_cor)

write.csv(cor_df, file.path(results_dir,
                            sprintf("%s_fold%d_cor.csv", scenario_name, i)),  
          row.names = FALSE)

# Write combined results file for compatibility
result_df <- data.frame(
  Dimension = N,
  Algorithm = "TopoLow",
  Scenario = scenario_name,
  Fold = i,
  Holdout_MAE = mean(abs(topolow_df$OutSampleError), na.rm = TRUE),
  Validation_Coverage = topolow_coverage,
  Validation_Cor = topolow_cor
)

write.csv(result_df, file.path(results_dir,
                               sprintf("%s_fold%d_results.csv", scenario_name, i)),
          row.names = FALSE)