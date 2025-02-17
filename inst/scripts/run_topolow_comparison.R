# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# inst/scripts/run_topolow_comparison.R

library(topolow)
library(tidyverse)
library(data.table)
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

matrix_list <- readRDS(matrix_list_file_name)
truth_matrix <- matrix_list[[i]][[1]]
input_matrix <- matrix_list[[i]][[2]]

result <- topolow_full(distance_matrix=input_matrix, 
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