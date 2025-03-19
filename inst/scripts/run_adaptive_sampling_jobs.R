# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# Licensed under MIT License

# inst/scripts/run_adaptive_sampling_jobs.R

# Check and install required packages if needed
source(system.file("scripts", "check_dependencies.R", package = "topolow"))

# Load required packages
library(topolow)
library(data.table)
library(dplyr)
library(parallel)
library(reshape2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

samples_file <- as.character(args[1])
distance_matrix_file <- as.character(args[2])
max_iter <- as.integer(args[3])
relative_epsilon <- as.numeric(args[4])
num_cores <- as.numeric(args[5])
scenario_name <- as.character(args[6])
i <- as.numeric(args[7])
n_iter <- as.numeric(args[8])
batch_size <- as.numeric(args[9])
output_dir <- as.character(args[10])  
folds <- as.integer(args[11]) 

# Add debug output
cat("Loading samples from:", samples_file, "\n")
cat("Loading distance matrix from:", distance_matrix_file, "\n")
cat("Job index:", i, "\n")

# Load input files and check them
if (!file.exists(distance_matrix_file)) {
  stop("Distance matrix file not found: ", distance_matrix_file)
}
if (!file.exists(samples_file)) {
  stop("Samples file not found: ", samples_file)
}

distance_matrix <- readRDS(distance_matrix_file)
current_samples <- read.csv(samples_file)

# Validate loaded data
if (!is.matrix(distance_matrix)) {
  stop("Invalid distance matrix format")
}

required_cols <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion", "NLL")
if (!all(required_cols %in% names(current_samples))) {
  stop("Missing required columns in samples file: ",
       paste(setdiff(required_cols, names(current_samples)), collapse=", "))
}

# Run adaptive sampling
adaptive_MC_sampling(
  samples_file = samples_file,
  distance_matrix = distance_matrix,
  n_iter = n_iter,
  batch_size = batch_size,
  max_iter = max_iter,
  relative_epsilon = relative_epsilon,
  folds = folds,
  num_cores = num_cores,
  scenario_name = scenario_name,
  output_dir = output_dir, 
  verbose = TRUE
)