# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# Licensed under Pre-Publication Academic License https://github.com/omid-arhami/topolow/blob/main/LICENSE

# inst/scripts/run_adaptive_sampling_jobs.R
# Script for running adaptive Monte Carlo sampling jobs submitted via SLURM

# Check and install required packages if needed
source(system.file("scripts", "check_dependencies.R", package = "topolow"))

# Load required packages
library(topolow)
library(data.table)
library(dplyr)
library(parallel)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

samples_file <- as.character(args[1])
distance_matrix_file <- as.character(args[2])
mapping_max_iter <- as.integer(args[3])
relative_epsilon <- as.numeric(args[4])
num_cores <- as.numeric(args[5]) # 1 but Still needed for internal implementation
scenario_name <- as.character(args[6])
i <- as.numeric(args[7])
iterations <- as.numeric(args[8]) # = num_samples , always = 1 in jobs submitted to slurm
output_dir <- as.character(args[9])  
folds <- as.integer(args[10]) 

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

# Run adaptive sampling - using original function call with batch_size=1
adaptive_MC_sampling(
  samples_file = samples_file,
  distance_matrix = distance_matrix,
  iterations = iterations, # Always set to 1
  batch_size = 1, # Always set to 1
  mapping_max_iter = mapping_max_iter,
  relative_epsilon = relative_epsilon,
  folds = folds,
  num_cores = num_cores, # Use 1 core per job (from script arg)
  scenario_name = scenario_name,
  output_dir = output_dir, 
  verbose = FALSE
)