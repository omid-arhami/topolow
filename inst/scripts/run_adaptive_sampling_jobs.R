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

# --- START OF FIX ---
initial_samples_file <- as.character(args[1]) # This is the MASTER file
distance_matrix_file <- as.character(args[2])
mapping_max_iter <- as.integer(args[3])
relative_epsilon <- as.numeric(args[4])
num_cores <- as.integer(args[5])
scenario_name <- as.character(args[6])
job_id <- as.integer(args[7]) # This is the SLURM_ARRAY_TASK_ID
iterations <- as.integer(args[8])
output_dir <- as.character(args[9])
folds <- as.integer(args[10])

cat("Running job ID:", job_id, "for scenario:", scenario_name, "\n")

# Load the distance matrix
if (!file.exists(distance_matrix_file)) stop("Distance matrix file not found: ", distance_matrix_file)
distance_matrix <- readRDS(distance_matrix_file)
if (!is.matrix(distance_matrix)) stop("Invalid distance matrix format")

# --- THIS IS THE CRITICAL NEW LOGIC ---
# Each job must create its OWN temporary file to work on.
# This prevents race conditions and allows the 'gather' script to find the results.
adaptive_dir <- file.path(output_dir, "adaptive_sampling_jobs")

# Construct the unique temporary file path for this specific job
job_output_file <- file.path(adaptive_dir, sprintf("job_%02d_%s.csv", job_id, scenario_name))

# Copy the initial samples from the master file to this job's unique file
if (!file.exists(initial_samples_file)) stop("Initial samples file not found: ", initial_samples_file)
file.copy(initial_samples_file, job_output_file, overwrite = TRUE)

cat("Job", job_id, "will write results to:", job_output_file, "\n")
# --- END OF FIX ---

# Run adaptive sampling on the NEWLY CREATED job-specific file
adaptive_MC_sampling(
  samples_file       = job_output_file, # Use the unique file for this job
  distance_matrix    = distance_matrix,
  iterations         = iterations,
  batch_size         = 1,
  mapping_max_iter   = mapping_max_iter,
  relative_epsilon   = relative_epsilon,
  folds              = folds,
  num_cores          = num_cores,
  scenario_name      = scenario_name,
  output_dir         = output_dir,
  verbose            = TRUE # Set to TRUE for better debugging on SLURM nodes
)

cat("Job", job_id, "has completed.\n")