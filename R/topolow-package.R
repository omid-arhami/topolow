# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# Licensed under MIT License
# R/topolow-package.R

#' @title TopoLow: Topological Optimization for Low-Dimensional Mapping
#'
#' @description
#' An implementation of the TopoLow algorithm for antigenic cartography mapping and analysis. 
#' The package provides tools for optimizing point configurations in high-dimensional spaces,
#' handling missing and thresholded measurements, processing antigenic assay data, and 
#' visualizing antigenic maps.
#'
#' @details
#' The package implements a physics-inspired approach combining spring forces and repulsive
#' interactions to find optimal point configurations. Key features include:
#' \itemize{
#'   \item Optimization of point configurations in high-dimensional spaces
#'   \item Handling of missing and thresholded measurements
#'   \item Processing of antigenic assay data
#'   \item Interactive visualization of antigenic maps
#'   \item Cross-validation and error analysis
#'   \item Network structure analysis
#'   \item Support for parallel processing and high-performance computing environments
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{topolow_full}}: Core optimization algorithm
#'   \item \code{\link{topolow_Smith_obj}}: Smith variant for HI assay data
#'   \item \code{\link{process_antigenic_data}}: Process raw antigenic data
#'   \item \code{\link{run_parameter_optimization}}: Optimize algorithm parameters
#'   \item \code{\link{plot_temporal_mapping}}: Create temporal visualizations
#'   \item \code{\link{plot_cluster_mapping}}: Create cluster-based visualizations
#' }
#'
#' @section Output Files:
#' Functions that generate output files (like parameter optimization results) will 
#' create subdirectories in either:
#' \itemize{
#'   \item The current working directory (if output_dir = NULL)
#'   \item A user-specified directory (via output_dir parameter)
#' }
#'
#' The following subdirectories may be created:
#' \itemize{
#'   \item hyperparam_results/: Contains optimization results and hyperparameter evaluations
#'   \item run_topolow_arguments/: Contains SLURM job files and outputs when using SLURM
#' }
#'
#' @section Citation:
#' If you use this package, please cite:
#' [Citation information to be added]
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom data.table data.table setDT fread
#' @importFrom dplyr %>% select filter mutate group_by summarise ungroup
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom stats cor median sd var
#' @importFrom parallel mclapply detectCores
#' @importFrom utils write.csv read.csv
#' @importFrom coda mcmc mcmc.list gelman.diag effectiveSize
NULL