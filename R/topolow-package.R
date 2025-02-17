# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# R/topolow-package.R

#' @title Topolow: A mapping algorithm for antigenic cross-reactivity and binding affinity assay results
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
#'   \item model_parameters/: Contains optimization results and parameter evaluations
#'   \item init_param_optimization/: Contains SLURM job files and outputs when using SLURM
#' }
#'
#' @section Citation:
#' If you use this package, please cite:
#' Arhami O, Rohani P (2025).
#' Topolow: A mapping algorithm for antigenic cross-reactivity and binding affinity assay results.
#' \emph{https://doi.org/10.1101/2025.02.09.637307}.
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