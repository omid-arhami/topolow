# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# License BSD 3-clause https://github.com/omid-arhami/topolow/blob/main/LICENSE

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
#'   \item Support for parallel processing
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{create_topolow_map}}: Core optimization algorithm
#'   \item \code{\link{process_antigenic_data}}: Process raw antigenic data
#'   \item \code{\link{initial_parameter_optimization}}: Optimize algorithm parameters
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
#'   \item init_param_optimization/: Contains files and outputs when using initial_parameter_optimization
#' }
#'
#' @section Citation:
#' If you use this package, please cite:
#' Omid Arhami, Pejman Rohani, 
#' Topolow: A mapping algorithm for antigenic cross-reactivity and binding affinity assays, 
#' Bioinformatics, 2025;, btaf372,
#' https://doi.org/10.1093/bioinformatics/btaf372
#' \doi{10.1093/bioinformatics/btaf372}.
#'
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # Check for required packages
  required_pkgs <- c("reshape2", "data.table", "dplyr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    warning(
      "The following required packages are missing: ", 
      paste(missing_pkgs, collapse = ", "), 
      "\nPlease install them with: install.packages(c('", 
      paste(missing_pkgs, collapse = "', '"), "'))"
    )
  }
}

#' @importFrom data.table data.table setDT fread
#' @importFrom dplyr %>% select filter mutate group_by summarise ungroup
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom parallel mclapply detectCores
#' @importFrom coda mcmc mcmc.list gelman.diag effectiveSize
#' @importFrom stats median sd var aggregate approx bw.nrd bw.nrd0 bw.ucv bw.bcv bw.SJ coef complete.cases cor cov dist lm na.omit qt residuals setNames mad
#' @importFrom utils read.csv write.csv write.table tail globalVariables
#' @importFrom grDevices colorRampPalette colors
#' @importFrom graphics hist
#' @importFrom vegan protest procrustes
#' @importFrom ape is.rooted unroot dist.nodes nodepath extract.clade
#' @importFrom Racmacs save.coords
NULL