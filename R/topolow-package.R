# Copyright (c) 2024 Omid Arhami o.arhami@gmail.com
# License BSD 3-clause https://github.com/omid-arhami/topolow/blob/main/LICENSE
# R/topolow-package.R
#' @title Robust Euclidean Embedding of Dissimilarity Data
#'
#' @description
#' The `topolow` package provides a robust implementation of the Topolow algorithm. It 
#' is designed to embed objects into a low-dimensional Euclidean space from a matrix of
#' pairwise dissimilarities, even when the data do not satisfy metric or Euclidean 
#' axioms. The package is particularly well-suited for sparse or incomplete datasets 
#' and includes methods for handling censored (thresholded) data. The package provides 
#' tools for processing antigenic assay data, and visualizing antigenic maps.
#'
#'
#' @details
#' The core of the package is a physics-inspired, gradient-free optimization framework.
#' It models objects as particles in a physical system, where observed dissimilarities
#' define spring rest lengths and unobserved pairs exert repulsive forces. Key features include:
#' \itemize{
#'   \item Quantitative reconstruction of metric space from non-metric data.
#'   \item Robustness against local optima, especially for sparse data, due to a
#'     stochastic pairwise optimization scheme.
#'   \item A statistically grounded approach based on maximizing the likelihood under a
#'     Laplace error model.
#'   \item Tools for parameter optimization, cross-validation, and convergence diagnostics.
#'   \item Support for parallel processing
#'   \item Cross-validation and error analysis
#'   \item A comprehensive suite of visualization functions for network analysis and results.
#'   \item Processing and visualization of antigenic maps
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{Euclidify}}: Wizard function to run all steps of the Topolow algorithm automatically
#'   \item \code{\link{euclidean_embedding}}: Core embedding algorithm
#'   \item \code{\link{initial_parameter_optimization}}: Find optimal parameters using Latin Hypercube Sampling.
#'   \item \code{\link{run_adaptive_sampling}}: Refine parameter estimates with adaptive Monte Carlo sampling.
#' }
#'
#' @section Output Files:
#' Functions that generate output files (like parameter optimization results) will 
#' create subdirectories in a user-specified directory (via output_dir parameter)
#'
#' The following subdirectories may be created:
#' \itemize{
#'   \item model_parameters/: Contains optimization results and parameter evaluations
#'   \item init_param_optimization/: Contains files and outputs when using initial_parameter_optimization
#' }
#'
#' @section Citation:
#' If you use this package, please cite the Bioinformatics paper:
#' Omid Arhami, Pejman Rohani, 
#' Topolow: A mapping algorithm for antigenic cross-reactivity and binding affinity assays, 
#' Bioinformatics, 2025;, btaf372,
#' https://doi.org/10.1093/bioinformatics/btaf372
#' \doi{10.1093/bioinformatics/btaf372}.
#' 
#' `bibtex` entry:
#'  title=\{Topolow: a mapping algorithm for antigenic cross-reactivity and binding affinity assays\},
#'  author=\{Arhami, Omid and Rohani, Pejman\},
#'  journal=\{Bioinformatics\},
#'  volume=\{41\},
#'  number=\{7\},
#'  pages=\{btaf372\},
#'  year=\{2025\},
#'  issn = \{1367-4811\},
#'  doi = \{10.1093/bioinformatics/btaf372\},
#'  url = \{https://doi.org/10.1093/bioinformatics/btaf372\},
#'  eprint = \{https://academic.oup.com/bioinformatics/article-pdf/41/7/btaf372/63582086/btaf372.pdf\},
#'  publisher=\{Oxford University Press\}
#'
#' And/or the preprint on mathematical properties:
#' Omid Arhami, Pejman Rohani, 
#' Topolow: Force-Directed Euclidean Embedding of Dissimilarity Data with Robustness Against Non-Metricity and Sparsity, 
#' arXiv:2508.01733,
#' https://doi.org/10.48550/arXiv.2508.01733
#' \doi{10.48550/arXiv.2508.01733}.
#' 
#' `bibtex` entry:
#'  title=\{Topolow: Force-Directed Euclidean Embedding of Dissimilarity Data with Robustness Against Non-Metricity and Sparsity\},
#'  author=\{Arhami, Omid and Rohani, Pejman\},
#'  year=\{2025\},
#'  doi = \{10.48550/arXiv.2508.01733\},
#'  url = \{https://arxiv.org/abs/2508.01733\},
#'  publisher=\{arXiv\}
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
#' @importFrom stats qchisq setNames
#' @importFrom rlang .data
#' @importFrom future availableCores
#' @importFrom reshape2 melt
#' @importFrom dplyr %>% select filter mutate group_by summarise ungroup
#' @importFrom filelock lock unlock
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal coord_fixed ggsave theme element_text margin geom_line geom_density labs element_blank element_rect
#' @importFrom lhs maximinLHS
#' @importFrom parallel mclapply detectCores makeCluster clusterExport clusterEvalQ parLapply stopCluster
#' @importFrom reshape2 melt
#' @importFrom rlang sym
#' @importFrom stats dist na.omit sd qunif complete.cases aggregate bw.nrd0 approx dnorm cov hclust lm coef qt
#' @importFrom utils read.csv write.csv write.table tail
#' @importFrom data.table setDT
## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL