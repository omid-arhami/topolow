# R/euclidify_diagnostics.R
# Copyright (c) 2025 Omid Arhami
# diagnostic and trace plotting functions for Euclidify

#' Plot Euclidify Optimization Diagnostics
#'
#' @description
#' Creates comprehensive diagnostic plots for the Euclidify optimization process,
#' including parameter search trajectory and embedding quality metrics.
#'
#' @param euclidify_result A result object from \code{Euclidify()} with create_diagnostic_plots=TRUE.
#' @param plot_types Character vector specifying which plots to create. Options:
#'   "all", "parameter_search", "convergence", "quality", "cv_errors".
#'   Default: "all"
#' @param save_plots Logical. Whether to save plots to files. Default: TRUE
#' @param output_dir Character. Directory for saving plots. Required if save_plots=TRUE.
#' @param width Numeric. Plot width in inches. Default: 12
#' @param height Numeric. Plot height in inches. Default: 8
#' @param dpi Numeric. Resolution for saved plots. Default: 300
#' @param return_plots Logical. Whether to return plot objects. Default: TRUE
#'
#' @return A list of ggplot objects if return_plots=TRUE, otherwise NULL invisibly.
#'
#' @details
#' This function creates several diagnostic visualizations:
#' \itemize{
#'   \item Parameter Search: Shows explored parameter space in 2D projections
#'   \item Convergence Trace: Plots MAE and parameter evolution during final embedding
#'   \item Quality Metrics: Scatter plots of predicted vs true distances
#'   \item CV Errors: Distribution of cross-validation errors across folds
#' }
#'
#' @examples
#' \dontrun{
#' # Run Euclidify with diagnostic plots
#' result <- Euclidify(
#'   dissimilarity_matrix = my_data,
#'   output_dir = "output",
#'   create_diagnostic_plots = TRUE
#' )
#' 
#' # Plots are automatically saved to output/diagnostics/
#' # View diagnostic report
#' report <- create_diagnostic_report(result)
#' cat(report, sep = "\n")
#' 
#' # Access specific diagnostic plots
#' print(result$diagnostic_plots$parameter_search)
#' print(result$diagnostic_plots$quality)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_density facet_wrap
#'   labs theme_minimal ggsave geom_hline scale_color_viridis_c
#' @importFrom gridExtra grid.arrange
#' @export
plot_euclidify_diagnostics <- function(euclidify_result,
                                       plot_types = "all",
                                       save_plots = TRUE,
                                       output_dir = NULL,
                                       width = 12,
                                       height = 8,
                                       dpi = 300,
                                       return_plots = TRUE) {
  
  # Check if all_samples data is available
  if (is.null(euclidify_result$all_samples)) {
    stop("No parameter samples found in result. Run Euclidify with create_diagnostic_plots=TRUE")
  }
  
  if (save_plots && is.null(output_dir)) {
    stop("output_dir must be specified when save_plots=TRUE")
  }
  
  if (save_plots) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Determine which plots to create
  available_types <- c("parameter_search", "quality", "cv_errors")
  if ("all" %in% plot_types) {
    plot_types <- available_types
  } else {
    plot_types <- intersect(plot_types, available_types)
  }
  
  plots <- list()
  all_samples <- euclidify_result$all_samples
  
  # 1. Parameter Search Trajectory
  if ("parameter_search" %in% plot_types) {
    p1 <- plot_parameter_search(all_samples, euclidify_result$optimal_params)
    plots$parameter_search <- p1
    
    if (save_plots) {
      ggsave(
        file.path(output_dir, "parameter_search_trajectory.png"),
        p1, width = width, height = height, dpi = dpi, bg = "white"
      )
    }
  }
  
  # 2. Embedding Quality
  if ("quality" %in% plot_types && !is.null(euclidify_result$dissimilarity_matrix)) {
    p2 <- plot_embedding_quality(
      euclidify_result$dissimilarity_matrix,
      euclidify_result$est_distances
    )
    plots$quality <- p2
    
    if (save_plots) {
      ggsave(
        file.path(output_dir, "embedding_quality.png"),
        p2, width = width, height = height, dpi = dpi, bg = "white"
      )
    }
  }
  
  # 3. Cross-Validation Errors (if CV data available)
  if ("cv_errors" %in% plot_types && "Holdout_MAE" %in% names(all_samples)) {
    p3 <- plot_cv_errors(all_samples)
    plots$cv_errors <- p3
    
    if (save_plots) {
      ggsave(
        file.path(output_dir, "cv_errors.png"),
        p3, width = width, height = height, dpi = dpi, bg = "white"
      )
    }
  }
  
  if (return_plots) {
    return(plots)
  } else {
    return(invisible(NULL))
  }
}


#' Plot Parameter Search Trajectory
#'
#' @description
#' Visualizes the parameter search process as 2D projections of the parameter space.
#'
#' @param param_history Data frame with parameter evaluations and performance.
#' @param best_params List with optimal parameter values.
#'
#' @return A ggplot object showing parameter space exploration.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_point scale_color_viridis_c
#'   facet_grid labs theme_minimal theme element_text
#' @keywords internal
plot_parameter_search <- function(param_history, best_params) {
  
  # Create all pairwise parameter plots
  params <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  
  # Add iteration number and phase
  param_history$iteration <- seq_len(nrow(param_history))
  
  # Create compound plot
  plots <- list()
  
  # Main scatter: k0 vs cooling_rate (most important interaction)
  p1 <- ggplot(param_history, aes(x = .data$log_k0, y = .data$log_cooling_rate)) +
    geom_point(aes(color = .data$NLL, size = -.data$Holdout_MAE), alpha = 0.6) +
    geom_point(
      data = data.frame(
        log_k0 = log(best_params$k0),
        log_cooling_rate = log(best_params$cooling_rate)
      ),
      color = "red", size = 5, shape = 18
    ) +
    scale_color_viridis_c(option = "plasma", direction = -1) +
    labs(
      title = "Parameter Search Trajectory",
      subtitle = "Red diamond = optimal parameters",
      x = "log(Initial Spring Constant)",
      y = "log(Cooling Rate)",
      color = "Negative Log\nLikelihood",
      size = "Negative MAE"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  # Marginal distributions
  p2 <- ggplot(param_history, aes(x = .data$iteration, y = .data$log_k0)) +
    geom_line(alpha = 0.3) +
    geom_point(aes(color = -.data$NLL), size = 2, alpha = 0.6) +
    geom_hline(yintercept = log(best_params$k0), 
               color = "red", linetype = "dashed") +
    scale_color_viridis_c(option = "plasma") +
    labs(
      y = "log(k0)",
      x = "Iteration",
      color = "-NLL"
    ) +
    theme_minimal()
  
  p3 <- ggplot(param_history, aes(x = .data$iteration, y = .data$log_cooling_rate)) +
    geom_line(alpha = 0.3) +
    geom_point(aes(color = -.data$NLL), size = 2, alpha = 0.6) +
    geom_hline(yintercept = log(best_params$cooling_rate), 
               color = "red", linetype = "dashed") +
    scale_color_viridis_c(option = "plasma") +
    labs(
      y = "log(Cooling Rate)",
      x = "Iteration",
      color = "-NLL"
    ) +
    theme_minimal()
  
  p4 <- ggplot(param_history, aes(x = .data$iteration, y = .data$log_c_repulsion)) +
    geom_line(alpha = 0.3) +
    geom_point(aes(color = -.data$NLL), size = 2, alpha = 0.6) +
    geom_hline(yintercept = log(best_params$c_repulsion), 
               color = "red", linetype = "dashed") +
    scale_color_viridis_c(option = "plasma") +
    labs(
      y = "log(Repulsion)",
      x = "Iteration",
      color = "-NLL"
    ) +
    theme_minimal()
  
  # Combine plots
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
}


#' Plot Cross-Validation Errors
#'
#' @description
#' Visualizes the distribution of CV errors across folds and parameter sets.
#'
#' @param all_samples Data frame with parameter evaluations and CV errors.
#' @return A ggplot object showing CV error distributions.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin labs theme_minimal
#' @keywords internal
plot_cv_errors <- function(all_samples) {
  
  # Plot MAE distribution across all parameter evaluations
  p <- ggplot(all_samples, aes(x = seq_along(.data$Holdout_MAE), y = .data$Holdout_MAE)) +
    geom_line(color = "steelblue", alpha = 0.5) +
    geom_point(aes(color = .data$NLL), size = 2, alpha = 0.7) +
    scale_color_viridis_c(option = "plasma", direction = -1) +
    geom_hline(yintercept = min(all_samples$Holdout_MAE, na.rm = TRUE),
               linetype = "dashed", color = "red") +
    labs(
      title = "Cross-Validation MAE Across Parameter Evaluations",
      subtitle = "Red line = best MAE achieved",
      x = "Parameter Evaluation",
      y = "Cross-Validation MAE",
      color = "NLL"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  return(p)
}


#' Create Summary Diagnostic Report
#'
#' @description
#' Generates a text summary of the Euclidify optimization process.
#'
#' @param euclidify_result A result object from Euclidify with diagnostics.
#' @param output_file Character. Optional path to save report as text file.
#'
#' @return Character vector with report lines (invisibly).
#'
#' @examples
#' \dontrun{
#' result <- Euclidify(..., create_diagnostic_plots = TRUE)
#' report <- create_diagnostic_report(result, "diagnostics/report.txt")
#' cat(report, sep = "\n")
#' }
#'
#' @export
create_diagnostic_report <- function(euclidify_result, output_file = NULL) {
  
  all_samples <- euclidify_result$all_samples
  
  # Calculate some statistics
  n_evaluations <- if (!is.null(all_samples)) nrow(all_samples) else 0
  best_mae <- if (!is.null(all_samples)) min(all_samples$Holdout_MAE, na.rm = TRUE) else NA
  best_nll <- if (!is.null(all_samples)) min(all_samples$NLL, na.rm = TRUE) else NA
  worst_nll <- if (!is.null(all_samples)) max(all_samples$NLL, na.rm = TRUE) else NA
  
  report <- c(
    "====================================",
    "   EUCLIDIFY DIAGNOSTIC REPORT",
    "====================================",
    "",
    "DATA CHARACTERISTICS:",
    sprintf("  Objects: %d", euclidify_result$data_characteristics$n_objects),
    sprintf("  Missing data: %.1f%%", euclidify_result$data_characteristics$missing_percentage),
    sprintf("  Non-Euclidean score: %.4f", euclidify_result$data_characteristics$non_euclidean_score),
    sprintf("  Dimensions for 90%% variance: %d", euclidify_result$data_characteristics$eigenvalue_dims_90),
    "",
    "OPTIMIZATION SUMMARY:",
    sprintf("  Parameter evaluations: %d", n_evaluations),
    sprintf("  Total runtime: %.1f seconds", as.numeric(euclidify_result$runtime)),
    sprintf("  Subsampling used: %s", euclidify_result$optimization_summary$used_subsampling),
    if (!is.na(best_mae)) sprintf("  Best CV MAE achieved: %.4f", best_mae) else "",
    "",
    "OPTIMAL PARAMETERS:",
    sprintf("  Dimensions: %d", euclidify_result$optimal_params$ndim),
    sprintf("  Initial spring constant (k0): %.4f", euclidify_result$optimal_params$k0),
    sprintf("  Cooling rate: %.6f", euclidify_result$optimal_params$cooling_rate),
    sprintf("  Repulsion constant: %.6f", euclidify_result$optimal_params$c_repulsion),
    if (!is.null(euclidify_result$optimal_params$cv_mae)) {
      sprintf("  CV MAE at optimal params: %.4f", euclidify_result$optimal_params$cv_mae)
    } else "",
    "",
    "FINAL EMBEDDING QUALITY:",
    sprintf("  Final MAE: %.4f", euclidify_result$mae),
    "",
    "PARAMETER SEARCH STATISTICS:",
    if (!is.na(best_nll)) sprintf("  Best NLL: %.2f", best_nll) else "",
    if (!is.na(worst_nll)) sprintf("  Worst NLL: %.2f", worst_nll) else "",
    if (!is.na(best_nll) && !is.na(worst_nll)) {
      sprintf("  NLL improvement: %.1f%%", 
              100 * (worst_nll - best_nll) / worst_nll)
    } else "",
    "===================================="
  )
  
  # Remove empty strings
  report <- report[report != ""]
  
  if (!is.null(output_file)) {
    writeLines(report, output_file)
    message("Diagnostic report saved to: ", output_file)
  }
  
  invisible(report)
}