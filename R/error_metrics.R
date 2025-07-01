# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/error_metrics.R

#' Error calculation and validation metrics for topolow
#' 
#' @description
#' This file contains functions for calculating error metrics and validation statistics
#' between predicted and true distance matrices. Functions handle missing values and
#' special cases like threshold measurements.
#'
#' @keywords internal
"_PACKAGE"

#' Calculate comprehensive error metrics between predicted and true distances
#'
#' @description
#' Computes various error metrics including in-sample and out-of-sample errors,
#' and Completeness statistics for model evaluation.
#'
#' @details
#' Input requirements and constraints:
#' * Matrices must have matching dimensions
#' * Row and column names must be consistent between matrices
#' * NAs are allowed and handled appropriately
#' * Threshold indicators (< or >) in input matrix are processed correctly
#'
#' @param p_dist_mat Matrix of predicted distances
#' @param truth_matrix Matrix of true distances 
#' @param input_matrix Matrix of input distances (may contain NAs and is used to find the NAs' pattern)
#' @return A list containing:
#'   \item{report_df}{A `data.frame` with detailed error metrics for each point-pair, including `InSampleError`, `OutSampleError`, and their percentage-based counterparts.}
#'   \item{Completeness}{A single numeric value representing the completeness statistic, which is the fraction of validation points for which a prediction could be made.}
#'   
#' @importFrom dplyr %>% select
#' @export
error_calculator_comparison <- function(p_dist_mat, truth_matrix, input_matrix) {
  # Validate inputs
  if (!is.matrix(p_dist_mat) || !is.matrix(truth_matrix) || !is.matrix(input_matrix)) {
    stop("All inputs must be matrices")
  }
  
  if (!all(dim(p_dist_mat) == dim(truth_matrix)) || 
      !all(dim(p_dist_mat) == dim(input_matrix))) {
    stop("All matrices must have same dimensions")
  }
  
  # Reorder rows and columns of matrix2 to match matrix1
  row_order <- order(factor(rownames(p_dist_mat), levels = rownames(truth_matrix)))
  col_order <- order(factor(colnames(p_dist_mat), levels = colnames(truth_matrix)))
  
  # Apply the orders to matrix2
  p_dist_mat <- p_dist_mat[row_order, col_order]
  
  if (!all(rownames(truth_matrix) == rownames(p_dist_mat)) ||
      !all(colnames(truth_matrix) == colnames(p_dist_mat))) {
    stop("Row/column names must match between matrices")
  }
  
  # Convert to numeric, handling threshold values
  input_matrix <- suppressWarnings(as.numeric(input_matrix))
  truth_matrix <- suppressWarnings(as.numeric(truth_matrix))
  
  # Find validation set (missing values in input)
  missing_pattern <- is.na(input_matrix)
  
  # Create matrices for observed vs predicted comparisons
  pred_for_observeds <- p_dist_mat
  pred_for_observeds[missing_pattern] <- NA
  
  pred_for_missing <- p_dist_mat
  pred_for_missing[!missing_pattern] <- NA
  
  # Construct report dataframe
  report_df <- data.frame(
    True_distance = as.vector(truth_matrix),
    Pred_for_observeds = as.vector(pred_for_observeds),
    Pred_for_missing = as.vector(pred_for_missing)
  )
  
  # Calculate absolute errors 
  report_df$InSampleError <- report_df$True_distance - report_df$Pred_for_observeds
  report_df$OutSampleError <- report_df$True_distance - report_df$Pred_for_missing
  
  # Calculate percentage errors - WITH FIX
  non_zero_mask <- !is.na(report_df$True_distance) & report_df$True_distance > 0
  
  report_df$InSamplePercentageError <- NA
  report_df$InSamplePercentageError[non_zero_mask] <- 
    (report_df$InSampleError[non_zero_mask] / report_df$True_distance[non_zero_mask]) * 100

  report_df$OutSamplePercentageError <- NA
  report_df$OutSamplePercentageError[non_zero_mask] <- 
    (report_df$OutSampleError[non_zero_mask] / report_df$True_distance[non_zero_mask]) * 100

  report_df <- report_df %>% 
    dplyr::select(-c(True_distance, Pred_for_observeds, Pred_for_missing))
  
  # Calculate Completeness on validation set
  validation_count <- sum(!is.na(truth_matrix[missing_pattern]))
  pred_for_validation_count <- sum(!is.na(truth_matrix) & !is.na(pred_for_missing))
  Completeness <- pred_for_validation_count / validation_count
  
  return(list(
    report_df = report_df,
    Completeness = Completeness
  ))
}


#' Calculate prediction interval for distance estimates
#'
#' Computes prediction intervals for the estimated distances based on
#' residual variation between true and predicted values.
#'
#' @param distance_matrix Matrix of true distances
#' @param p_dist_mat Matrix of predicted distances 
#' @param confidence_level Confidence level for interval (default: 0.95)
#' @return A single numeric value representing the margin of error for the prediction interval.
#' @importFrom stats sd qt
#' @export
calculate_prediction_interval <- function(distance_matrix, p_dist_mat, confidence_level = 0.95) {
  # Validate inputs
  if (!is.matrix(distance_matrix) || !is.matrix(p_dist_mat)) {
    stop("Both inputs must be matrices")
  }
  
  if (!identical(dim(distance_matrix), dim(p_dist_mat))) {
    stop("Matrices must have same dimensions")
  }
  
  if (!is.numeric(confidence_level) || 
      confidence_level <= 0 || 
      confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1")
  }
  
  # simple numeric conversion to remove NA and thresholded values
  true_distances <- suppressWarnings(as.numeric(as.vector(distance_matrix)))
  predicted_distances <- as.vector(p_dist_mat)
  
  # Remove missing values
  valid_indices <- !is.na(true_distances) & !is.na(predicted_distances)
  if(sum(valid_indices) == 0) {
    stop("No valid paired observations found")
  }
  
  true_distances <- true_distances[valid_indices]
  predicted_distances <- predicted_distances[valid_indices]
  
  # Calculate residuals and their standard deviation
  residuals <- true_distances - predicted_distances
  residual_sd <- sd(residuals)
  
  # Calculate critical value and margin of error
  critical_value <- qt((1 + confidence_level) / 2, df = length(true_distances) - 1)
  margin_of_error <- critical_value * residual_sd
  
  return(margin_of_error)
}
