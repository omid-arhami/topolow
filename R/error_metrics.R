# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/error_metrics.R

#' Error calculation and validation metrics for topolow


# Newed
#' Calculate Comprehensive Error Metrics
#'
#' @description
#' Computes a comprehensive set of error metrics (in-sample, out-of-sample, completeness)
#' between predicted and true dissimilarities for model evaluation.
#'
#' @details
#' Input requirements and constraints:
#' * All input matrices must have matching dimensions.
#' * Row and column names must be consistent across matrices.
#' * NAs are allowed and handled appropriately.
#' * Threshold indicators (< or >) in the input matrix are processed correctly.
#'
#' When `input_dissimilarities` is provided, it represents the training data where some
#' values have been set to NA to create a holdout set. This allows calculation of:
#' * In-sample errors: for data available during training
#' * Out-of-sample errors: for data held out during training
#'
#' When `input_dissimilarities` is NULL (default), all errors are treated as in-sample
#' since no data was held out.
#'
#' @param predicted_dissimilarities Matrix of predicted dissimilarities from the model.
#' @param true_dissimilarities Matrix of true, ground-truth dissimilarities.
#' @param input_dissimilarities Matrix of input dissimilarities, which may contain NAs
#'        and is used to identify the pattern of missing values for out-of-sample error calculation.
#'        Optional - if not provided, defaults to `true_dissimilarities` (no holdout set).
#' @return A list containing:
#'   \item{report_df}{A `data.frame` with detailed error metrics for each point-pair, including `InSampleError`, `OutSampleError`, and their percentage-based counterparts.}
#'   \item{Completeness}{A single numeric value representing the completeness statistic, which is the fraction of validation points for which a prediction could be made.}
#'
#' @examples
#' # Example 1: Normal evaluation (no cross-validation)
#' true_mat <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, 3)
#' pred_mat <- true_mat + rnorm(9, 0, 0.1)  # Add some noise
#'
#' # Evaluate all predictions (input_dissimilarities defaults to true_dissimilarities)
#' errors1 <- error_calculator_comparison(pred_mat, true_mat)
#'
#' # Example 2: Cross-validation evaluation
#' input_mat <- true_mat
#' input_mat[1, 3] <- input_mat[3, 1] <- NA  # Create holdout set
#'
#' # Evaluate with train/test split
#' errors2 <- error_calculator_comparison(pred_mat, true_mat, input_mat)
#'
#' @importFrom dplyr %>% select
#' @export
error_calculator_comparison <- function(predicted_dissimilarities, true_dissimilarities, input_dissimilarities = NULL) {
  # Validate inputs
  if (!is.matrix(predicted_dissimilarities) || !is.matrix(true_dissimilarities)) {
    stop("predicted_dissimilarities and true_dissimilarities must be matrices")
  }

  # Default input_dissimilarities to true_dissimilarities if not provided
  if (is.null(input_dissimilarities)) {
    input_dissimilarities <- true_dissimilarities
  }

  if (!is.matrix(input_dissimilarities)) {
    stop("input_dissimilarities must be a matrix")
  }

  if (!all(dim(predicted_dissimilarities) == dim(true_dissimilarities)) ||
      !all(dim(predicted_dissimilarities) == dim(input_dissimilarities))) {
    stop("All matrices must have the same dimensions")
  }

  # Ensure consistent ordering of matrices only if they have row/column names
  if (!is.null(rownames(predicted_dissimilarities)) && !is.null(rownames(true_dissimilarities)) &&
      !is.null(colnames(predicted_dissimilarities)) && !is.null(colnames(true_dissimilarities))) {

    row_order <- order(factor(rownames(predicted_dissimilarities), levels = rownames(true_dissimilarities)))
    col_order <- order(factor(colnames(predicted_dissimilarities), levels = colnames(true_dissimilarities)))
    predicted_dissimilarities <- predicted_dissimilarities[row_order, col_order]

    if (!all(rownames(true_dissimilarities) == rownames(predicted_dissimilarities)) ||
        !all(colnames(true_dissimilarities) == colnames(predicted_dissimilarities))) {
      stop("Row and column names must match between matrices after ordering")
    }
  }

  # Convert matrices to numeric vectors, handling potential threshold strings
  input_vec <- suppressWarnings(as.numeric(input_dissimilarities))
  truth_vec <- suppressWarnings(as.numeric(true_dissimilarities))
  pred_vec <- as.vector(predicted_dissimilarities)

  # Identify the pattern of missing values from the input matrix
  missing_pattern <- is.na(input_vec)

  # Separate predictions for observed (in-sample) and missing (out-of-sample) points
  pred_for_observed <- pred_vec
  pred_for_observed[missing_pattern] <- NA

  pred_for_missing <- pred_vec
  pred_for_missing[!missing_pattern] <- NA

  # Construct a report data frame
  report_df <- data.frame(
    True_dissimilarity = truth_vec,
    Pred_for_observed = pred_for_observed,
    Pred_for_missing = pred_for_missing
  )

  # Calculate absolute errors
  report_df$InSampleError <- report_df$True_dissimilarity - report_df$Pred_for_observed
  report_df$OutSampleError <- report_df$True_dissimilarity - report_df$Pred_for_missing

  # Calculate percentage errors, avoiding division by zero
  non_zero_mask <- !is.na(report_df$True_dissimilarity) & report_df$True_dissimilarity > 0

  report_df$InSamplePercentageError <- NA
  report_df$InSamplePercentageError[non_zero_mask] <-
    (report_df$InSampleError[non_zero_mask] / report_df$True_dissimilarity[non_zero_mask]) * 100

  report_df$OutSamplePercentageError <- NA
  report_df$OutSamplePercentageError[non_zero_mask] <-
    (report_df$OutSampleError[non_zero_mask] / report_df$True_dissimilarity[non_zero_mask]) * 100

  # Clean up the final data frame
  report_df <- report_df %>%
    dplyr::select(-c(True_dissimilarity, Pred_for_observed, Pred_for_missing))

  # Calculate Completeness on the validation set (out-of-sample)
  validation_count <- sum(!is.na(truth_vec[missing_pattern]))
  pred_for_validation_count <- sum(!is.na(report_df$OutSampleError))
  Completeness <- if (validation_count > 0) pred_for_validation_count / validation_count else {
    # When no holdout set exists, completeness should reflect prediction availability
    total_possible <- sum(!is.na(truth_vec))
    total_predictions <- sum(!is.na(pred_vec))
    if (total_possible > 0) total_predictions / total_possible else 0
  }

  return(list(
    report_df = report_df,
    Completeness = Completeness
  ))
}


# Newed
#' Calculate Prediction Interval for Dissimilarity Estimates
#'
#' @description
#' Computes prediction intervals for the estimated dissimilarities based on
#' residual variation between true and predicted values.
#'
#' @param dissimilarity_matrix Matrix of true dissimilarities.
#' @param predicted_dissimilarity_matrix Matrix of predicted dissimilarities.
#' @param confidence_level The confidence level for the interval (default: 0.95).
#' @return A single numeric value representing the margin of error for the prediction interval.
#' @importFrom stats sd qt
#' @export
calculate_prediction_interval <- function(dissimilarity_matrix, predicted_dissimilarity_matrix, confidence_level = 0.95) {
  # Validate inputs
  if (!is.matrix(dissimilarity_matrix) || !is.matrix(predicted_dissimilarity_matrix)) {
    stop("Both inputs must be matrices")
  }

  if (!identical(dim(dissimilarity_matrix), dim(predicted_dissimilarity_matrix))) {
    stop("Matrices must have the same dimensions")
  }

  if (!is.numeric(confidence_level) ||
      confidence_level <= 0 ||
      confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1")
  }

  if (nrow(dissimilarity_matrix) < 2) {
    stop("Matrices must have at least 2 rows/columns for meaningful prediction intervals")
  }

  # Convert matrices to numeric vectors, handling threshold indicators
  true_dissimilarities <- extract_numeric_values(as.vector(dissimilarity_matrix))
  predicted_dissimilarities <- as.vector(predicted_dissimilarity_matrix)

  # Keep only valid, paired observations
  valid_indices <- !is.na(true_dissimilarities) & !is.na(predicted_dissimilarities)
  if(sum(valid_indices) < 2) { # Need at least 2 points to calculate sd
    stop("Fewer than 2 valid paired observations found. Cannot calculate prediction interval.")
  }

  true_dissimilarities <- true_dissimilarities[valid_indices]
  predicted_dissimilarities <- predicted_dissimilarities[valid_indices]

  # Calculate residuals and their standard deviation
  residuals <- true_dissimilarities - predicted_dissimilarities
  residual_sd <- stats::sd(residuals)

  # Calculate critical value from t-distribution and the margin of error
  critical_value <- stats::qt((1 + confidence_level) / 2, df = length(true_dissimilarities) - 1)
  margin_of_error <- critical_value * residual_sd

  return(margin_of_error)
}
