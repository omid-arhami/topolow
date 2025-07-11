% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error_metrics.R
\name{calculate_prediction_interval}
\alias{calculate_prediction_interval}
\title{Calculate prediction interval for distance estimates}
\usage{
calculate_prediction_interval(
  distance_matrix,
  p_dist_mat,
  confidence_level = 0.95
)
}
\arguments{
\item{distance_matrix}{Matrix of true distances}

\item{p_dist_mat}{Matrix of predicted distances}

\item{confidence_level}{Confidence level for interval (default: 0.95)}
}
\value{
A single numeric value representing the margin of error for the prediction interval.
}
\description{
Computes prediction intervals for the estimated distances based on
residual variation between true and predicted values.
}
