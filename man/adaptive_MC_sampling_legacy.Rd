% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptive_sampling.R
\name{adaptive_MC_sampling_legacy}
\alias{adaptive_MC_sampling_legacy}
\title{Perform Adaptive Monte Carlo Sampling}
\usage{
adaptive_MC_sampling_legacy(
  samples_file,
  distance_matrix,
  n_iter = 1,
  batch_size = 1,
  max_iter,
  relative_epsilon,
  folds = 20,
  num_cores,
  scenario_name,
  replace_csv
)
}
\arguments{
\item{samples_file}{Path to CSV with initial samples}

\item{distance_matrix}{Distance matrix to fit}

\item{n_iter}{Number of sampling iterations}

\item{batch_size}{Samples per iteration}

\item{max_iter}{Maximum optimization iterations}

\item{relative_epsilon}{Convergence threshold}

\item{folds}{Number of CV folds}

\item{num_cores}{Number of cores for parallel processing}

\item{scenario_name}{Name for output files}

\item{replace_csv}{Whether to replace existing CSV}
}
\value{
Data frame of samples with evaluated likelihoods
}
\description{
Main function implementing adaptive Monte Carlo sampling (https://www.sciencedirect.com/science/article/pii/0167473088900203)
to explore parameter space. Updates sampling distribution based on evaluated likelihoods.
}
