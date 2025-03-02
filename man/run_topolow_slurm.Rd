% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algorithm_comparison_helpers.R
\name{run_topolow_slurm}
\alias{run_topolow_slurm}
\title{Run Topolow via SLURM  (LEGACY)}
\usage{
run_topolow_slurm(
  truth_matrix,
  input_matrix,
  best_params,
  opt_params,
  scenario_name,
  fold,
  cider,
  time = "8:00:00",
  memory = "18G",
  matrix_dir,
  results_dir
)
}
\arguments{
\item{truth_matrix}{Matrix of true distances}

\item{input_matrix}{Matrix for optimization}

\item{best_params}{List of optimal parameters}

\item{opt_params}{List of optimization parameters}

\item{scenario_name}{Character scenario identifier}

\item{fold}{Integer fold number}

\item{cider}{Logical; whether to use cider queue}
}
\value{
List of performance metrics
}
\description{
Submits Topolow job to SLURM cluster and collects results.
}
\keyword{internal}
