% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptive_sampling.R
\name{combine_job_results}
\alias{combine_job_results}
\title{Combine Results from Individual Job Files}
\usage{
combine_job_results(scenario_name, parallel_jobs, output_dir = NULL)
}
\arguments{
\item{scenario_name}{Character. Base scenario name used for all jobs.}

\item{parallel_jobs}{Integer. Number of parallel jobs that were run.}

\item{output_dir}{Character. Directory for output files. If NULL, uses current directory.}
}
\value{
Logical indicating success (TRUE) or failure (FALSE).
}
\description{
Helper function to combine results from multiple parallel job output files
into a single consolidated file.
}
\keyword{internal}
