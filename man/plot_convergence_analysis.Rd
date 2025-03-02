% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_convergence_analysis}
\alias{plot_convergence_analysis}
\title{Plot Convergence Analysis Results}
\usage{
plot_convergence_analysis(conv_results, param_names)
}
\arguments{
\item{conv_results}{List output from check_gaussian_convergence()}

\item{param_names}{Character vector of parameter names}
}
\value{
A grid of plots showing convergence metrics
}
\description{
Visualizes convergence diagnostics including parameter mean trajectories
and covariance changes over iterations. Covariance norm changes measured by
Frobenius norm (also called Hilbert-Schmidt norm), the square root of the
sum of the absolute squares of all matrix elements = sqrt(sum|a_ij|²)
}
\examples{
\dontrun{
results <- check_gaussian_convergence(chain_data)
plot_convergence_analysis(results, c("mu", "sigma"))
}
}
