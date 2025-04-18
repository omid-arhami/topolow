% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{check_gaussian_convergence}
\alias{check_gaussian_convergence}
\title{Check Multivariate Gaussian Convergence}
\usage{
check_gaussian_convergence(data, window_size = 300, tolerance = 0.01)
}
\arguments{
\item{data}{Matrix or data frame of samples where columns are parameters}

\item{window_size}{Integer size of sliding window for statistics}

\item{tolerance}{Numeric convergence threshold for relative changes}
}
\value{
List containing:
\item{converged}{Logical indicating if convergence achieved}
\item{mean_converged}{Logical for mean convergence}
\item{cov_converged}{Logical for covariance convergence}
\item{final_mean}{Vector of final mean values}
\item{final_cov}{Final covariance matrix}
\item{mean_history}{Matrix of mean values over iterations}
\item{cov_changes}{Vector of covariance changes}
}
\description{
Assesses convergence of multivariate samples by monitoring changes in mean
vector and covariance matrix over a sliding window. Useful for checking
stability of parameter distributions in optimization or sampling.
}
\examples{
\dontrun{
data <- read.csv("chain_data.csv")
conv_results <- check_gaussian_convergence(data)
print(conv_results)  # Shows summary
plot(conv_results)   # Creates convergence plots
}
}
