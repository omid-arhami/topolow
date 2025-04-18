% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{calculate_diagnostics}
\alias{calculate_diagnostics}
\title{Calculate Adaptive Monte Carlo Sampling Diagnostics}
\usage{
calculate_diagnostics(chain_files, mutual_size = 500)
}
\arguments{
\item{chain_files}{Character vector of paths to CSV files containing chains}

\item{mutual_size}{Integer number of samples to use from end of each chain}
}
\value{
List containing:
\item{rhat}{R-hat statistic for each parameter}
\item{ess}{Effective sample size for each parameter}
}
\description{
Calculates standard Adaptive Monte Carlo Sampling diagnostics including R-hat (potential scale reduction)
and effective sample size for multiple chains. Can be used with any iterative
sampling or optimization procedure that produces chain-like output.
}
\examples{
\dontrun{
chain_files <- c("chain1.csv", "chain2.csv", "chain3.csv")
diag <- calculate_diagnostics(chain_files, mutual_size = 1000)
print(diag)  # Shows R-hat and ESS
plot(diag)   # Creates density plots
print(diag$rhat) # Should be close to 1
print(diag$ess)  # Should be large enough (>400) for reliable inference
}
}
