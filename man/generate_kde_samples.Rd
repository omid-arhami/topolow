% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptive_sampling.R
\name{generate_kde_samples}
\alias{generate_kde_samples}
\title{Generate New Parameter Samples Using KDE}
\usage{
generate_kde_samples(samples, n, epsilon = 0)
}
\arguments{
\item{samples}{Data frame of previous samples with parameters and NLL}

\item{n}{Integer number of samples to generate}

\item{epsilon}{Numeric probability of wider bandwidth sampling}
}
\value{
Data frame of n new samples
}
\description{
Generates new parameter samples using weighted kernel density estimation
for each parameter independently.
}
\keyword{internal}
