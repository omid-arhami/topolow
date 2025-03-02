% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptive_sampling.R
\name{safe_likelihood_function}
\alias{safe_likelihood_function}
\title{Safe Wrapper for Likelihood Evaluation}
\usage{
safe_likelihood_function(...)
}
\arguments{
\item{...}{Arguments passed to likelihood_function}
}
\value{
Same as likelihood_function or NA if error
}
\description{
Wraps likelihood calculation in error handler to return NA on failure.
}
\keyword{internal}
