% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{coordinates_to_matrix}
\alias{coordinates_to_matrix}
\title{Convert coordinates to distance matrix}
\usage{
coordinates_to_matrix(positions)
}
\arguments{
\item{positions}{Matrix of coordinates where rows are points and columns are dimensions;
It can be a matrix or a data frame.}
}
\value{
A symmetric \code{matrix} of pairwise Euclidean distances between points.
}
\description{
Calculates pairwise Euclidean distances between points in coordinate space
}
