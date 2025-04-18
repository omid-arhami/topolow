% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{create_topolow_map}
\alias{create_topolow_map}
\title{Main TopoLow algorithm implementation}
\usage{
create_topolow_map(
  distance_matrix,
  ndim,
  mapping_max_iter,
  k0,
  cooling_rate,
  c_repulsion,
  relative_epsilon = 1e-04,
  convergence_counter = 5,
  initial_positions = NULL,
  write_positions_to_csv = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{distance_matrix}{Matrix. Square, symmetric distance matrix. Can contain NA values
for missing measurements and character strings with < or > prefixes for thresholded
measurements.}

\item{ndim}{Integer. Number of dimensions for the embedding space.}

\item{mapping_max_iter}{Integer. Maximum number of map optimization iterations.}

\item{k0}{Numeric. Initial spring constant controlling spring forces.}

\item{cooling_rate}{Numeric. Rate of spring constant decay per iteration (0 < cooling_rate < 1).}

\item{c_repulsion}{Numeric. Repulsion constant controlling repulsive forces.}

\item{relative_epsilon}{Numeric. Convergence threshold for relative change in error.
Default is 1e-4.}

\item{convergence_counter}{Integer. Number of iterations below threshold before declaring
convergence. Default is 5.}

\item{initial_positions}{Matrix or NULL. Optional starting coordinates. If NULL,
random initialization is used. Matrix should have nrow = nrow(distance_matrix)
and ncol = ndim.}

\item{write_positions_to_csv}{Logical. Whether to save point positions to CSV file.
Default is TRUE.}

\item{verbose}{Logical. Whether to print progress messages. Default is TRUE.}
}
\value{
A list with class "topolow" containing:
\itemize{
\item positions: Matrix of optimized point coordinates
\item est_distances: Matrix of distances in the optimized configuration
\item mae: Mean absolute error between target and optimized distances
\item iter: Number of iterations performed
\item parameters: List of input parameters used
\item convergence: List with convergence status and final error
}
}
\description{
TopoLow (Topological Optimization for Low-Dimensional Mapping) optimizes point positions in n-dimensional
space to match a target distance matrix. The algorithm uses a physics-inspired approach with
spring and repulsive forces to find optimal point configurations while handling missing
and thresholded measurements.
}
\details{
The algorithm iteratively updates point positions using:
\itemize{
\item Spring forces between points with measured distances
\item Repulsive forces between points without measurements
\item Modified forces for thresholded measurements (< or >)
\item Adaptive spring constant that decays over iterations
\item Convergence monitoring based on relative error change
}

Valid parameter ranges and constraints:
\itemize{
\item ndim: Positive integer, typically 2-20.
\item k0: Initial spring constant, positive numeric > 0. Typical range: 0.1-30
Controls initial force strength
\item cooling_rate: Spring and repulsion decay rate, numeric between 0 and 1. Typical range: 0.0001-0.1
Controls how quickly spring forces weaken
\item c_repulsion: Repulsion constant, positive numeric > 0. Typical range: 0.00001-0.1
Controls strength of repulsive forces
\item relative_epsilon: Positive numeric, typically 1e-9 to 1e-3
Smaller values require more iterations but give higher precision
\item convergence_counter: Positive integer, typically 5-20
Higher values do not necessarily lead to a better convergence
}
}
\examples{
# Create a simple distance matrix
dist_mat <- matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), nrow=3)

# Run TopoLow in 2D
result <- create_topolow_map(dist_mat, ndim=2, mapping_max_iter=1000, 
                      k0=1.0, cooling_rate=0.001, c_repulsion=0.1)
                      
# Plot results
plot(result$positions)

}
