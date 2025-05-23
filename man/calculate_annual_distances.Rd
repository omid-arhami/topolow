% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{calculate_annual_distances}
\alias{calculate_annual_distances}
\title{Calculate Annual Distance Metrics}
\usage{
calculate_annual_distances(df_coords, ndim, na.rm = TRUE)
}
\arguments{
\item{df_coords}{Data frame containing:
- V1...Vn coordinate columns
- year: Numeric years
- name: Point identifiers (will use rownames if missing)}

\item{ndim}{Number of coordinate dimensions}

\item{na.rm}{Logical indicating whether to remove NA values}
}
\value{
List containing:
\item{dist_data}{Data frame with columns:
\itemize{
\item year: Collection year
\item distance: Distance from previous year mean
}
}
\item{summary}{List with:
\itemize{
\item overall_mean: Mean distance across all years
\item overall_sd: Standard deviation of distances
}
}
}
\description{
Calculates year-over-year antigenic distances and statistics. Compares each point
to the mean coordinates of the previous year.
}
\examples{
\dontrun{
annual_stats <- calculate_annual_distances(coords, ndim=2)
print(annual_stats$summary$overall_mean)
}
}
