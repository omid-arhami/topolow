% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{long_to_matrix}
\alias{long_to_matrix}
\title{Convert Long Format Data to Distance Matrix}
\usage{
long_to_matrix(
  data,
  chnames,
  chorder = NULL,
  rnames,
  rorder = NULL,
  values_column,
  rc = FALSE,
  sort = FALSE
)
}
\arguments{
\item{data}{Data frame in long format}

\item{chnames}{Character. Name of column holding the challenge point names.}

\item{chorder}{Character. Optional name of column for challenge point ordering.}

\item{rnames}{Character. Name of column holding reference point names.}

\item{rorder}{Character. Optional name of column for reference point ordering.}

\item{values_column}{Character. Name of column containing distance/difference values. It should be from the nature of "distance" (e.g., antigenic distance or IC50), not "similarity" (e.g., HI Titer.)}

\item{rc}{Logical. If TRUE, reference points are treated as a subset of challenge
points. If FALSE, they are treated as distinct sets. Default is FALSE.}

\item{sort}{Logical. Whether to sort rows/columns by chorder/rorder. Default FALSE.}
}
\value{
A symmetric \code{matrix} of distances with row and column names corresponding
to the unique points in the data. \code{NA} values represent unmeasured pairs.
}
\description{
Converts a dataset from long format to a symmetric distance matrix. The function
handles antigenic cartography data where measurements may exist between antigens
and antisera points. Row and column names can be optionally sorted by a time
variable.
}
\details{
The function expects data in long format with at least three columns:
\itemize{
\item A column for challenge point names
\item A column for reference point names
\item A column containing the distance/difference values
}

Optionally, ordering columns can be provided to sort the output matrix.
The 'rc' parameter determines how to handle shared names between references
and challenges.
}
\examples{
data <- data.frame(
  antigen = c("A", "B", "A"),
  serum = c("X", "X", "Y"), 
  distance = c(2.5, 1.8, 3.0),
  year = c(2000, 2001, 2000)
)

# Basic conversion
mat <- long_to_matrix(data, 
                     chnames = "antigen",
                     rnames = "serum",
                     values_column = "distance")
                     
# With sorting by year
mat_sorted <- long_to_matrix(data,
                            chnames = "antigen",
                            chorder = "year",
                            rnames = "serum", 
                            rorder = "year",
                            values_column = "distance",
                            sort = TRUE)
}
