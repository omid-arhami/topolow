% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{new_dim_reduction_config}
\alias{new_dim_reduction_config}
\title{Dimension Reduction Configuration Class}
\usage{
new_dim_reduction_config(
  method = "pca",
  n_components = 2,
  scale = TRUE,
  center = TRUE,
  pca_params = list(tol = sqrt(.Machine$double.eps), rank. = NULL),
  umap_params = list(n_neighbors = 15, min_dist = 0.1, metric = "euclidean", n_epochs =
    200),
  tsne_params = list(perplexity = 30, mapping_max_iter = 1000, theta = 0.5),
  compute_loadings = FALSE,
  random_state = NULL
)
}
\arguments{
\item{method}{Dimension reduction method ("pca", "umap", "tsne")}

\item{n_components}{Number of components to compute}

\item{scale}{Scale the data before reduction}

\item{center}{Center the data before reduction}

\item{pca_params}{List of PCA-specific parameters}

\item{umap_params}{List of UMAP-specific parameters}

\item{tsne_params}{List of t-SNE-specific parameters}

\item{compute_loadings}{Compute and return loadings}

\item{random_state}{Random seed for reproducibility}
}
\value{
An S3 object of class \code{dim_reduction_config}, which is a list
containing the specified configuration parameters for dimensionality reduction.
}
\description{
S3 class for configuring dimension reduction parameters including
method selection and algorithm-specific parameters.
}
