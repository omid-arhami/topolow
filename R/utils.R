# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/utils.R

#' Utility functions for the topolow package

# New
#' Extract Numeric Values from Mixed Data
#'
#' @description
#' Extracts numeric values from data that may contain threshold indicators
#' (e.g., "<10", ">1280") or regular numeric values.
#'
#' @param x A vector that may contain numeric values, character strings with
#'   threshold indicators, or a mix of both.
#' @return A numeric vector with threshold indicators converted to their
#'   numeric equivalents.
#' @examples
#' # Mixed data with threshold indicators
#' mixed_data <- c(10, 20, "<5", ">100", 50)
#' extract_numeric_values(mixed_data)
#'
#' @export
extract_numeric_values <- function(x) {
  sapply(x, function(val) {
    if (is.character(val)) {
      if (grepl("^<", val)) {
        as.numeric(sub("<", "", val))
      } else if (grepl("^>", val)) {
        as.numeric(sub(">", "", val))
      } else {
        as.numeric(val)
      }
    } else {
      as.numeric(val)
    }
  }, USE.NAMES = FALSE)
}


# Newed
#' Create Cross-Validation Folds for a Dissimilarity Matrix
#'
#' @description
#' Creates k-fold cross-validation splits from a dissimilarity matrix while maintaining
#' symmetry. Each fold in the output consists of a training matrix (with some
#' values masked as `NA`) and a corresponding ground truth matrix for validation.
#'
#' @param dissimilarity_matrix The input dissimilarity matrix, which may contain noise.
#' @param ground_truth_matrix An optional, noise-free dissimilarity matrix to be used as the ground truth for evaluation. If `NULL`, the input `dissimilarity_matrix` is used as the truth.
#' @param n_folds The integer number of folds to create.
#' @param random_seed An optional integer to set the random seed for reproducibility.
#'
#' @return A list of length `n_folds`. Each element of the list is itself a list
#'   containing two matrices: `truth` (the ground truth for that fold) and `train`
#'   (the training matrix with `NA` values for validation).
#' @note This function has breaking changes from previous versions:
#'   \itemize{
#'     \item Parameter \code{truth_matrix} renamed to \code{dissimilarity_matrix}
#'     \item Parameter \code{no_noise_truth} renamed to \code{ground_truth_matrix}  
#'     \item Return structure now uses named elements (\code{$truth}, \code{$train})
#'   }
#' @examples
#' # Create a sample dissimilarity matrix
#' d_mat <- matrix(runif(100), 10, 10)
#' diag(d_mat) <- 0
#'
#' # Create 5-fold cross-validation splits
#' folds <- create_cv_folds(d_mat, n_folds = 5, random_seed = 123)
#' @export
create_cv_folds <- function(dissimilarity_matrix, ground_truth_matrix = NULL,
                            n_folds = 10, random_seed = NULL) {
  # --- Input Validation and Setup ---
  if (!is.null(random_seed)) {
    if (!is.numeric(random_seed) || random_seed != round(random_seed)) {
      stop("`random_seed` must be an integer.")
    }
    set.seed(random_seed)
  }

  if (!is.matrix(dissimilarity_matrix)) {
    stop("`dissimilarity_matrix` must be a matrix.")
  }

  if (!is.null(ground_truth_matrix)) {
    if (!is.matrix(ground_truth_matrix)) {
      stop("`ground_truth_matrix` must be NULL or a matrix.")
    }
    if (!identical(dim(dissimilarity_matrix), dim(ground_truth_matrix))) {
      stop("`dissimilarity_matrix` and `ground_truth_matrix` must have the same dimensions.")
    }
  }

  if (!is.numeric(n_folds) || n_folds < 2 || n_folds != round(n_folds)) {
    stop("`n_folds` must be an integer greater than or equal to 2.")
  }

  if (n_folds > nrow(dissimilarity_matrix)) {
    stop("`n_folds` cannot be larger than the number of rows in the matrix.")
  }

  # Use the perfect and full matrix as the ground truth for evaluation if provided
  eval_truth <- if (!is.null(ground_truth_matrix)) ground_truth_matrix else dissimilarity_matrix

  # --- Fold Creation ---
  # Initialize a list to store the [truth, train] pairs for each fold
  matrix_list <- vector("list", n_folds)
  for(i in 1:n_folds) {
    matrix_list[[i]] <- list(truth = eval_truth, train = NULL)
  }

  # Determine the number of elements to hold out in each fold
  num_elements <- sum(!is.na(dissimilarity_matrix))
  holdout_size <- floor(num_elements / (n_folds * 2)) # Factor of 2 for symmetry

  # Create a temporary copy to track available elements for sampling
  sampling_pool <- dissimilarity_matrix

  for(i in 1:n_folds) {
    # Ensure there are enough elements left to sample for the holdout set
    if (sum(!is.na(sampling_pool)) < holdout_size) {
        warning("Could not create all requested folds due to data sparsity. Returning fewer folds.")
        return(matrix_list[1:(i-1)])
    }

    # Sample validation indices from the remaining available elements
    random_indices <- sample(which(!is.na(sampling_pool)), size = holdout_size)

    # Create the training matrix for this fold by masking the holdout set
    train_matrix <- dissimilarity_matrix
    for(index in random_indices) {
      # Convert the linear index to row/column coordinates
      row <- (index - 1) %/% nrow(dissimilarity_matrix) + 1
      col <- (index - 1) %% ncol(dissimilarity_matrix) + 1

      # Mask the validation entries symmetrically
      train_matrix[row, col] <- NA
      train_matrix[col, row] <- NA

      # Remove the selected elements from the sampling pool for subsequent folds
      sampling_pool[row, col] <- NA
      # Ensure symmetry is maintained in the sampling pool as well
      sampling_pool[col, row] <- NA
    }

    matrix_list[[i]]$train <- train_matrix


  }

  return(matrix_list)
}


# Newed
#' Save ggplot with white background
#'
#' @description
#' Wrapper around `ggplot2::ggsave` that ensures a white background by default.
#' @importFrom ggplot2 ggsave
#' @inheritParams ggplot2::ggsave
#' @return No return value, called for side effects.
#' @export
ggsave_white_bg <- function(..., bg = 'white') {
  ggplot2::ggsave(..., bg = bg)
}


# Newed
#' Color Palettes
#'
#' @description
#' Predefined color palettes optimized for visualization.
#'
#' @name color_palettes
NULL

#' @rdname color_palettes
#' @export
c25 <- c(
  "purple", "green1", "blue1", "gold1", "red",
  "darkturquoise", "darkorange", "skyblue2", "green4", "maroon",
  "yellow3", "gray40", "hotpink", "darkorange4", "deeppink1",
  "khaki2", "palegreen2", "dodgerblue2", "brown", "orchid1"
)
