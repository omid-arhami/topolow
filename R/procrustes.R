# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# Licensed under MIT License 
# R/procrustes.R

#' Procrustes Analysis Functions
#' 
#' @description
#' Functions for comparing different map configurations using Procrustes analysis.
#' These functions help assess:
#' - Statistical significance of differences between maps
#' - Quantitative measures of map differences
#' - Stability of mapping solutions
#'
#' @keywords internal
"_PACKAGE"

#' Calculate Statistical Significance Between Maps Using Procrustes Analysis
#'
#' @description
#' Performs Procrustes analysis between two maps and calculates statistical 
#' significance of their differences using permutation tests. Handles common 
#' data cleaning steps like removing missing values and ensuring comparable 
#' point sets.
#'
#' @param map1 Data frame with coordinates from first map (must have X, X.1 columns)
#' @param map2 Data frame with coordinates from second map (must have X, X.1 columns)
#' @return Numeric p-value from Procrustes permutation test
#' @examples
#' \dontrun{
#' map1 <- read.csv("map1_coords.csv")
#' map2 <- read.csv("map2_coords.csv")
#' p_val <- calculate_procrustes_significance(map1, map2)
#' }
#' @importFrom vegan protest
#' @export
calculate_procrustes_significance <- function(map1, map2) {
  # Ensure no missing or infinite values
  map1 <- map1 %>% filter(!is.na(X) & !is.na(X.1) & is.finite(X) & is.finite(X.1))
  map2 <- map2 %>% filter(!is.na(X) & !is.na(X.1) & is.finite(X) & is.finite(X.1))
  
  # Ensure both maps have the same number of points
  common_names <- intersect(map1$name, map2$name)
  map1 <- map1 %>% filter(name %in% common_names)
  map2 <- map2 %>% filter(name %in% common_names)

  procrustes_result <- protest(map1[, c("X", "X.1")], map2[, c("X", "X.1")])
  return(procrustes_result$signif)
}

#' Calculate Procrustes Difference Between Maps
#'
#' @description 
#' Computes the quantitative difference between two maps using Procrustes analysis.
#' The difference is calculated as the sum of squared differences after optimal
#' rotation and scaling.
#'
#' @param map1 Data frame with coordinates from first map (must have X, X.1 columns)
#' @param map2 Data frame with coordinates from second map (must have X, X.1 columns) 
#' @return Numeric sum of squared differences after Procrustes transformation
#' @examples
#' \dontrun{
#' map1 <- read.csv("map1_coords.csv")
#' map2 <- read.csv("map2_coords.csv")
#' diff <- calculate_procrustes_difference(map1, map2)
#' }
#' @importFrom vegan procrustes
#' @export
calculate_procrustes_difference <- function(map1, map2) {
  # Validate inputs are data frames
  if (!is.data.frame(map1) || !is.data.frame(map2)) {
    stop("Both inputs must be data frames")
  }
  
  # Validate coordinate columns
  req_cols <- c("X", "X.1")
  if (!all(req_cols %in% names(map1)) || !all(req_cols %in% names(map2))) {
    stop("Maps must have X and X.1 coordinate columns")
  }
  
  # Validate numeric coordinates
  if (!all(sapply(map1[req_cols], is.numeric)) || 
      !all(sapply(map2[req_cols], is.numeric))) {
    stop("Coordinates must be numeric")
  }

  # Ensure no missing or infinite values
  map1 <- map1 %>% filter(!is.na(X) & !is.na(X.1) & is.finite(X) & is.finite(X.1))
  map2 <- map2 %>% filter(!is.na(X) & !is.na(X.1) & is.finite(X) & is.finite(X.1))
  
  # Ensure both maps have the same number of points, or find the intersection
  common_names <- intersect(map1$name, map2$name)
  map1 <- map1 %>% filter(name %in% common_names)
  map2 <- map2 %>% filter(name %in% common_names)

  procrustes_result <- procrustes(map1[, c("X", "X.1")], map2[, c("X", "X.1")], 
                                 scale = TRUE)
  return(procrustes_result$ss)
}