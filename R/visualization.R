# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# License: free of charge access granted to any academic researcher to use this software for non-commercial, academic research purposes **only**.  Nobody may modify, distribute, sublicense, or publicly share the Software or any derivative works, until the paper is published by the original authors.  The Software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.  In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

# R/visualization.R

#' Visualization functions for the topolow package
#' 
#' @description
#' This file contains functions for visualizing topolow results including
#' dimension reduction plots and cluster visualizations. Supports multiple
#' plotting methods and customization options.
#' 
#' Functions handle:
#' - Temporal mapping visualizations
#' - Cluster mapping visualizations  
#' - 2D and 3D projections
#' - Multiple dimension reduction methods
#' - Interactive and static plots
#' - Diagnostic visualizations
#' - Monte Carlo analysis visualizations
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_minimal coord_fixed 
#'             ggsave theme element_text margin unit scale_shape_manual guides
#' @importFrom dplyr %>% select filter
#' @importFrom stats prcomp dist optimize
#' @importFrom umap umap
#' @importFrom rgl open3d plot3d axes3d rgl.snapshot rgl.close
#' @importFrom gridExtra grid.arrange
#' @importFrom scales comma
#' @keywords internal
"_PACKAGE"


#' Plot Aesthetic Configuration Class
#'
#' @description
#' S3 class for configuring plot visual aesthetics including points, colors,
#' labels and text elements.
#'
#' @param point_size Base point size
#' @param point_alpha Point transparency 
#' @param point_shapes Named vector of shapes for different point types
#' @param color_palette Color palette name or custom palette
#' @param gradient_colors List with low and high colors for gradients
#' @param show_labels Whether to show point labels
#' @param show_title Whether to show plot title (default: TRUE)
#' @param label_size Label text size 
#' @param title_size Title text size
#' @param subtitle_size Subtitle text size
#' @param axis_title_size Axis title text size
#' @param axis_text_size Axis text size
#' @param legend_text_size Legend text size
#' @param legend_title_size Legend title text size
#' @param show_legend Whether to show the legend
#' @param legend_position Legend position ("none", "right", "left", "top", "bottom")
#' @return An aesthetic_config object
#' @export
new_aesthetic_config <- function(
    point_size = 3.5,
    point_alpha = 0.8,
    point_shapes = c(antigen = 16, antiserum = 0),
    color_palette = c25,
    gradient_colors = list(low = "blue", high = "red"),
    show_labels = FALSE,
    show_title = TRUE,
    label_size = 3,
    title_size = 14,
    subtitle_size = 12,
    axis_title_size = 12,
    axis_text_size = 10,
    legend_text_size = 10,
    legend_title_size = 12,
    show_legend = TRUE,
    legend_position = "right"
) {
  # Validate point_shapes
  if (!is.numeric(point_shapes) || !all(c("antigen", "antiserum") %in% names(point_shapes))) {
    stop("point_shapes must be a named numeric vector with 'antigen' and 'antiserum' elements")
  }
  # Validate color palette
  if (is.character(color_palette) && length(color_palette) == 1) {
    # If a single string is provided, check if it's a valid palette name
    if (!exists(color_palette)) {
      stop("Invalid color palette name: ", color_palette, 
           ". Must be a valid palette object or vector of colors.")
    }
    color_palette <- get(color_palette)  # Get the actual palette
  }
  
  if (!is.vector(color_palette)) {
    stop("color_palette must be a vector of colors or a valid palette name")
  }
  
  config <- list(
    point_size = point_size,
    point_alpha = point_alpha,
    point_shapes = point_shapes,
    color_palette = color_palette,
    gradient_colors = gradient_colors,
    show_labels = show_labels,
    show_title = show_title,
    label_size = label_size,
    title_size = title_size,
    subtitle_size = subtitle_size,
    axis_title_size = axis_title_size,
    axis_text_size = axis_text_size,
    legend_text_size = legend_text_size,
    legend_title_size = legend_title_size,
    show_legend = show_legend,
    legend_position = legend_position
  )
  
  # Validate inputs
  stopifnot(
    is.numeric(point_size), point_size > 0,
    is.numeric(point_alpha), point_alpha >= 0, point_alpha <= 1,
    is.numeric(point_shapes),
    is.character(color_palette),
    is.list(gradient_colors), 
    all(c("low", "high") %in% names(gradient_colors)),
    is.logical(show_labels),
    is.logical(show_title),
    is.numeric(label_size), label_size > 0,
    is.numeric(title_size), title_size > 0,
    is.numeric(subtitle_size), subtitle_size > 0,
    is.numeric(axis_title_size), axis_title_size > 0,
    is.numeric(axis_text_size), axis_text_size > 0,
    is.numeric(legend_text_size), legend_text_size > 0,
    is.numeric(legend_title_size), legend_title_size > 0,
    is.logical(show_legend),
    legend_position %in% c("none", "right", "left", "top", "bottom")
  )
  
  structure(config, class = "aesthetic_config")
}


#' Plot Layout Configuration Class
#'
#' @description
#' S3 class for configuring plot layout including dimensions, margins,
#' grids and coordinate systems.
#'
#' @param width Plot width in inches
#' @param height Plot height in inches 
#' @param dpi Plot resolution
#' @param aspect_ratio Plot aspect ratio
#' @param show_grid Show plot grid
#' @param grid_type Grid type ("none", "major", "minor", "both")
#' @param grid_color Grid color
#' @param grid_linetype Grid line type
#' @param show_axis Show axes
#' @param axis_lines Show axis lines
#' @param plot_margin Plot margins in cm
#' @param coord_type Coordinate type ("fixed", "equal", "flip", "polar")
#' @param background_color Plot background color
#' @param panel_background_color Panel background color
#' @param panel_border Show panel border
#' @param panel_border_color Panel border color
#' @param save_format Plot save format ("png", "pdf", "svg", "eps")
#' @param reverse_x Numeric multiplier for x-axis direction (1 or -1)
#' @param reverse_y Numeric multiplier for y-axis direction (1 or -1)
#' @param x_limits Numeric vector of length 2 specifying c(min, max) for x-axis. If NULL, limits are set automatically.
#' @param y_limits Numeric vector of length 2 specifying c(min, max) for y-axis. If NULL, limits are set automatically.
#' @return A layout_config object
#' @export
new_layout_config <- function(
    width = 8,
    height = 8,
    dpi = 300,
    aspect_ratio = 1,
    show_grid = TRUE,
    grid_type = "major",
    grid_color = "grey80",
    grid_linetype = "dashed",
    show_axis = TRUE,
    axis_lines = TRUE,
    plot_margin = margin(1, 1, 1, 1, "cm"),
    coord_type = "fixed",
    background_color = "white",
    panel_background_color = "white",
    panel_border = TRUE,
    panel_border_color = "black",
    save_format = "png",
    reverse_x = 1,
    reverse_y = 1,
    x_limits = NULL,
    y_limits = NULL
) {
  config <- list(
    width = width,
    height = height,
    dpi = dpi,
    aspect_ratio = aspect_ratio,
    show_grid = show_grid,
    grid_type = grid_type,
    grid_color = grid_color,
    grid_linetype = grid_linetype,
    show_axis = show_axis,
    axis_lines = axis_lines,
    plot_margin = plot_margin,
    coord_type = coord_type,
    background_color = background_color,
    panel_background_color = panel_background_color,
    panel_border = panel_border,
    panel_border_color = panel_border_color,
    save_format = save_format,
    reverse_x = reverse_x,
    reverse_y = reverse_y,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  # Validate inputs
  stopifnot(
    is.numeric(width), width > 0,
    is.numeric(height), height > 0,
    is.numeric(dpi), dpi > 0,
    is.numeric(aspect_ratio), aspect_ratio > 0,
    is.logical(show_grid),
    grid_type %in% c("none", "major", "minor", "both"),
    is.character(grid_color),
    is.character(grid_linetype),
    is.logical(show_axis),
    is.logical(axis_lines),
    inherits(plot_margin, "margin"),
    coord_type %in% c("fixed", "equal", "flip", "polar"),
    is.character(background_color),
    is.character(panel_background_color),
    is.logical(panel_border),
    is.character(panel_border_color),
    save_format %in% c("png", "pdf", "svg", "eps"),
    reverse_x %in% c(1, -1),
    reverse_y %in% c(1, -1)
  )
  
  # Validate axis limits if provided
  if (!is.null(x_limits)) {
    if (!is.numeric(x_limits) || length(x_limits) != 2 || x_limits[1] >= x_limits[2]) {
      stop("x_limits must be a numeric vector of length 2 with min < max")
    }
  }
  if (!is.null(y_limits)) {
    if (!is.numeric(y_limits) || length(y_limits) != 2 || y_limits[1] >= y_limits[2]) {
      stop("y_limits must be a numeric vector of length 2 with min < max")
    }
  }
  
  structure(config, class = "layout_config")
}


#' Dimension Reduction Configuration Class
#' 
#' @description
#' S3 class for configuring dimension reduction parameters including
#' method selection and algorithm-specific parameters.
#'
#' @param method Dimension reduction method ("pca", "umap", "tsne")
#' @param n_components Number of components to compute
#' @param scale Scale the data before reduction
#' @param center Center the data before reduction
#' @param pca_params List of PCA-specific parameters
#' @param umap_params List of UMAP-specific parameters
#' @param tsne_params List of t-SNE-specific parameters
#' @param compute_loadings Compute and return loadings
#' @param random_state Random seed for reproducibility
#' @return A dim_reduction_config object
#' @export
new_dim_reduction_config <- function(
    method = "pca",
    n_components = 2,
    scale = FALSE,
    center = TRUE,
    pca_params = list(
      tol = sqrt(.Machine$double.eps),
      rank. = NULL
    ),
    umap_params = list(
      n_neighbors = 15,
      min_dist = 0.1,
      metric = "euclidean",
      n_epochs = 200
    ),
    tsne_params = list(
      perplexity = 30,
      max_iter = 1000,
      theta = 0.5
    ),
    compute_loadings = FALSE,
    random_state = NULL
) {
  config <- list(
    method = method,
    n_components = n_components,
    scale = scale,
    center = center,
    pca_params = pca_params,
    umap_params = umap_params,
    tsne_params = tsne_params,
    compute_loadings = compute_loadings,
    random_state = random_state
  )
  
  # Validate inputs
  stopifnot(
    method %in% c("pca", "umap", "tsne"),
    is.numeric(n_components), n_components > 0,
    is.logical(scale),
    is.logical(center),
    is.list(pca_params),
    is.list(umap_params),
    is.list(tsne_params),
    is.logical(compute_loadings),
    is.null(random_state) || is.numeric(random_state)
  )
  
  structure(config, class = "dim_reduction_config")
}


#' Validate Input Data Frame
#'
#' @description
#' Validates input data frame for visualization functions, checking required 
#' columns and data types.
#'
#' @param df Data frame to validate
#' @param ndim Number of dimensions expected in coordinate columns. Names of coordinate columns must start with a "V".
#' @param require_clusters Whether cluster column is required
#' @param require_temporal Whether year column is required
#' @return Validated data frame or throws error if invalid
#' @keywords internal
validate_topolow_df <- function(df, ndim, require_clusters = FALSE, require_temporal = FALSE) {
  # Check if df is a data frame
  if (!is.data.frame(df)) {
    stop("Input must be a data frame")
  }
  
  # Check coordinate columns
  coord_cols <- paste0("V", 1:ndim)
  if (!all(coord_cols %in% colnames(df))) {
    stop(sprintf("Missing coordinate columns. Required: %s", 
                 paste(coord_cols, collapse = ", ")))
  }
  
  # Check coordinate data types
  coord_data <- df[, coord_cols]
  if (!all(sapply(coord_data, is.numeric))) {
    stop("All coordinate columns must be numeric")
  }
  
  # Check required type columns
  if (!all(c("antigen", "antiserum") %in% colnames(df))) {
    stop("Missing antigen/antiserum indicator columns")
  }
  
  # Check cluster column if required
  if (require_clusters && !("cluster" %in% colnames(df))) {
    stop("Missing cluster column")
  }
  
  # Check temporal column if required
  if (require_temporal && !("year" %in% colnames(df))) {
    stop("Missing year column")
  }
  
  # Remove NA rows
  df <- na.omit(df)
  
  # Verify we still have data after NA removal
  if (nrow(df) == 0) {
    stop("No valid data remains after removing NA values")
  }
  
  return(df)
}


#' Perform Dimension Reduction
#'
#' @description
#' Applies configured dimension reduction method to input data.
#'
#' @param df Data frame containing coordinate data
#' @param config Dimension reduction configuration object
#' @return Data frame with reduced dimensions
#' @keywords internal
reduce_dimensions <- function(df, config) {
  if (!inherits(config, "dim_reduction_config")) {
    stop("config must be a dim_reduction_config object")
  }
  
  # Extract coordinate columns
  coord_cols <- grep("^V\\d+$", names(df), value = TRUE)
  coords <- df[, coord_cols]
  
  # Set random seed if provided
  if (!is.null(config$random_state)) {
    set.seed(config$random_state)
  }
  
  # Apply dimension reduction
  result <- switch(config$method,
                   "pca" = {
                     pca_result <- do.call(prcomp,
                                           c(list(x = coords, 
                                                  scale. = config$scale,
                                                  center = config$center),
                                             config$pca_params))
                     
                     # Scale components if requested
                     if (config$scale) {
                       # Calculate original distances
                       orig_dist <- as.matrix(dist(coords))
                       # Get first n_components
                       comps <- pca_result$x[, 1:config$n_components]
                       # Scale to match original distances
                       scaled_comps <- scale_to_original_distances(comps, orig_dist)
                       data.frame(scaled_comps)
                     } else {
                       data.frame(pca_result$x[, 1:config$n_components])
                     }
                   },
                   "umap" = {
                     umap_result <- do.call(umap,
                                            c(list(d = coords,
                                                   n_components = config$n_components),
                                              config$umap_params))
                     data.frame(umap_result$layout)
                   },
                   "tsne" = {
                     if (!requireNamespace("Rtsne", quietly = TRUE)) {
                       stop("Rtsne package required for t-SNE reduction")
                     }
                     tsne_result <- do.call(Rtsne::Rtsne,
                                            c(list(X = coords,
                                                   dims = config$n_components),
                                              config$tsne_params))
                     data.frame(tsne_result$Y)
                   },
                   stop("Unsupported dimension reduction method")
  )
  
  # Name columns consistently
  names(result) <- paste0("dim", 1:ncol(result))
  
  # Add original metadata columns
  metadata_cols <- setdiff(names(df), coord_cols)
  result[metadata_cols] <- df[metadata_cols]
  
  return(result)
}


#' Scale Reduced Dimensions to Match Original Distances
#'
#' @description
#' Helper function to scale reduced dimensions to better match original distances.
#'
#' @param reduced_coords Matrix of reduced coordinates
#' @param orig_dist Original distance matrix
#' @return Scaled coordinate matrix
#' @keywords internal
scale_to_original_distances <- function(reduced_coords, orig_dist) {
  # Calculate reduced distance matrix
  reduced_dist <- as.matrix(dist(reduced_coords))
  
  # Define optimization function for scaling factor
  sum_squared_diff <- function(scale_factor) {
    scaled_dist <- reduced_dist * scale_factor
    sum((orig_dist - scaled_dist)^2)
  }
  
  # Find optimal scaling factor
  optimal_scale <- optimize(sum_squared_diff, 
                            interval = c(0, 10))$minimum
  
  # Return scaled coordinates
  reduced_coords * optimal_scale
}


#' Create Base Theme
#'
#' @description
#' Creates a ggplot2 theme based on aesthetic and layout configurations.
#'
#' @param aesthetic_config Aesthetic configuration object
#' @param layout_config Layout configuration object
#' @return ggplot2 theme object
#' @keywords internal
create_base_theme <- function(aesthetic_config, layout_config) {
  theme_minimal() +
    theme(
      # Title 
      plot.title = if(aesthetic_config$show_title) 
        element_text(
          size = aesthetic_config$title_size,
          face = "bold",
          hjust = 0.5
        ) else element_blank(),
      
      # Text sizes  
      axis.title = element_text(size = aesthetic_config$axis_title_size),
      axis.text = element_text(size = aesthetic_config$axis_text_size),
      
      # Legend
      legend.position = if(aesthetic_config$show_legend) 
        aesthetic_config$legend_position else "none",
      legend.text = element_text(size = aesthetic_config$legend_text_size),
      legend.title = element_text(size = aesthetic_config$legend_title_size),
      
      # Grid
      panel.grid.major = if(layout_config$show_grid && 
                            layout_config$grid_type %in% c("major", "both"))
        element_line(color = layout_config$grid_color,
                     linetype = layout_config$grid_linetype) else element_blank(),
      panel.grid.minor = if(layout_config$show_grid && 
                            layout_config$grid_type %in% c("minor", "both"))
        element_line(color = layout_config$grid_color,
                     linetype = layout_config$grid_linetype) else element_blank(),
      
      # Backgrounds - Remove outer border by setting linewidth=0
      plot.background = element_rect(fill = layout_config$background_color, 
                                     color = NA, linewidth = 0),
      panel.background = element_rect(fill = layout_config$panel_background_color,
                                      color = NA, linewidth = 0),
      
      # Panel border and axis lines
      panel.border = if(layout_config$panel_border)
        element_rect(color = layout_config$panel_border_color,
                     fill = NA,
                     linewidth = 0.5) else element_blank(),
      axis.line = if(!layout_config$panel_border) 
        element_line(color = "black", size = 0.5) else element_blank(),
      
      # Margin
      plot.margin = layout_config$plot_margin
    )
}


#' Create Temporal Mapping Plot
#'
#' @description
#' Creates a visualization of points colored by time (year) using dimension reduction.
#' Points are colored on a gradient scale based on their temporal values, with
#' different shapes for antigens and antisera.
#' @param df Data frame containing:
#'        - V1, V2, ... Vn: Coordinate columns
#'        - antigen: Binary indicator for antigen points 
#'        - antiserum: Binary indicator for antiserum points
#'        - year: Numeric year values for temporal coloring
#' @param ndim Number of dimensions in input coordinates
#' @param dim_config Dimension reduction configuration object specifying method and parameters
#' @param aesthetic_config Aesthetic configuration object controlling plot appearance
#' @param layout_config Layout configuration object controlling plot dimensions and style.
#'        Use x_limits and y_limits in layout_config to set axis limits.
#' @param output_dir Character. Directory for output files. If NULL, uses current directory
#' 
#' @details
#' The function performs these steps:
#' 1. Validates input data structure and types
#' 2. Applies dimension reduction if ndim > 2
#' 3. Creates visualization with temporal color gradient
#' 4. Applies specified aesthetic and layout configurations
#' 5. Applies custom axis limits if specified in layout_config
#'
#' Different shapes distinguish between antigens and antisera points, while
#' color represents temporal progression.
#'
#' @return ggplot object containing the temporal mapping visualization
#'
#' @examples
#' \dontrun{
#' # Basic usage with default configurations
#' data <- data.frame(
#'   V1 = rnorm(100),
#'   V2 = rnorm(100),
#'   V3 = rnorm(100),
#'   antigen = rep(c(0,1), 50),
#'   antiserum = rep(c(1,0), 50),
#'   year = rep(2000:2009, each=10)
#' )
#' # Default axis limits
#' p1 <- plot_temporal_mapping(data, ndim=3)
#'
#' # Custom axis limits via layout configuration
#' layout_config <- new_layout_config(
#'   x_limits = c(-10, 10),
#'   y_limits = c(-8, 8)
#' )
#' p2 <- plot_temporal_mapping(data, ndim=3, 
#'                            layout_config=layout_config)
#' }
#'
#' @seealso 
#' \code{\link{plot_cluster_mapping}} for cluster-based visualization
#' \code{\link{plot_3d_mapping}} for 3D visualization
#' \code{\link{new_dim_reduction_config}} for dimension reduction options
#' \code{\link{new_aesthetic_config}} for aesthetic options
#' \code{\link{new_layout_config}} for layout options
#'
#' @export
plot_temporal_mapping <- function(df, ndim, 
                                  dim_config = new_dim_reduction_config(),
                                  aesthetic_config = new_aesthetic_config(),
                                  layout_config = new_layout_config(),
                                  output_dir = NULL) {
  # Validate input data
  df <- validate_topolow_df(df, ndim, require_temporal = TRUE)
  
  # Perform dimension reduction
  reduced_df <- reduce_dimensions(df, dim_config)
  
  # Apply axis reversals
  reduced_df$plot_x <- reduced_df$dim2 * layout_config$reverse_x
  reduced_df$plot_y <- reduced_df$dim1 * layout_config$reverse_y
  
  # Create base theme
  base_theme <- create_base_theme(aesthetic_config, layout_config)
  
  # Create point type with explicit factor levels
  reduced_df$point_type <- NA_character_  # Initialize
  reduced_df$point_type[reduced_df$antigen] <- "antigen"    # Use lowercase to match names
  reduced_df$point_type[reduced_df$antiserum] <- "antiserum"
  reduced_df$point_type <- factor(reduced_df$point_type, 
                                  levels = names(aesthetic_config$point_shapes))
  
  # Create plot
  p <- ggplot(reduced_df, aes(x = plot_x, y = plot_y, 
                              colour = year,
                              shape = point_type)) +
    geom_point(size = aesthetic_config$point_size,
               alpha = aesthetic_config$point_alpha) +
    scale_colour_gradient(low = aesthetic_config$gradient_colors$low,
                          high = aesthetic_config$gradient_colors$high,
                          na.value = "gray50") +
    scale_shape_manual(name = "Type",
                       values = aesthetic_config$point_shapes,
                       labels = c(antigen = "Antigen", 
                                  antiserum = "Antiserum")) +
    base_theme +
    labs(title = if(aesthetic_config$show_title) "Temporal Mapping" else "",
         x = "Dimension 1",
         y = "Dimension 2",
         colour = "Year")
  
  # Add axis limits if specified in layout_config
  if (!is.null(layout_config$x_limits)) {
    p <- p + scale_x_continuous(limits = layout_config$x_limits)
  }
  if (!is.null(layout_config$y_limits)) {
    p <- p + scale_y_continuous(limits = layout_config$y_limits)
  }
  
  # Add fixed coordinates if specified
  if(layout_config$coord_type == "fixed") {
    p <- p + coord_fixed(ratio = layout_config$aspect_ratio)
  }
  
  # Save plot if save format is specified
  if (!is.null(layout_config$save_format)) {
    filename <- sprintf("%s_temporal_mapping_ndim_%d.%s", 
                        dim_config$method, ndim, layout_config$save_format)
    save_plot(p, filename, layout_config, output_dir)
  }
  
  return(p)
}


#' Create Clustered Mapping Plots
#'
#' @description
#' Creates a visualization of points colored by cluster assignment using dimension 
#' reduction. Points are colored by cluster with different shapes for antigens and 
#' antisera.
#'
#' @param df_coords Data frame containing:
#'        - V1, V2, ... Vn: Coordinate columns
#'        - antigen: Binary indicator for antigen points 
#'        - antiserum: Binary indicator for antiserum points
#'        - cluster: Factor or integer cluster assignments
#' @param ndim Number of dimensions in input coordinates
#' @param dim_config Dimension reduction configuration object specifying method and parameters
#' @param aesthetic_config Aesthetic configuration object controlling plot appearance
#' @param layout_config Layout configuration object controlling plot dimensions and style.
#'        Use x_limits and y_limits in layout_config to set axis limits.
#' @param output_dir Character. Directory for output files. If NULL, uses current directory
#' 
#' @details
#' The function performs these steps:
#' 1. Validates input data structure and types
#' 2. Applies dimension reduction if ndim > 2
#' 3. Creates visualization with cluster-based coloring
#' 4. Applies specified aesthetic and layout configurations
#' 5. Applies custom axis limits if specified in layout_config
#'
#' Different shapes distinguish between antigens and antisera points, while
#' color represents cluster assignment. The color palette can be customized
#' through the aesthetic_config.
#'
#' @return ggplot object containing the cluster mapping visualization
#'
#' @examples
#' \dontrun{
#' # Basic usage with default configurations
#' data <- data.frame(
#'   V1 = rnorm(100),
#'   V2 = rnorm(100),
#'   V3 = rnorm(100),
#'   antigen = rep(c(0,1), 50),
#'   antiserum = rep(c(1,0), 50),
#'   cluster = rep(1:5, each=20)
#' )
#' p1 <- plot_cluster_mapping(data, ndim=3)
#'
#' # Custom configurations with specific color palette and axis limits
#' aesthetic_config <- new_aesthetic_config(
#'   point_size = 4,
#'   point_alpha = 0.7,
#'   color_palette = c("red", "blue", "green", "purple", "orange"),
#'   show_labels = TRUE,
#'   label_size = 3
#' )
#'
#' layout_config <- new_layout_config(
#'   width = 10,
#'   height = 8,
#'   coord_type = "fixed",
#'   show_grid = TRUE,
#'   grid_type = "major",
#'   x_limits = c(-10, 10),
#'   y_limits = c(-8, 8)
#' )
#'
#' p2 <- plot_cluster_mapping(
#'   data, 
#'   ndim = 3,
#'   aesthetic_config = aesthetic_config,
#'   layout_config = layout_config
#' )
#' }
#'
#' @seealso 
#' \code{\link{plot_temporal_mapping}} for temporal visualization
#' \code{\link{plot_3d_mapping}} for 3D visualization
#' \code{\link{plot_combined}} for creating multiple visualizations
#'
#' @export
plot_cluster_mapping <- function(df_coords, ndim,
                                 dim_config = new_dim_reduction_config(),
                                 aesthetic_config = new_aesthetic_config(),
                                 layout_config = new_layout_config(),
                                 output_dir = NULL) {
  # Validate input data
  df_coords <- validate_topolow_df(df_coords, ndim, require_clusters = TRUE)
  
  # Perform dimension reduction
  reduced_df <- reduce_dimensions(df_coords, dim_config)
  
  # Apply axis reversals
  reduced_df$plot_x <- reduced_df$dim2 * layout_config$reverse_x
  reduced_df$plot_y <- reduced_df$dim1 * layout_config$reverse_y
  
  # Create base theme
  base_theme <- create_base_theme(aesthetic_config, layout_config)
  
  # Get color palette
  n_clusters <- length(unique(reduced_df$cluster))
  colors <- aesthetic_config$color_palette[1:min(n_clusters, length(aesthetic_config$color_palette))]
  
  if (n_clusters > length(aesthetic_config$color_palette)) {
    warning("More clusters than available colors. Colors will be recycled.")
    colors <- rep_len(aesthetic_config$color_palette, n_clusters)
  }
  
  # Calculate optimal number of legend columns if legend is shown
  legend_theme <- if(aesthetic_config$show_legend) {
    n_legend_cols <- min(max(1, ceiling(n_clusters / 15)), 3)
    theme(
      legend.box = "vertical",
      legend.margin = margin(l = 5, r = 5),
      legend.spacing.y = unit(0.1, "cm"),
      legend.key.size = unit(0.8, "lines")
    )
  } else {
    theme()  # Empty theme if no legend
  }
  
  # Create point type with explicit factor levels
  reduced_df$point_type <- NA_character_  # Initialize
  reduced_df$point_type[reduced_df$antigen] <- "antigen"    # Use lowercase to match names
  reduced_df$point_type[reduced_df$antiserum] <- "antiserum"
  reduced_df$point_type <- factor(reduced_df$point_type, 
                                  levels = names(aesthetic_config$point_shapes))
  
  # Create plot
  p <- ggplot(reduced_df, aes(x = plot_x, y = plot_y,
                              colour = cluster,
                              shape = point_type)) +
    geom_point(size = aesthetic_config$point_size,
               alpha = aesthetic_config$point_alpha) +
    scale_colour_manual(values = colors) +
    scale_shape_manual(name = "Type",
                       values = aesthetic_config$point_shapes,
                       labels = c(antigen = "Antigen", 
                                  antiserum = "Antiserum")) +
    guides(colour = if(aesthetic_config$show_legend) 
      guide_legend(ncol = n_legend_cols,
                   title.position = "top",
                   byrow = TRUE)
      else "none") +
    base_theme +
    legend_theme +
    labs(title = if(aesthetic_config$show_title) "Cluster Mapping" else "",
         x = "Dimension 1", 
         y = "Dimension 2",
         colour = "Cluster")
  
  # Add axis limits if specified in layout_config
  if (!is.null(layout_config$x_limits)) {
    p <- p + scale_x_continuous(limits = layout_config$x_limits)
  }
  if (!is.null(layout_config$y_limits)) {
    p <- p + scale_y_continuous(limits = layout_config$y_limits)
  }
  
  # Add fixed coordinates if specified
  if(layout_config$coord_type == "fixed") {
    p <- p + coord_fixed(ratio = layout_config$aspect_ratio)
  }
  
  # Save plot if save format is specified
  if (!is.null(layout_config$save_format)) {
    filename <- sprintf("%s_cluster_mapping_ndim_%d.%s", 
                        dim_config$method, ndim, layout_config$save_format)
    save_plot(p, filename, layout_config, output_dir)
  }
  
  return(p)
}


#' Create 3D Visualization
#'
#' @description
#' Creates an interactive or static 3D visualization using rgl. Supports both 
#' temporal and cluster-based coloring schemes with configurable point 
#' appearances and viewing options.
#'
#' @param df Data frame containing:
#'        - V1, V2, ... Vn: Coordinate columns
#'        - antigen: Binary indicator for antigen points 
#'        - antiserum: Binary indicator for antiserum points
#'        - cluster: (Optional) Factor or integer cluster assignments
#'        - year: (Optional) Numeric year values for temporal coloring
#' @param ndim Number of dimensions in input coordinates (must be >= 3)
#' @param dim_config Dimension reduction configuration object
#' @param aesthetic_config Aesthetic configuration object
#' @param layout_config Layout configuration object
#' @param interactive Logical; whether to create an interactive plot
#' @param output_dir Character. Directory for output files. If NULL, uses current directory
#' 
#' @details
#' The function supports two main visualization modes:
#' 1. Interactive mode: Creates a manipulatable 3D plot window
#' 2. Static mode: Generates a static image from a fixed viewpoint
#'
#' Color schemes are automatically selected based on available data:
#' - If cluster data is present: Uses discrete colors per cluster
#' - If year data is present: Uses continuous color gradient
#' - Otherwise: Uses default point colors
#'
#' For data with more than 3 dimensions, dimension reduction is applied first.
#'
#' @return Invisibly returns rgl scene ID for further manipulation
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' data <- data.frame(
#'   V1 = rnorm(100),
#'   V2 = rnorm(100),
#'   V3 = rnorm(100),
#'   V4 = rnorm(100),
#'   antigen = rep(c(0,1), 50),
#'   antiserum = rep(c(1,0), 50),
#'   cluster = rep(1:5, each=20),
#'   year = rep(2000:2009, each=10)
#' )
#'
#' # Basic interactive plot
#' plot_3d_mapping(data, ndim=4)
#'
#' # Custom configuration for temporal visualization
#' aesthetic_config <- new_aesthetic_config(
#'   point_size = 5,
#'   point_alpha = 0.8,
#'   gradient_colors = list(
#'     low = "blue",
#'     high = "red"
#'   )
#' )
#'
#' layout_config <- new_layout_config(
#'   width = 12,
#'   height = 12,
#'   background_color = "black",
#'   show_axis = TRUE
#' )
#'
#' # Create customized static plot
#' plot_3d_mapping(data, ndim=4,
#'   aesthetic_config = aesthetic_config,
#'   layout_config = layout_config,
#'   interactive = FALSE
#' )
#'
#' # Dimension reduction with UMAP
#' dim_config <- new_dim_reduction_config(
#'   method = "umap",
#'   n_components = 3,
#'   umap_params = list(
#'     n_neighbors = 20,
#'     min_dist = 0.2
#'   )
#' )
#'
#' plot_3d_mapping(data, ndim=4,
#'   dim_config = dim_config,
#'   interactive = TRUE
#' )
#' }
#'
#' @seealso 
#' \code{\link{plot_temporal_mapping}} for 2D temporal visualization
#' \code{\link{plot_cluster_mapping}} for 2D cluster visualization
#' \code{\link{make_interactive}} for converting 2D plots to interactive versions
#'
#' @export
plot_3d_mapping <- function(df, ndim, 
                            dim_config = new_dim_reduction_config(),
                            aesthetic_config = new_aesthetic_config(),
                            layout_config = new_layout_config(),
                            interactive = TRUE,
                            output_dir = NULL) {
  # Validate input data and dimensions
  if(ndim < 3) {
    stop("3D visualization requires at least 3 dimensions in input data")
  }
  df <- validate_topolow_df(df, ndim)
  
  # Perform dimension reduction if needed
  if(ndim > 3) {
    dim_config$n_components <- 3
    reduced_df <- reduce_dimensions(df, dim_config)
  } else {
    reduced_df <- df
  }
  
  # Set up RGL window
  if(interactive) {
    open3d(windowRect = c(100, 100, 
                          100 + layout_config$width * 100,
                          100 + layout_config$height * 100))
  }
  
  # Set point colors
  if("cluster" %in% names(reduced_df)) {
    n_clusters <- length(unique(reduced_df$cluster))
    if(aesthetic_config$color_palette == "c25") {
      colors <- c25[1:n_clusters]
    } else {
      colors <- aesthetic_config$color_palette
    }
    point_colors <- colors[as.numeric(factor(reduced_df$cluster))]
  } else if("year" %in% names(reduced_df)) {
    # Create color gradient for years
    year_range <- range(reduced_df$year)
    year_colors <- colorRampPalette(c(aesthetic_config$gradient_colors$low,
                                      aesthetic_config$gradient_colors$high))
    year_normalized <- (reduced_df$year - year_range[1]) / diff(year_range)
    point_colors <- year_colors(100)[ceiling(year_normalized * 99) + 1]
  } else {
    point_colors <- "black"
  }
  
  # Plot points
  # Open a new 3D plotting device with a larger window size
  open3d(windowRect = c(100, 100, 1200, 1200))
  plot3d(reduced_df$dim1, reduced_df$dim2, reduced_df$dim3,
         col = point_colors,
         size = aesthetic_config$point_size * 0.8, # Adjust for 3D
         type = "s",  # spheres
         alpha = aesthetic_config$point_alpha)
  
  # Add axes and labels
  if(layout_config$show_axis) {
    axes3d(edges = "bbox",
           labels = TRUE,
           tick = TRUE,
           box = TRUE)
  }
  
  # Save if not interactive
  if(!interactive) {
    if (is.null(output_dir)) {
      output_dir <- getwd()
    }
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    filename <- sprintf("3d_mapping_ndim_%d.%s", 
                        ndim, layout_config$save_format)
    full_path <- file.path(output_dir, filename)
    
    rgl.snapshot(filename = full_path)
    close3d()
    return(invisible(full_path))
  }
  
  # Return scene ID invisibly
  invisible(rgl.cur())
}


#' Save Plot to File
#'
#' @description
#' Saves a plot (ggplot or rgl scene) to file with specified configuration.
#' Supports multiple output formats and configurable dimensions.
#'
#' @param plot ggplot or rgl scene object to save
#' @param filename Output filename (with or without extension)
#' @param layout_config Layout configuration object controlling output parameters
#' @param output_dir Character. Directory for output files. If NULL, uses current directory
#' 
#' @return Invisible NULL
#' 
#' @details
#' Supported file formats:
#' - PNG: Best for web and general use
#' - PDF: Best for publication quality vector graphics
#' - SVG: Best for web vector graphics
#' - EPS: Best for publication quality vector graphics
#'
#' The function will:
#' 1. Auto-detect plot type (ggplot or rgl)
#' 2. Use appropriate saving method
#' 3. Apply layout configuration settings
#' 4. Add file extension if not provided
#'
#' @examples
#' \dontrun{
#' # Create sample plot
#' data <- data.frame(
#'   V1 = rnorm(100),
#'   V2 = rnorm(100),
#'   antigen = rep(c(0,1), 50),
#'   antiserum = rep(c(1,0), 50),
#'   year = rep(2000:2009, each=10)
#' )
#' p <- plot_temporal_mapping(data, ndim=2)
#'
#' # Basic save
#' save_plot(p, "temporal_plot.png")
#'
#' # Save with custom layout
#' layout_config <- new_layout_config(
#'   width = 12,
#'   height = 8,
#'   dpi = 600,
#'   save_format = "pdf"
#' )
#'
#' save_plot(p, "high_res_plot", layout_config)
#'
#' # Save 3D plot
#' p3d <- plot_3d_mapping(data, ndim=3, interactive=FALSE)
#' save_plot(p3d, "3d_plot.png", layout_config)
#' }
#'
#' @export
save_plot <- function(plot, filename, layout_config = new_layout_config(),
                      output_dir = NULL) {
  # Handle output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Extract file extension
  ext <- tools::file_ext(filename)
  if(ext == "") {
    filename <- paste0(filename, ".", layout_config$save_format)
    ext <- layout_config$save_format
  }
  
  # Validate format
  if(!ext %in% c("png", "pdf", "svg", "eps")) {
    stop("Unsupported file format: ", ext)
  }
  
  # Create full file path
  full_path <- file.path(output_dir, filename)
  
  # Save based on plot type
  if(inherits(plot, "ggplot")) {
    ggsave(filename = full_path,
           plot = plot,
           width = layout_config$width,
           height = layout_config$height,
           dpi = layout_config$dpi)
  } else if(inherits(plot, c("arrangement", "grob"))) {
    # For arranged plots from gridExtra
    ggsave(filename = full_path,
           plot = plot,
           width = layout_config$width,
           height = layout_config$height,
           dpi = layout_config$dpi)
  } else if(is.character(plot) && file.exists(plot)) {
    # For 3D plots, copy the temporary file
    file.copy(plot, full_path, overwrite = TRUE)
    file.remove(plot)  # Clean up temporary file
  } else {
    stop("Unsupported plot type")
  }
  
  invisible(NULL)
}


#' Create Interactive Plot
#'
#' @description
#' Converts a static ggplot visualization to an interactive plotly visualization
#' with customizable tooltips and interactive features.
#'
#' @param plot ggplot object to convert
#' @param tooltip_vars Vector of variable names to include in tooltips
#'
#' @details
#' The function enhances static plots by adding:
#' - Hover tooltips with data values
#' - Zoom capabilities
#' - Pan capabilities
#' - Click interactions
#' - Double-click to reset
#'
#' If tooltip_vars is NULL, the function attempts to automatically determine
#' relevant variables from the plot's mapping.
#'
#' @return plotly object with interactive features
#'
#' @examples
#' \dontrun{
#' # Create sample data and plot
#' data <- data.frame(
#'   V1 = rnorm(100),
#'   V2 = rnorm(100),
#'   antigen = rep(c(0,1), 50),
#'   antiserum = rep(c(1,0), 50),
#'   year = rep(2000:2009, each=10),
#'   cluster = rep(1:5, each=20)
#' )
#'
#' # Create temporal plot
#' p1 <- plot_temporal_mapping(data, ndim=2)
#'
#' # Make interactive with default tooltips
#' p1_interactive <- make_interactive(p1)
#'
#' # Create cluster plot with custom tooltips
#' p2 <- plot_cluster_mapping(data, ndim=2)
#' p2_interactive <- make_interactive(p2,
#'   tooltip_vars = c("cluster", "year", "antigen")
#' )
#' }
#'
#' @export
make_interactive <- function(plot, tooltip_vars = NULL) {
  if(!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly package is required for interactive plots")
  }
  
  if(is.null(tooltip_vars)) {
    # Try to determine relevant variables from the plot
    mapping <- plot$mapping
    tooltip_vars <- sapply(mapping, function(x) as.character(x)[2])
    tooltip_vars <- tooltip_vars[!tooltip_vars %in% c("x", "y")]
  }
  
  # Ensure 'text' is included for the name tooltips
  if(!"text" %in% tooltip_vars) {
    tooltip_vars <- c(tooltip_vars, "text")
  }
  
  p <- plotly::ggplotly(plot, tooltip = tooltip_vars)
  
  # Configure layout
  p <- plotly::layout(p,
                      hoverlabel = list(
                        bgcolor = "white",
                        font = list(family = "Arial", size = 12)
                      ))
  
  return(p)
}


#' Create Combined Visualization
#'
#' @description
#' Creates multiple coordinated visualizations of the same data using different 
#' methods and arrangements. Supports combining temporal, cluster, and 3D 
#' visualizations in flexible layouts.
#'
#' @param df_coords Data frame containing:
#'        - V1, V2, ... Vn: Coordinate columns
#'        - antigen: Binary indicator for antigen points 
#'        - antiserum: Binary indicator for antiserum points
#'        - cluster: (Optional) Factor or integer cluster assignments
#'        - year: (Optional) Numeric year values for temporal coloring
#' @param ndim Number of dimensions in input coordinates
#' @param plot_types Vector of plot types to create ("temporal", "cluster", "3d")
#' @param dim_config Dimension reduction configuration object
#' @param aesthetic_config Aesthetic configuration object
#' @param layout_config Layout configuration object
#' @param arrange How to arrange multiple plots ("grid", "vertical", "horizontal")
#' @param output_dir Character. Directory for output files. If NULL, uses current directory
#' 
#' @details
#' This function provides a high-level interface for creating multiple coordinated
#' views of the same data. It supports:
#' 
#' Plot Types:
#' - temporal: Time-based color gradients
#' - cluster: Cluster-based discrete colors
#' - 3d: Three-dimensional interactive or static views
#'
#' Arrangement Options:
#' - grid: Automatic square-like arrangement
#' - vertical: Plots stacked vertically
#' - horizontal: Plots arranged horizontally
#'
#' All plots share consistent:
#' - Color schemes
#' - Point styles
#' - Axis scales
#' - Theme elements
#'
#' @return Combined plot object (grid arrangement of plots)
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' data <- data.frame(
#'   V1 = rnorm(100),
#'   V2 = rnorm(100),
#'   V3 = rnorm(100),
#'   V4 = rnorm(100),
#'   antigen = rep(c(0,1), 50),
#'   antiserum = rep(c(1,0), 50),
#'   cluster = rep(1:5, each=20),
#'   year = rep(2000:2009, each=10)
#' )
#'
#' # Basic combined plot
#' p1 <- plot_combined(data, ndim=4,
#'   plot_types = c("temporal", "cluster")
#' )
#'
#' # Advanced configuration
#' dim_config <- new_dim_reduction_config(
#'   method = "umap",
#'   n_components = 2,
#'   scale = TRUE,
#'   umap_params = list(
#'     n_neighbors = 15,
#'     min_dist = 0.1
#'   )
#' )
#'
#' aesthetic_config <- new_aesthetic_config(
#'   point_size = 3,
#'   point_alpha = 0.7,
#'   point_shapes = c(antigen = 17, antiserum = 1),
#'   gradient_colors = list(
#'     low = "navy",
#'     high = "red"
#'   ),
#'   show_labels = TRUE,
#'   label_size = 3
#' )
#'
#' layout_config <- new_layout_config(
#'   width = 12,
#'   height = 8,
#'   aspect_ratio = 1,
#'   show_grid = TRUE,
#'   grid_type = "major",
#'   background_color = "white",
#'   panel_border = TRUE
#' )
#'
#' # Create comprehensive visualization
#' p2 <- plot_combined(data, ndim=4,
#'   plot_types = c("temporal", "cluster", "3d"),
#'   dim_config = dim_config,
#'   aesthetic_config = aesthetic_config,
#'   layout_config = layout_config,
#'   arrange = "grid"
#' )
#'
#' # Save combined plot
#' save_plot(p2, "combined_visualization.pdf")
#'
#' # Create interactive versions
#' p3 <- plot_combined(data, ndim=4,
#'   plot_types = c("temporal", "cluster"),
#'   arrange = "horizontal"
#' )
#'
#' p3_interactive <- make_interactive(p3,
#'   tooltip_vars = c("year", "cluster", "antigen")
#' )
#'
#' # Example with different layouts
#' # Vertical arrangement
#' p4 <- plot_combined(data, ndim=4,
#'   plot_types = c("temporal", "cluster", "3d"),
#'   arrange = "vertical"
#' )
#'
#' # Horizontal arrangement with temporal and cluster only
#' p5 <- plot_combined(data, ndim=4,
#'   plot_types = c("temporal", "cluster"),
#'   arrange = "horizontal"
#' )
#'
#' # Grid arrangement with custom layout
#' layout_config$width <- 15
#' layout_config$height <- 15
#' p6 <- plot_combined(data, ndim=4,
#'   plot_types = c("temporal", "cluster", "3d"),
#'   layout_config = layout_config,
#'   arrange = "grid"
#' )
#'
#' # Example workflow for publication-quality figures
#' # 1. Create base visualization
#' p7 <- plot_combined(data, ndim=4,
#'   plot_types = c("temporal", "cluster")
#' )
#'
#' # 2. Customize for publication
#' layout_config <- new_layout_config(
#'   width = 8,
#'   height = 6,
#'   dpi = 600,
#'   save_format = "pdf",
#'   background_color = "white",
#'   panel_border = TRUE,
#'   grid_type = "major"
#' )
#'
#' # 3. Save high-resolution version
#' save_plot(p7, "publication_figure.pdf", layout_config)
#' }
#'
#' @seealso 
#' \code{\link{plot_temporal_mapping}} for individual temporal plots
#' \code{\link{plot_cluster_mapping}} for individual cluster plots
#' \code{\link{plot_3d_mapping}} for individual 3D plots
#' \code{\link{make_interactive}} for creating interactive versions
#' \code{\link{save_plot}} for saving plots to files
#'
#' @export
plot_combined <- function(df_coords, ndim,
                          plot_types = c("temporal", "cluster"),
                          dim_config = new_dim_reduction_config(),
                          aesthetic_config = new_aesthetic_config(),
                          layout_config = new_layout_config(),
                          arrange = "grid",
                          output_dir = NULL) {
  
  # Create requested plots
  plots <- list()
  
  if("temporal" %in% plot_types) {
    plots$temporal <- plot_temporal_mapping(df_coords, ndim, dim_config,
                                            aesthetic_config, layout_config)
  }
  
  if("cluster" %in% plot_types) {
    plots$cluster <- plot_cluster_mapping(df_coords, ndim, dim_config,
                                          aesthetic_config, layout_config)
  }
  
  if("3d" %in% plot_types && ndim >= 3) {
    plots$three_d <- plot_3d_mapping(df_coords, ndim, dim_config,
                                     aesthetic_config, layout_config,
                                     interactive = FALSE)
  }
  
  # Arrange plots
  if(length(plots) == 1) {
    return(plots[[1]])
  } else {
    if(!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("gridExtra package is required for combined plots")
    }
    
    n_plots <- length(plots)
    
    # Calculate grid dimensions based on arrangement type
    switch(arrange,
           "grid" = {
             ncol <- ceiling(sqrt(n_plots))
             nrow <- ceiling(n_plots/ncol)
           },
           "vertical" = {
             nrow <- n_plots
             ncol <- 1
           },
           "horizontal" = {
             nrow <- 1
             ncol <- n_plots
           },
           stop("Invalid arrange value"))
    
    # Adjust layout_config dimensions for the arrangement
    if(arrange == "vertical") {
      plot_heights <- rep(layout_config$height/n_plots, n_plots)
      plot_widths <- layout_config$width
    } else if(arrange == "horizontal") {
      plot_heights <- layout_config$height
      plot_widths <- rep(layout_config$width/n_plots, n_plots)
    } else { # grid
      plot_heights <- rep(layout_config$height/nrow, nrow)
      plot_widths <- rep(layout_config$width/ncol, ncol)
    }
    
    combined <- do.call(gridExtra::grid.arrange,
                        c(plots,
                          list(nrow = nrow,
                               ncol = ncol,
                               widths = plot_widths,
                               heights = plot_heights)))
    
    # Save combined plot if requested
    if (!is.null(layout_config$save_format)) {
      filename <- sprintf("combined_mapping_ndim_%d.%s", 
                          ndim, layout_config$save_format)
      save_plot(combined, filename, layout_config, output_dir)
    }
    
    return(combined)
  }
}


#' Create Diagnostic Plots for Multiple Chains
#'
#' @description
#' Creates trace and density plots for multiple Adaptive Monte Carlo Sampling or optimization chains to assess
#' convergence and mixing. Displays parameter trajectories and distributions across chains.
#'
#' @param chain_files Character vector of paths to CSV files containing chain data
#' @param mutual_size Integer number of samples to use from end of each chain
#' @param output_file Character path for saving plot
#' @param output_dir Character. Directory for output files. If NULL, uses current directory
#' @param save_plot Logical. Whether to save plots to files. Default: TRUE
#' @param width,height,res Plot dimensions and resolution for saving
#' @return Invisible NULL, saves plot to file
#' @examples
#' \dontrun{
#' chain_files <- c("chain1.csv", "chain2.csv", "chain3.csv")
#' create_diagnostic_plots(chain_files, mutual_size = 2000,
#'   output_file = "chain_diagnostics.png")
#' }
#' @export
create_diagnostic_plots <- function(chain_files, 
                                  mutual_size = 2000,
                                  output_file = "diagnostic_plots.png",
                                  output_dir = NULL,
                                  save_plot = TRUE,
                                  width = 3000, height = 3000, res = 300) {
  
  # Read chains
  chains <- lapply(chain_files, read.csv)
  
  # Process samples 
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  for(i in seq_along(chains)) {
    chains[[i]] <- chains[[i]][
      (nrow(chains[[i]]) - mutual_size - 1):nrow(chains[[i]]), 
      par_names]
  }
  
  n_params <- ncol(chains[[1]])
  n_chains <- length(chains)
  
  # Create plots for each parameter
  plot_list <- vector("list", n_params * 2)
  for (i in seq_len(n_params)) {
    # Prepare trace data
    trace_data <- do.call(rbind, lapply(seq_len(n_chains), function(j) {
      data.frame(
        Chain = j, 
        Iteration = seq_len(nrow(chains[[j]])),
        Value = chains[[j]][,i]
      )
    }))
    
    # Trace plot
    plot_list[[i*2-1]] <- ggplot(trace_data, 
      aes(x = Iteration, y = Value, color = factor(Chain))) +
      geom_line(size = 0.5) +
      labs(title = paste("Trace Plot:", par_names[i]),
           x = "Iteration", y = "Value") +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
      )
    
    # Density plot  
    plot_list[[i*2]] <- ggplot(trace_data,
      aes(x = Value, color = factor(Chain))) +
      geom_density(alpha = 0.3) +
      labs(title = paste("Density:", par_names[i]),
           x = "Value", y = "Density") +
      theme_minimal() +
      theme(
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
      )
  }
  
  # Arrange and save
  combined_plot <- gridExtra::grid.arrange(
    grobs = plot_list, 
    ncol = 2
  )
  
  # Optionally save
  if(save_plot) {
    if (is.null(output_dir)) {
      output_dir <- getwd()
    }
    output_file <- file.path(output_dir, output_file)
    
    ggsave(output_file, combined_plot,
           width = width/300, height = height/300, 
           dpi = res, limitsize = FALSE)
  }
  
  # Return plot object
  return(combined_plot)
}


#' Create Profile Likelihood Plot (Legacy Version)
#'
#' @description
#' Creates a visualization of profile likelihood for a parameter showing maximum
#' likelihood estimates and confidence intervals. For legacy data formats.
#' Consider using the S3 method plot.profile_likelihood() instead.
#'
#' @param LL_list_param Data frame with parameter values and log-likelihoods
#' @param param_name Character name of parameter being profiled
#' @param LL_max Numeric maximum log-likelihood value
#' @return A ggplot object
#' @examples
#' \dontrun{
#' LL_data <- data.frame(
#'   param = seq(0, 1, 0.1),
#'   LL = dnorm(seq(0, 1, 0.1), 0.5, 0.2)
#' )
#' plot_profile_likelihood(LL_data, "mu", max(LL_data$LL))
#' }
#' @export 
plot_profile_likelihood <- function(LL_list_param, param_name, LL_max) {
  CI_95_LL <- LL_max - 3.84 / 2
  ggplot(LL_list_param, aes_string(x = "param", y = "LL")) +
    geom_line(color = "steelblue", size = 0.5) +
    geom_hline(yintercept = CI_95_LL, 
               linetype = "dashed", color = "black", size = 0.4) +
    geom_text(aes(x = min(param), y = CI_95_LL + 0.02, 
                  label = "95% CI"),
              color = "black", vjust = -0.5, hjust = -0.05) +
    labs(title = paste("Profile Likelihood:", param_name),
         x = param_name,
         y = "Log Likelihood") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    scale_y_continuous(labels = scales::comma)
}


#' Plot Convergence Analysis Results
#'
#' @description
#' Visualizes convergence diagnostics including parameter mean trajectories
#' and covariance changes over iterations. Covariance norm changes measured by 
#' Frobenius norm (also called Hilbert-Schmidt norm), the square root of the 
#' sum of the absolute squares of all matrix elements = sqrt(sum|a_ij|²)
#'
#' @param conv_results List output from check_gaussian_convergence()
#' @param param_names Character vector of parameter names
#' @return A grid of plots showing convergence metrics
#' @examples
#' \dontrun{
#' results <- check_gaussian_convergence(chain_data)
#' plot_convergence_analysis(results, c("mu", "sigma"))
#' }
#' @export
plot_convergence_analysis <- function(conv_results, param_names) {
  # Parameter mean plots
  mean_plots <- lapply(seq_along(param_names), function(i) {
    plot_data <- data.frame(
      Iteration = seq_len(nrow(conv_results$mean_history)),
      Value = conv_results$mean_history[,i]
    )
    
    ggplot(plot_data, aes(x = Iteration, y = Value)) +
      geom_line(color = "steelblue") +
      labs(title = paste("Parameter Mean:", param_names[i]),
           x = "Iteration", y = "Value") +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
      )
  })
  
  # Covariance change plot
  cov_plot <- {
    plot_data <- data.frame(
      Iteration = seq_along(conv_results$cov_changes),
      Change = conv_results$cov_changes
    )
    
    ggplot(plot_data, aes(x = Iteration, y = Change)) +
      geom_line(color = "steelblue") +
      labs(title = "Covariance Changes",
           x = "Iteration", y = "Relative Change") +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
      )
  }
  
  # Combine plots
  gridExtra::grid.arrange(
    grobs = c(mean_plots, list(cov_plot)),
    ncol = 2
  )
}


#' Plot Method for Convergence Diagnostics
#' 
#' @description
#' Plots convergence diagnostics including parameter mean trajectories and 
#' covariance changes over iterations.
#'
#' @param x A topolow_convergence object from check_gaussian_convergence()
#' @param ... Additional arguments passed to underlying plot functions
#' @return A grid of plots showing convergence metrics
#' @method plot topolow_convergence
#' @export
plot.topolow_convergence <- function(x, ...) {
  plot_convergence_analysis(x, x$param_names)
}

#' Print Method for Convergence Diagnostics
#'
#' @param x A topolow_convergence object
#' @param ... Additional arguments passed to print
#' @method print topolow_convergence
#' @export
print.topolow_convergence <- function(x, ...) {
  cat("TopoLow Convergence Diagnostics:\n")
  cat(sprintf("Convergence achieved: %s\n", x$converged))
  cat(sprintf("Mean convergence: %s\n", x$mean_converged))
  cat(sprintf("Covariance convergence: %s\n", x$cov_converged))
  cat("\nFinal parameter means:\n")
  print(setNames(x$final_mean, x$param_names))
}

#' Plot Method for Adaptive Monte Carlo Sampling Diagnostics
#'
#' @description
#' Creates trace and density plots for multiple chains to assess convergence
#' and mixing.
#'
#' @param x A topolow_amcs_diagnostics object
#' @param output_file Character path for saving plot 
#' @param width,height,res Plot dimensions and resolution
#' @param ... Additional arguments passed to plot functions
#' @return Invisible NULL, saves plot to file
#' @method plot topolow_amcs_diagnostics
#' @export
plot.topolow_amcs_diagnostics <- function(x, 
                                        output_file = "mc_diagnostics.png",
                                        width = 3000, height = 3000, 
                                        res = 300, ...) {
  create_diagnostic_plots(x$chains, x$mutual_size, 
                         output_file = output_file,
                         width = width, height = height, 
                         res = res, ...)
}

#' Print Method for Adaptive Monte Carlo Sampling Diagnostics
#'
#' @param x A topolow_amcs_diagnostics object
#' @param ... Additional arguments passed to print
#' @method print topolow_amcs_diagnostics
#' @export
print.topolow_amcs_diagnostics <- function(x, ...) {
  cat("TopoLow Adaptive Monte Carlo Sampling Diagnostics:\n")
  cat("\nR-hat values (should be close to 1):\n")
  print(setNames(x$rhat, x$param_names))
  cat("\nEffective Sample Sizes:\n")
  print(setNames(x$ess, x$param_names))
}


#' Plot Network Structure Analysis
#'
#' @description
#' Creates visualization of distance matrix network structure showing data
#' availability patterns and connectivity.
#'
#' @param network_results List output from analyze_network_structure()
#' @param scenario_name Character string for output file naming
#' @param aesthetic_config Plot aesthetic configuration object
#' @param layout_config Plot layout configuration object
#' @return ggplot object
#' @examples
#' \dontrun{
#' net_analysis <- analyze_network_structure(dist_mat)
#' p <- plot_network_structure(net_analysis, "scenario1")
#' }
#' @importFrom igraph graph_from_adjacency_matrix layout_with_fr get.edgelist
#' @importFrom ggplot2 coord_fixed geom_segment

#' @export
plot_network_structure <- function(network_results, scenario_name,
                                 aesthetic_config = new_aesthetic_config(),
                                 layout_config = new_layout_config()) {
  # Create graph layout
  graph <- igraph::graph_from_adjacency_matrix(
    network_results$adjacency,
    mode = "undirected"
  )
  
  layout <- igraph::layout_with_fr(graph)
  layout_df <- data.frame(
    x = layout[,1],
    y = layout[,2],
    node = rownames(network_results$adjacency)
  )
  
  # Create edge data
  edges <- as.data.frame(igraph::get.edgelist(graph))
  edges <- merge(
    edges,
    layout_df,
    by.x = "V1",
    by.y = "node"
  )
  edges <- merge(
    edges,
    layout_df,
    by.x = "V2",
    by.y = "node",
    suffixes = c(".from", ".to")
  )
  
  # Create plot
  plot <- ggplot() +
    geom_segment(
      data = edges,
      aes(
        x = x.from,
        y = y.from,
        xend = x.to,
        yend = y.to
      ),
      alpha = 0.3,
      color = "grey50"
    ) +
    geom_point(
      data = layout_df,
      aes(x = x, y = y),
      color = "red",
      alpha = aesthetic_config$point_alpha
    ) +
    coord_fixed() +
    theme_void() +
    theme(
      plot.title = element_text(
        size = aesthetic_config$title_size,
        hjust = 0.5
      ),
      legend.position = if(aesthetic_config$show_legend) 
        aesthetic_config$legend_position else "none"
    ) +
    labs(
      title = sprintf(
        "Network Structure (%.1f%% Complete)",
        network_results$summary$completeness * 100
      ),
      size = "Connections"
    )
  
  # Save plot
  filename <- sprintf(
    "%s_network.%s",
    scenario_name,
    layout_config$save_format
  )
  
  ggsave(
    filename = filename,
    plot = plot,
    width = layout_config$width,
    height = layout_config$height,
    dpi = layout_config$dpi
  )
  
  return(plot)
}



#' Plot Distance Matrix Heatmap
#'
#' @description
#' Creates heatmap visualization of distance matrix showing patterns and
#' structure in the measurements.
#'
#' @param heatmap_data List output from prepare_heatmap_data()
#' @param scenario_name Character string for output file naming
#' @param aesthetic_config Plot aesthetic configuration object
#' @param layout_config Plot layout configuration object
#'
#' @return A ggplot object containing:
#'   - Heatmap visualization of the distance matrix
#'   - Color gradient representing distance values
#'   - Title showing matrix completeness percentage
#'
#' @examples
#' \dontrun{
#' # Basic heatmap
#' hmap_data <- prepare_heatmap_data(dist_mat)
#' p <- plot_distance_heatmap(hmap_data, "scenario1")
#'
#' # Heatmap with clustering
#' hmap_data <- prepare_heatmap_data(dist_mat, cluster_rows = TRUE)
#' p2 <- plot_distance_heatmap(hmap_data, "scenario1")
#'
#' # Custom configuration
#' aesthetic_config <- new_aesthetic_config(
#'   gradient_colors = list(low = "navy", high = "red")
#' )
#' p3 <- plot_distance_heatmap(hmap_data, "scenario1",
#'                            aesthetic_config = aesthetic_config)
#' }
#'
#' @importFrom reshape2 melt
#' @export
plot_distance_heatmap <- function(heatmap_data, scenario_name,
                                  aesthetic_config = new_aesthetic_config(),
                                  layout_config = new_layout_config()) {
  
  # Prepare data for plotting
  plot_data <- reshape2::melt(heatmap_data$matrix_data)
  plot_data$value <- as.numeric(as.character(plot_data$value))
  
  # Remove NAs and non-numeric values
  plot_data <- plot_data[!is.na(plot_data$value) & is.finite(plot_data$value),]
  
  # Create base theme  
  base_theme <- theme_minimal() +
    theme(
      plot.title = element_text(
        size = aesthetic_config$title_size,
        hjust = 0.5
      ),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = aesthetic_config$axis_text_size * 0.7
      ),
      axis.text.y = element_text(
        size = aesthetic_config$axis_text_size * 0.7
      ),
      plot.margin = layout_config$plot_margin
    )
  
  # Create plot
  plot <- ggplot(
    plot_data,
    aes(Var1, Var2, fill = value)
  ) +
    geom_tile() +
    scale_fill_gradient(
      low = aesthetic_config$gradient_colors$low,
      high = aesthetic_config$gradient_colors$high,
      na.value = "grey"
    ) +
    base_theme +
    labs(
      x = NULL,
      y = NULL,
      title = sprintf(
        "Distance Matrix (%.1f%% Complete)",
        heatmap_data$stats$completeness * 100
      )
    )
  
  # Save plot if required
  if (!is.null(scenario_name)) {
    filename <- sprintf(
      "%s_heatmap.%s",
      scenario_name,
      layout_config$save_format
    )
    
    ggsave(
      filename = filename,
      plot = plot,
      width = layout_config$width,
      height = layout_config$height,
      dpi = layout_config$dpi
    )
  }
  
  return(plot)
}


#' Plot Fitted vs True Distances
#'
#' @description
#' Creates diagnostic plots comparing fitted distances from a model against true distances.
#' Generates both a scatter plot with prediction intervals and a residuals plot.
#'
#' @param distance_matrix Matrix of true distances
#' @param p_dist_mat Matrix of predicted/fitted distances
#' @param scenario_name Character string for output file naming
#' @param ndim Integer number of dimensions used in the model
#' @param confidence_level Numeric confidence level for prediction intervals (default: 0.95)
#' @param save_plot Logical. Whether to save plots to files. Default: TRUE
#' @param output_dir Character. Directory for output files. If NULL, uses current directory
#' @return Invisibly returns NULL, creates two plot files:
#'   - \{scenario_name\}_prediction_scatter_dim_\{ndim\}.png
#'   - \{scenario_name\}_residuals_vs_fitted_dim_\{ndim\}.png
#' @examples
#' \dontrun{
#' # Create scatter and residual plots
#' scatterplot_fitted_vs_true(truth_matrix, predicted_matrix, 
#'                           scenario_name = "example",
#'                           ndim = 5)
#' }
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_smooth labs theme_classic
#' @importFrom scales comma
#' @export
# In visualization.R
scatterplot_fitted_vs_true <- function(distance_matrix, p_dist_mat, 
                                       scenario_name = NA, ndim = NA,
                                       save_plot = TRUE,
                                       output_dir = NULL,
                                       confidence_level = 0.95) {
  # Create evaluation data frame with simple numeric conversion to remove NA and thresholded values
  evaldf <- data.frame(
    distance_matrix = suppressWarnings(as.numeric(as.vector(distance_matrix))),
    p_dist_mat = as.vector(p_dist_mat)
  )
  evaldf <- na.omit(evaldf)
  
  # Calculate performance metrics  
  mae <- mean(abs(evaldf$distance_matrix - evaldf$p_dist_mat), na.rm = TRUE)
  Pearson_corr <- cor(evaldf$p_dist_mat, evaldf$distance_matrix, 
                      method = "pearson", use = "pairwise.complete.obs")
  
  # Calculate prediction interval
  margin_of_error <- calculate_prediction_interval(distance_matrix, p_dist_mat,
                                                 confidence_level = confidence_level)
  
  # Get regression line parameters
  reg_line <- lm(p_dist_mat ~ distance_matrix, data = evaldf)
  slope <- coef(reg_line)[2]
  intercept <- coef(reg_line)[1]
  
  # Create scatter plot 
  scatter_plot <- ggplot(evaldf, aes(x = distance_matrix, y = p_dist_mat)) +
    geom_point(alpha = 0.6, color = "#3366FF", size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "green", size = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "red", 
                linetype = "dashed", level = 0.999,
                fill = "pink", alpha = 0.6, size = 0.3) +
    geom_abline(intercept = intercept + mean(margin_of_error), 
                slope = slope, linetype = "dashed", 
                color = "black", size = 0.7) +
    geom_abline(intercept = intercept - mean(margin_of_error), 
                slope = slope, linetype = "dashed", 
                color = "black", size = 0.7) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    scale_x_continuous(limits = c(0, max(evaldf$distance_matrix)), 
                      expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, max(evaldf$p_dist_mat)), 
                      expand = c(0, 0)) +
    labs(
      x = "Original Distance",
      y = "Estimated Distance",
      title = sprintf("Mean Absolute Error = %.2f, Pearson's r = %.2f\nPrediction interval = ±%.2f",
                     mae, Pearson_corr, margin_of_error)
    )
  
  # Create residuals plot
  evaldf$residuals <- evaldf$distance_matrix - evaldf$p_dist_mat
  
  residuals_plot <- ggplot(evaldf, aes(x = p_dist_mat, y = residuals)) +
    geom_point(alpha = 0.6, color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_classic() +
    theme(plot.title = element_text(size = 10)) +
    labs(x = "Fitted Values", y = "Residuals",
         title = "Residuals vs. Fitted Values") +
    annotate("text", x = max(evaldf$p_dist_mat) * 0.7, y = Inf, vjust = 1.5,
             label = sprintf("Mean absolute errors = %.2f", 
                           mean(abs(evaldf$residuals))),
             color = "blue", size = 3) +
    annotate("text", x = max(evaldf$p_dist_mat) * 0.7, y = Inf, vjust = 3,
             label = sprintf("SD of errors = %.2f", 
                           sd(evaldf$residuals)),
             color = "blue", size = 3)
  
  # Save plots if requested
  if(save_plot && !is.na(scenario_name) && !is.na(ndim)) {
    if (is.null(output_dir)) {
      output_dir <- getwd()
    }
    
    # Save scatter plot
    scatter_file <- file.path(output_dir,
                              sprintf("%s_prediction_scatter_dim_%d.png", 
                                      scenario_name, ndim))
    ggsave(scatter_file, scatter_plot, width = 8, height = 8, dpi = 300)
    
    # Save residuals plot
    residuals_file <- file.path(output_dir,
                                sprintf("%s_residuals_vs_fitted_dim_%d.png",
                                        scenario_name, ndim))
    ggsave(residuals_file, residuals_plot, width = 8, height = 8, dpi = 300)
  }
  
  # Return both plots
  return(list(
    scatter_plot = scatter_plot,
    residuals_plot = residuals_plot
  ))
}