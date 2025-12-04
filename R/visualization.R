# Copyright (c) 2024 Omid Arhami omid.arhami@uga.edu
# R/visualization.R

#' Visualization functions for the topolow package


#' Plot Annotation Configuration Class
#'
#' @description
#' S3 class for configuring point annotations in plots, including labels,
#' connecting lines, and visual properties.
#'
#' @param notable_points Character vector of notable points to highlight
#' @param size Numeric. Size of annotations for notable points
#' @param color Character. Color of annotations for notable points
#' @param alpha Numeric. Alpha transparency of annotations
#' @param fontface Character. Font face of annotations ("plain", "bold", "italic", etc.)
#' @param box Logical. Whether to draw a box around annotations
#' @param segment_size Numeric. Size of segments connecting annotations to points
#' @param segment_alpha Numeric. Alpha transparency of connecting segments
#' @param min_segment_length Numeric. Minimum length of connecting segments
#' @param max_overlaps Numeric. Maximum number of overlaps allowed for annotations
#' @param outline_size Numeric. Size of the outline for annotations
#' @return An S3 object of class \code{annotation_config}, which is a list 
#'         containing the specified configuration parameters for plot annotations.
#' @export
new_annotation_config <- function(
    notable_points = NULL,
    size = 4.9,
    color = "black",
    alpha = 0.9,
    fontface = "plain",
    box = FALSE,
    segment_size = 0.3,
    segment_alpha = 0.6,
    min_segment_length = 0,
    max_overlaps = Inf,
    outline_size = 0.4
) {
    config <- list(
        notable_points = notable_points,
        size = size,
        color = color,
        alpha = alpha,
        fontface = fontface, 
        box = box,
        segment_size = segment_size,
        segment_alpha = segment_alpha,
        min_segment_length = min_segment_length,
        max_overlaps = max_overlaps,
        outline_size = outline_size
    )
    
    # Validate inputs
    stopifnot(
        is.null(notable_points) || is.character(notable_points),
        is.numeric(size), size > 0,
        is.character(color),
        is.numeric(alpha), alpha >= 0, alpha <= 1,
        is.character(fontface),
        is.logical(box),
        is.numeric(segment_size), segment_size >= 0,
        is.numeric(segment_alpha), segment_alpha >= 0, segment_alpha <= 1,
        is.numeric(min_segment_length), min_segment_length >= 0,
        is.numeric(max_overlaps), max_overlaps >= 0,
        is.numeric(outline_size), outline_size >= 0
    )
    
    structure(config, class = "annotation_config")
}


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
#' @param show_title Whether to show plot title (default: FALSE)
#' @param label_size Label text size 
#' @param title_size Title text size
#' @param subtitle_size Subtitle text size
#' @param axis_title_size Axis title text size
#' @param axis_text_size Axis text size
#' @param legend_text_size Legend text size
#' @param legend_title_size Legend title text size
#' @param show_legend Whether to show the legend
#' @param legend_position Legend position ("none", "right", "left", "top", "bottom")
#' @param arrow_head_size Size of the arrow head for velocity arrows (in cm)
#' @param arrow_alpha Transparency of arrows (0 = invisible, 1 = fully opaque)
#' @return An S3 object of class \code{aesthetic_config}, which is a list 
#'         containing the specified configuration parameters for plot aesthetics.
#' @export
new_aesthetic_config <- function(
    point_size = 3.5,
    point_alpha = 0.8,
    point_shapes = c(antigen = 16, antiserum = 0),
    color_palette = c25,
    gradient_colors = list(low = "blue", high = "red"),
    show_labels = FALSE,
    show_title = FALSE,
    label_size = 3,
    title_size = 14,
    subtitle_size = 12,
    axis_title_size = 12,
    axis_text_size = 10,
    legend_text_size = 10,
    legend_title_size = 12,
    show_legend = TRUE,
    legend_position = "right",
    arrow_head_size=0.2,
    arrow_alpha=0.6
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
    legend_position %in% c("none", "right", "left", "top", "bottom"),
    is.numeric(arrow_head_size), arrow_head_size >= 0, arrow_head_size <= 1,
    is.numeric(arrow_alpha), arrow_alpha >= 0, arrow_alpha <= 1
  )
    
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
    legend_position = legend_position,
    arrow_head_size= arrow_head_size,
    arrow_alpha= arrow_alpha
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
#' @param save_plot Logical. Whether to save the plot to a file.
#' @param save_format Plot save format ("png", "pdf", "svg", "eps")
#' @param reverse_x Numeric multiplier for x-axis direction (1 or -1)
#' @param reverse_y Numeric multiplier for y-axis direction (1 or -1)
#' @param x_limits Numeric vector of length 2 specifying c(min, max) for x-axis. If NULL, limits are set automatically.
#' @param y_limits Numeric vector of length 2 specifying c(min, max) for y-axis. If NULL, limits are set automatically.
#' @param arrow_plot_threshold Threshold for velocity arrows to be drawn in the same antigenic distance unit (default: 0.10)
#' @return An S3 object of class \code{layout_config}, which is a list containing 
#'         the specified configuration parameters for plot layout.
#' @importFrom ggplot2 margin
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
  save_plot = FALSE,
  save_format = "png",
  reverse_x = 1,
  reverse_y = 1,
  x_limits = NULL,
  y_limits = NULL,
  arrow_plot_threshold   = 0.1  # velocity arrows larger than this are shown
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
    save_plot = save_plot,
    save_format = save_format,
    reverse_x = reverse_x,
    reverse_y = reverse_y,
    x_limits = x_limits,
    y_limits = y_limits,
    arrow_plot_threshold   = arrow_plot_threshold
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
    ggplot2::is_margin(plot_margin),
    coord_type %in% c("fixed", "equal", "flip", "polar"),
    is.character(background_color),
    is.character(panel_background_color),
    is.logical(panel_border),
    is.character(panel_border_color),
    is.logical(save_plot),
    save_format %in% c("png", "pdf", "svg", "eps"),
    reverse_x %in% c(1, -1),
    reverse_y %in% c(1, -1),
    is.numeric(arrow_plot_threshold), arrow_plot_threshold >= 0
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
#' @return An S3 object of class \code{dim_reduction_config}, which is a list 
#'         containing the specified configuration parameters for dimensionality reduction.
#' @export
new_dim_reduction_config <- function(
    method = "pca",
    n_components = 2,
    scale = TRUE,
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
      mapping_max_iter = 1000,
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
#' @importFrom stats na.omit
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
#' @importFrom stats prcomp dist
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
                                                  scale. = FALSE,
                                                  center = TRUE),
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
                      if (!requireNamespace("umap", quietly = TRUE)) {
                        stop("umap package required for UMAP reduction. Please install it.")
                      }
                     umap_result <- do.call(umap::umap,
                                            c(list(d = coords,
                                                   n_components = config$n_components),
                                              config$umap_params))
                     data.frame(umap_result$layout)
                   },
                   "tsne" = {
                     if (!requireNamespace("Rtsne", quietly = TRUE)) {
                       stop("Rtsne package required for t-SNE reduction. Please install it.")
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
#' @importFrom stats dist optimize
#' @keywords internal
scale_to_original_distances <- function(reduced_coords, orig_dist) {
  # Calculate reduced distance matrix
  reduced_dist <- as.matrix(dist(reduced_coords))
  
  # Define optimization function for scaling factor
  sum_squared_diff <- function(scale_factor) {
    scaled_dist <- reduced_dist * scale_factor
    sum(abs(orig_dist - scaled_dist))
  }
  
  # Find optimal scaling factor
  optimal_scale <- optimize(sum_squared_diff, 
                            interval = c(0.01, 40))$minimum
  
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
#' @importFrom ggplot2 theme_minimal theme element_text element_blank element_line element_rect
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
#' Antigenic Mapping and Antigenic Velocity Function. Creates a visualization of points colored by time (year) using dimension reduction, with optional antigenic velocity arrows.
#' Points are colored on a gradient scale based on their temporal values, with
#' different shapes for antigens and antisera.
#'
#' @param df_coords Data frame containing:
#'        - V1, V2, ... Vn: Coordinate columns
#'        - antigen: Binary indicator for antigen points 
#'        - antiserum: Binary indicator for antiserum points
#'        - year: Numeric year values for temporal coloring
#' @param ndim Number of dimensions in input coordinates
#' @param draw_arrows logical; if TRUE, compute and draw antigenic drift vectors
#' @param annotate_arrows logical; if TRUE, show names of the points having arrows
#' @param phylo_tree Optional; phylo object in Newick format. Does not need to be rooted. If provided, used to compute antigenic velocity arrows.
#' @param clade_node_depth Optional; integer; number of levels of parent nodes to define clades. Antigens from different clades will be excluded from the calculation antigenic velocity arrows. (Default: Automatically calculated mode of leaf-to-backbone distance of the tree)
#' @param dim_config Dimension reduction configuration object specifying method and parameters
#' @param aesthetic_config Aesthetic configuration object controlling plot appearance
#' @param layout_config Layout configuration object controlling plot dimensions and style.
#'        Use x_limits and y_limits in layout_config to set axis limits.
#' @param output_dir Character. Directory for output files. Required if `layout_config$save_plot` is `TRUE`.
#' @param show_shape_legend Logical. Whether to show the shape legend (default: TRUE)
#' @param annotation_config Annotation configuration object for labeling notable points
#' @param sigma_t Optional; numeric; bandwidth for the Gaussian kernel discounting on time in years or the time unit of the data. If NULL, uses Silverman's rule of thumb.
#' @param sigma_x Optional; numeric; bandwidth for the Gaussian kernel discounting on antigenic distancein antigenic units. If NULL, uses Silverman's rule of thumb.
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
#' @return A \code{ggplot} object containing the temporal mapping visualization.
#'
#' @examples
#' # Basic usage with default configurations
#' data <- data.frame(
#'   V1 = rnorm(100), V2 = rnorm(100), V3 = rnorm(100), name = 1:100,
#'   antigen = rep(c(0,1), 50), antiserum = rep(c(1,0), 50),
#'   year = rep(2000:2009, each=10)
#' )
#' # Plot without saving
#' p1 <- plot_temporal_mapping(data, ndim=3)
#'
#' # Save plot to a temporary directory
#' temp_dir <- tempdir()
#' layout_config_save <- new_layout_config(save_plot = TRUE,
#'                        x_limits = c(-10, 10),
#'                        y_limits = c(-8, 8))
#' p_saved <- plot_temporal_mapping(data, ndim = 3, layout_config = layout_config_save, 
#'                                  output_dir = temp_dir)
#' list.files(temp_dir) # Check that file was created
#' unlink(temp_dir, recursive = TRUE) # Clean up
#'
#' @seealso 
#' \code{\link{plot_cluster_mapping}} for cluster-based visualization
#' \code{\link{plot_3d_mapping}} for 3D visualization
#' \code{\link{new_dim_reduction_config}} for dimension reduction options
#' \code{\link{new_aesthetic_config}} for aesthetic options
#' \code{\link{new_layout_config}} for layout options
#' \code{\link{new_annotation_config}} for annotation options
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_segment scale_colour_gradient scale_shape_manual scale_x_continuous scale_y_continuous labs
#' @importFrom grid unit arrow
#' @importFrom stats dist bw.nrd median
#' @export
plot_temporal_mapping <- function(df_coords, ndim, 
                                  dim_config = new_dim_reduction_config(),
                                  aesthetic_config = new_aesthetic_config(),
                                  layout_config = new_layout_config(),
                                  annotation_config = new_annotation_config(),
                                  output_dir,
                                  show_shape_legend = TRUE,
                                  draw_arrows = FALSE,
                                  annotate_arrows = TRUE,
                                  phylo_tree = NULL,
                                  sigma_t = NULL,
                                  sigma_x = NULL,
                                  clade_node_depth = NULL) {
  
  # Ensure ggrepel is available
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    warning("The ggrepel package is required for optimal label placement. Install with: install.packages('ggrepel')")
  }
  
  if (layout_config$save_plot && missing(output_dir)) {
    stop("An 'output_dir' must be provided when 'layout_config$save_plot' is TRUE.", call. = FALSE)
  }
  
  # Validate input data
  df_coords <- validate_topolow_df(df_coords, ndim, require_temporal = TRUE)
  
  # Perform dimension reduction
  reduced_df <- reduce_dimensions(df_coords, dim_config)
  
  # Apply axis reversals
  reduced_df$plot_x <- reduced_df$dim2 * layout_config$reverse_x
  reduced_df$plot_y <- reduced_df$dim1 * layout_config$reverse_y
  
  # Process row names to get strain names
  reduced_df$clean_name <- sub("^(V/|S/)", "", reduced_df$name)
  
  # Create base theme
  base_theme <- create_base_theme(aesthetic_config, layout_config)
  
  # Flag notable points
  if (!is.null(annotation_config$notable_points) && length(annotation_config$notable_points) > 0) {
    reduced_df$is_notable <- reduced_df$clean_name %in% annotation_config$notable_points
    annotation_df <- reduced_df[reduced_df$is_notable, ]
    
    if (nrow(annotation_df) == 0) {
      warning("None of the specified notable points found in the data")
    }
  } else {
    reduced_df$is_notable <- FALSE
    annotation_df <- reduced_df[0, ]  # Empty dataframe
  }
  
  # Create point type with explicit factor levels - this is just for regular points
  reduced_df$point_type <- NA_character_  # Initialize
  reduced_df$point_type[reduced_df$antigen & !reduced_df$is_notable] <- "antigen"
  reduced_df$point_type[reduced_df$antiserum & !reduced_df$is_notable] <- "antiserum"
  reduced_df$point_type <- factor(reduced_df$point_type, 
                                  levels = names(aesthetic_config$point_shapes))
  
  # Create plot
  p <- ggplot(
    data = reduced_df[!reduced_df$is_notable, ],
    aes(x = .data$plot_x, y = .data$plot_y, color = .data$year,
        shape = .data$point_type)) +
    geom_point(size = aesthetic_config$point_size,
               alpha = aesthetic_config$point_alpha) +
    # Plot notable points with filled shapes and outlines
    geom_point(
      data = reduced_df[reduced_df$is_notable & reduced_df$antigen, ],
      aes(x = .data$plot_x, y = .data$plot_y, 
          fill = .data$year),
      shape = 21,  # Filled circle with outline
      color = "black",  # Outline color
      size = aesthetic_config$point_size * 1.2,
      stroke = annotation_config$outline_size,
      alpha = aesthetic_config$point_alpha
    ) +
    # Notable antisera with different filled shape
    geom_point(
      data = reduced_df[reduced_df$is_notable & reduced_df$antiserum, ],
      aes(x = .data$plot_x, y = .data$plot_y, 
          fill = .data$year),
      shape = 22,  # Filled square with outline
      color = "black",  # Outline color
      size = aesthetic_config$point_size * 1.2,
      stroke = annotation_config$outline_size,
      alpha = aesthetic_config$point_alpha
    ) +
    scale_colour_gradient(low = aesthetic_config$gradient_colors$low,
                          high = aesthetic_config$gradient_colors$high,
                          na.value = "gray50") +
    scale_shape_manual(
      name = "Type",
      values = aesthetic_config$point_shapes,
      labels = c(antigen = "Antigen", antiserum = "Antiserum")) +
    base_theme +
    labs(title = if(aesthetic_config$show_title) "Temporal Mapping" else "",
         x = "Dimension 1",
         y = "Dimension 2",
         colour = "Year")
  
  # Add labels to notable points if any exist
  if (nrow(annotation_df) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      if (annotation_config$box) {
        # Use label boxes with background
        p <- p + ggrepel::geom_label_repel(
          data = annotation_df,
          aes(x = .data$plot_x, y = .data$plot_y, label = .data$clean_name),
          size = annotation_config$size / ggplot2::.pt,
          color = annotation_config$color,
          alpha = annotation_config$alpha,
          fontface = annotation_config$fontface,
          segment.size = annotation_config$segment_size,
          segment.alpha = annotation_config$segment_alpha,
          min.segment.length = annotation_config$min_segment_length,
          max.overlaps = annotation_config$max_overlaps,
          box.padding = unit(0.4, "lines"),
          point.padding = unit(0.3, "lines"),
          force = 1,
          direction = "both"
        )
      } else {
        # Use simple text without background
        p <- p + ggrepel::geom_text_repel(
          data = annotation_df,
          aes(x = .data$plot_x, y = .data$plot_y, label = .data$clean_name),
          size = annotation_config$size / ggplot2::.pt,
          color = annotation_config$color,
          alpha = annotation_config$alpha,
          fontface = annotation_config$fontface,
          segment.size = annotation_config$segment_size,
          segment.alpha = annotation_config$segment_alpha,
          min.segment.length = annotation_config$min_segment_length,
          max.overlaps = annotation_config$max_overlaps,
          box.padding = unit(0.4, "lines"),
          point.padding = unit(0.3, "lines"),
          force = 1,
          direction = "both"
        )
      }
    } else {
      # Fallback if ggrepel is not available - basic text labels
      warning("ggrepel package not available - using basic text labels without repulsion")
      p <- p + geom_text(
        data = annotation_df,
        aes(x = .data$plot_x, y = .data$plot_y, label = .data$clean_name),
        size = annotation_config$size / ggplot2::.pt,
        color = annotation_config$color,
        alpha = annotation_config$alpha,
        fontface = annotation_config$fontface,
        nudge_x = 0.1,
        nudge_y = 0.1,
        check_overlap = TRUE
      )
    }
  }
  
  # Initialize velocity data as NULL
  velocity_data <- NULL
  if (draw_arrows) {
    if (!"year" %in% names(reduced_df)) {
      stop("`year` column is required when draw_arrows = TRUE")
    }
    # limit to antigens only:
    positions <- reduced_df[reduced_df$antigen, ]
    positions$V1 <- positions$plot_x
    positions$V2 <- positions$plot_y
    
    if (!is.null(phylo_tree)) {
      if (!requireNamespace("ape", quietly = TRUE)) {
        stop("Package 'ape' is required for this function. Please install it.")
      }
      if (ape::is.rooted(phylo_tree)) {
        phylo_tree <- ape::unroot(phylo_tree)
      }
      #-- identify which tip names actually exist in the tree
      positions$name <- toupper(positions$name)
      all_points      <- positions$name
      tree_tips_up  <- toupper(phylo_tree$tip.label)
      tip_idx <- match(toupper(positions$name), tree_tips_up)
      # tip_idx[i] is the row in D_edge (DN) for positions$name[i], or NA if absent
      
      # compare in uppercase for consistency
      present_mask <- toupper(all_points) %in% tree_tips_up
      tree_present_points <- unique(all_points[present_mask])
      absent_tips <- setdiff(all_points, tree_present_points)
      
      if (length(absent_tips) > 0) {
        cat(
          "\nThe following antigens are not in the provided phylo_tree\n ",
          "Thus, they did not contribute to limiting antigenic velo.\n",
          paste(absent_tips, collapse = ", "), "\n"
          
        )
      }
      
      # phylogenetic clade depth cutoff via leaf-to-backbone (longest-path in terms of nodes, aka tree spine)
      # 1) force every edge to length 1
      tree_unit <- phylo_tree
      tree_unit$edge.length <- rep(1, nrow(tree_unit$edge))
      
      # 2) compute node-to-node distances
      DN <- ape::dist.nodes(tree_unit)
      
      # 3) find the two tips with maximum separation
      n_tip <- length(tree_unit$tip.label)
      tip_idx <- seq_len(n_tip)
      tip_dist_mat <- DN[tip_idx, tip_idx, drop = FALSE]
      # pick the first pair achieving the max
      max_pair <- which(tip_dist_mat == max(tip_dist_mat), arr.ind = TRUE)[1, ]
      
      # 4) get the (internal + tip) nodes along that path
      path_nodes <- ape::nodepath(tree_unit, max_pair[1], max_pair[2])
      
      # 5) for each tip, its distance to the nearest node on that path
      leaf_distances <- apply(DN[tip_idx, path_nodes, drop = FALSE], 1, min)
      
      # 6) clade_node_depth = median distance (unless overridden)
      clade_node_depth <- ifelse(
        is.null(clade_node_depth),
        ceiling(median(leaf_distances)),
        clade_node_depth
      )
      
      cat(sprintf("Average leaf-to-backbone distance of the tree (%.3f) was used as clades' depth cutoff.\n", clade_node_depth))
      
      # 7) get the clade nodes for each tip            
      get_clade_node <- function(phy, tip_label, depth) {
        # exact match in uppercase
        ix <- match(toupper(tip_label), tree_tips_up)
        if (is.na(ix)) {
          stop("Internal error: tip '", tip_label, "' lookup failed")
        }
        node <- ix
        for (k in seq_len(depth)) {
          parent <- phy$edge[phy$edge[,2] == node, 1]
          if (length(parent) != 1) break
          node <- parent
        }
        node
      }
      
      # only compute clades for points that are present in the tree
      clade_nodes  <- setNames(
        lapply(tree_present_points, get_clade_node, phy = phylo_tree, depth = clade_node_depth),
        tree_present_points
      )
      clade_members <- lapply(
        clade_nodes,
        function(nd) ape::extract.clade(phylo_tree, nd)$tip.label
      )
    }
    
    # Estimate kernel bandwidths based on Silverman's rule (bw.nrd)
    # Use 'dists' is of class "dist" and length n*(n-1)/2
    distmat_t <- dist(positions$year)
    sigma_t <- ifelse(is.null(sigma_t), bw.nrd(distmat_t), sigma_t)
    
    distmat_x <- dist(positions[, c("V1", "V2")], method = "euclidean")
    sigma_x <- ifelse(is.null(sigma_x),
                      sqrt(0.5*( (bw.nrd(positions$V1))^2 + (bw.nrd(positions$V2))^2 )),
                      sigma_x)
    
    cat(sprintf(
      "Kernel bandwidth for time = %.3f\n", sigma_t
    ))
    cat(sprintf(
      "Kernel bandwidth for antigenic distance = %.3f\n", sigma_x
    ))
    
    
    
    n   <- nrow(positions)
    v1  <- numeric(n)
    v2  <- numeric(n)
    
    for (i in seq_len(n)) {
      this_pt <- positions$name[i]
      # determine past indices, excluding only known "non-clade" tips
      if (!is.null(phylo_tree) && this_pt %in% names(clade_members)) {
        # those present but *not* in the same clade
        bad <- setdiff(tree_present_points, clade_members[[this_pt]])
        past_idx <- which(
          positions$year < positions$year[i] &
            !(positions$name %in% bad)
        )
      } else {
        # either no tree or tip absent from tree -> include *all* past points
        past_idx <- which(positions$year < positions$year[i])
      }
      
      if (length(past_idx)) {
        #  positive dt = current - past
        dt <- positions$year[i] - positions$year[past_idx]
        dx <- positions$V1[i] - positions$V1[past_idx]
        dy <- positions$V2[i] - positions$V2[past_idx]
        w  <- exp(-(dx^2 + dy^2)/(2*sigma_x^2)) * exp(- (dt^2)/(2*sigma_t^2))
        v1[i] <- sum(w * (dx / dt)) / sum(w)
        v2[i] <- sum(w * (dy / dt)) / sum(w)
      } else {
        v1[i] <- NA
        v2[i] <- NA
      }
    }
    positions$v1  <- v1
    positions$v2  <- v2
    positions$mag <- sqrt(v1^2 + v2^2)
    # Store velocity data for return
    velocity_data <- positions[, c("name", "V1", "V2", "year", "v1", "v2", "mag")]
    
    top_vel <- subset(positions, mag >= layout_config$arrow_plot_threshold)
    cat(sprintf("Showing only arrows with magnitude >= %.3f (figure unit)\n", layout_config$arrow_plot_threshold))
    

    # -- overlay top-velocity points with filled shape + black outline --
    p <- p +
      geom_point(
        data        = top_vel,
        inherit.aes = FALSE,
        aes(x = .data$V1, y = .data$V2),
        shape       = 21,   # same filled-circle shape you used for notable antigens
        color       = "black",
        size        = aesthetic_config$point_size * 1.2,
        stroke      = annotation_config$outline_size,
        alpha       = aesthetic_config$point_alpha
      )
    
    # add arrow layer
    p <- p +
      geom_segment(
        data      = top_vel,
        inherit.aes = FALSE,
        aes(x    = .data$V1 - .data$v1,
            y    = .data$V2 - .data$v2,
            xend = .data$V1,
            yend = .data$V2),
        arrow = arrow(length = unit(aesthetic_config$arrow_head_size, "cm")),
        alpha = aesthetic_config$arrow_alpha
      )
    
    
    if (annotate_arrows) {
      # Add arrow labels (length in parentheses) at midpoint
      p <- p + 
        geom_text(
          data = top_vel,
          aes(
            x = .data$V1 - .data$v1 / 2,
            y = .data$V2 - .data$v2 / 2,
            label = sprintf("(%.2f)", .data$mag)
          ),
          size = 0.8*(annotation_config$size / ggplot2::.pt),  # Small font size
          hjust = 0.5,
          vjust = 0.5,
          alpha = 0.8*annotation_config$alpha
        )
      # Annotate top-velocity points exactly like notable-point labels
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        if (annotation_config$box) {
          p <- p +
            ggrepel::geom_label_repel(
              data        = top_vel,
              inherit.aes = FALSE,
              aes(x = .data$V1, y = .data$V2, label = .data$name),
              size              = annotation_config$size / ggplot2::.pt,
              color             = annotation_config$color,
              alpha             = annotation_config$alpha,
              fontface          = annotation_config$fontface,
              segment.size      = annotation_config$segment_size,
              segment.alpha     = annotation_config$segment_alpha,
              min.segment.length= annotation_config$min_segment_length,
              max.overlaps      = annotation_config$max_overlaps,
              box.padding       = unit(0.4, "lines"),
              point.padding     = unit(0.3, "lines"),
              force             = 1,
              direction         = "both"
            )
        } else {
          p <- p +
            ggrepel::geom_text_repel(
              data        = top_vel,
              inherit.aes = FALSE,
              aes(x = .data$V1, y = .data$V2, label = .data$name),
              size              = annotation_config$size / ggplot2::.pt,
              color             = annotation_config$color,
              alpha             = annotation_config$alpha,
              fontface          = annotation_config$fontface,
              segment.size      = annotation_config$segment_size,
              segment.alpha     = annotation_config$segment_alpha,
              min.segment.length= annotation_config$min_segment_length,
              max.overlaps      = annotation_config$max_overlaps,
              box.padding       = unit(0.4, "lines"),
              point.padding     = unit(0.3, "lines"),
              force             = 1,
              direction         = "both"
            )
        }
      } else {
        warning("ggrepel package not available - using basic text labels without repulsion")
        p <- p + geom_text(
          data = top_vel,
          aes(x = .data$V1, y = .data$V2, label = .data$name),
          size           = annotation_config$size / ggplot2::.pt,
          color          = annotation_config$color,
          alpha          = annotation_config$alpha,
          fontface       = annotation_config$fontface,
          nudge_x        = - 0.15,
          nudge_y        = - 0.15,
          check_overlap  = TRUE
        )
      }
    }
  }
  
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
  
  # Save plot if requested
  if (layout_config$save_plot) {
    if (draw_arrows){
      if (!is.null(phylo_tree)) {
        filename <- sprintf(
          "temporal_map_ndim_%d_s_t_%g_s_x_%g_cladeDepth_%g_arrowthresh_%g.%s",
          ndim,
          sigma_t,
          sigma_x,
          clade_node_depth,
          layout_config$arrow_plot_threshold,
          layout_config$save_format
        )
      } else {
        filename <- sprintf(
          "temporal_map_ndim_%d_s_t_%g_s_x_%g_arrowthresh_%g.%s",
          ndim,
          sigma_t,
          sigma_x,
          layout_config$arrow_plot_threshold,
          layout_config$save_format
        )
      }
    } else {
      filename <- sprintf(
        "temporal_map_ndim_%d.%s",
        ndim,
        layout_config$save_format
      )
    }
    save_plot(p, filename, layout_config, output_dir)
  }
  
  # Return both plot and velocity data
  if (!is.null(velocity_data)) {
    return(list(
      plot = p,
      velocity_data = velocity_data
    ))
  } else {
    return(p)
  }
}


#' Create Clustered Mapping Plots
#'
#' @description
#' Antigenic Mapping and Antigenic Velocity Function. Creates a visualization of points colored by cluster assignment using dimension 
#' reduction, with optional antigenic velocity arrows. Points are colored by cluster with different shapes for antigens and 
#' antisera.
#'
#' @param df_coords Data frame containing:
#'        - V1, V2, ... Vn: Coordinate columns
#'        - antigen: Binary indicator for antigen points 
#'        - antiserum: Binary indicator for antiserum points
#'        - cluster: Factor or integer cluster assignments
#' @param ndim Number of dimensions in input coordinates
#' @param draw_arrows logical; if TRUE, compute and draw antigenic drift vectors
#' @param annotate_arrows logical; if TRUE, show names of the points having arrows
#' @param phylo_tree Optional; phylo object in Newick format. Does not need to be rooted. If provided, used to compute antigenic velocity arrows.
#' @param clade_node_depth Optional; integer; number of levels of parent nodes to define clades. Antigens from different clades will be excluded from the calculation antigenic velocity arrows. (Default: Automatically calculated mode of leaf-to-backbone distance of the tree)
#' @param dim_config Dimension reduction configuration object specifying method and parameters
#' @param aesthetic_config Aesthetic configuration object controlling plot appearance
#' @param layout_config Layout configuration object controlling plot dimensions and style.
#'        Use x_limits and y_limits in layout_config to set axis limits.
#' @param output_dir Character. Directory for output files. Required if `layout_config$save_plot` is `TRUE`.
#' @param show_shape_legend Logical. Whether to show the shape legend (default: TRUE)
#' @param cluster_legend_title Character. Custom title for the cluster legend (default: "Cluster")
#' @param annotation_config Annotation configuration object for labeling notable points
#' @param sigma_t Optional; numeric; bandwidth for the Gaussian kernel discounting on time in years or the time unit of the data. If NULL, uses Silverman's rule of thumb.
#' @param sigma_x Optional; numeric; bandwidth for the Gaussian kernel discounting on antigenic distance in antigenic units. If NULL, uses Silverman's rule of thumb.
#' @param show_one_arrow_per_cluster Shows only the largest antigenic velocity arrow in each cluster
#' @param cluster_legend_order in case you prefer a certain order for clusters in the legend, 
#'        provide a list with that order here; e.g., c("cluster 2", "cluster 1")
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
#' @return A \code{ggplot} object containing the cluster mapping visualization.
#'
#' @examples
#' # Basic usage with default configurations
#' data <- data.frame(
#'   V1 = rnorm(100), V2 = rnorm(100), V3 = rnorm(100), name = 1:100,
#'   antigen = rep(c(0,1), 50), antiserum = rep(c(1,0), 50),
#'   cluster = rep(1:5, each=20)
#' )
#' p1 <- plot_cluster_mapping(data, ndim=3)
#'
#' # Save plot to a temporary directory
#' temp_dir <- tempdir()
#' # Custom configurations with specific color palette and axis limits
#' aesthetic_config <- new_aesthetic_config(
#'   point_size = 4,
#'   point_alpha = 0.7,
#'   color_palette = c("red", "blue", "green", "purple", "orange"),
#'   show_labels = TRUE,
#'   label_size = 3
#' )
#'
#' layout_config_save <- new_layout_config(save_plot = TRUE,
#'   width = 10,
#'   height = 8,
#'   coord_type = "fixed",
#'   show_grid = TRUE,
#'   grid_type = "major",
#'   x_limits = c(-10, 10),
#'   y_limits = c(-8, 8)
#' )
#'
#' p_saved <- plot_cluster_mapping(data, ndim=3, 
#'   layout_config = layout_config_save, 
#'   aesthetic_config = aesthetic_config,
#'   output_dir = temp_dir
#' )
#' 
#' list.files(temp_dir)
#' unlink(temp_dir, recursive = TRUE)
#'
#' @seealso 
#' \code{\link{plot_temporal_mapping}} for temporal visualization
#' \code{\link{plot_3d_mapping}} for 3D visualization
#' \code{\link{new_dim_reduction_config}} for dimension reduction options
#' \code{\link{new_aesthetic_config}} for aesthetic options
#' \code{\link{new_layout_config}} for layout options
#' \code{\link{new_annotation_config}} for annotation options
#' @importFrom ggplot2 ggplot aes coord_fixed theme geom_point geom_text geom_segment scale_colour_manual scale_fill_manual scale_shape_manual guides guide_legend labs scale_x_continuous scale_y_continuous
#' @importFrom grid unit arrow
#' @importFrom dplyr group_by filter ungroup
#' @importFrom stats dist bw.nrd median
#' @export
plot_cluster_mapping <- function(df_coords, ndim,
                                 dim_config = new_dim_reduction_config(),
                                 aesthetic_config = new_aesthetic_config(),
                                 layout_config = new_layout_config(),
                                 annotation_config = new_annotation_config(),
                                 output_dir,
                                 show_shape_legend = TRUE,
                                 cluster_legend_title = "Cluster",
                                 draw_arrows = FALSE,
                                 annotate_arrows = TRUE,
                                 phylo_tree = NULL,
                                 sigma_t = NULL,
                                 sigma_x = NULL,
                                 clade_node_depth = NULL,
                                 show_one_arrow_per_cluster = FALSE,
                                 cluster_legend_order = NULL) {
  
  # Ensure ggrepel is available
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    warning("The ggrepel package is required for optimal label placement. Install with: install.packages('ggrepel')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The dplyr package is required for this function. Please install it.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting. Please install it.")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("The grid package is required for arrow specifications. Please install it.")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("The stats package is required for statistical operations. Please install it.")
  }

  if (layout_config$save_plot && missing(output_dir)) {
    stop("An 'output_dir' must be provided when 'layout_config$save_plot' is TRUE.", call. = FALSE)
  }
  
  # Validate input data
  df_coords <- validate_topolow_df(df_coords, ndim, require_clusters = TRUE)
  
  # Perform dimension reduction
  reduced_df <- reduce_dimensions(df_coords, dim_config)
  
  # Apply axis reversals
  reduced_df$plot_x <- reduced_df$dim2 * layout_config$reverse_x
  reduced_df$plot_y <- reduced_df$dim1 * layout_config$reverse_y
  
  # Process row names to get strain names
  reduced_df$clean_name <- sub("^(V/|S/)", "", reduced_df$name)
  
  #NOTABLE-POINT & ARROW FLAGS
  reduced_df$is_notable <- FALSE
  reduced_df$is_arrow   <- FALSE
  
  if (!is.null(annotation_config$notable_points) &&
      length(annotation_config$notable_points) > 0) {
    reduced_df$is_notable <- reduced_df$clean_name %in% annotation_config$notable_points
  }
  
  ## cluster legend housekeeping 
  # Apply custom cluster legend order if provided
  if (!is.null(cluster_legend_order)) {
    reduced_df$cluster <- factor(reduced_df$cluster, levels = cluster_legend_order)
  }
  
  if (!is.null(cluster_legend_order)) {
    cluster_levels <- cluster_legend_order
  } else {
    cluster_levels <- sort(unique(reduced_df$cluster))
  }
  
  n_clusters <- length(cluster_levels)
  colors <- aesthetic_config$color_palette[1:min(n_clusters, length(aesthetic_config$color_palette))]
  names(colors) <- cluster_levels
  
  if (n_clusters > length(aesthetic_config$color_palette)) {
    warning("More clusters than available colors. Colors will be recycled.")
    colors <- rep_len(aesthetic_config$color_palette, n_clusters)
  }
  
  reduced_df$cluster <- factor(reduced_df$cluster, levels = cluster_levels)
  
  ## --------------------------
  ##  BASE PLOT (without arrows / notable overlays yet)
  ## --------------------------
  # Create base theme
  base_theme <- create_base_theme(aesthetic_config, layout_config)
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
  
  # Flag notable points
  if (!is.null(annotation_config$notable_points) && length(annotation_config$notable_points) > 0) {
    reduced_df$is_notable <- reduced_df$clean_name %in% annotation_config$notable_points
    annotation_df <- reduced_df[reduced_df$is_notable, ]
    
    if (nrow(annotation_df) == 0) {
      warning("None of the specified notable points found in the data")
    }
  } else {
    reduced_df$is_notable <- FALSE
    annotation_df <- reduced_df[0, ]  # Empty dataframe
  }
  
  # Create point type with explicit factor levels - this is just for regular points
  reduced_df$point_type <- NA_character_  # Initialize
  reduced_df$point_type[reduced_df$antigen & !reduced_df$is_notable] <- "antigen"
  reduced_df$point_type[reduced_df$antiserum & !reduced_df$is_notable] <- "antiserum"
  reduced_df$point_type <- factor(reduced_df$point_type, 
                                  levels = names(aesthetic_config$point_shapes))
  
  # Create plot
  p <- ggplot()
  
  ## Regular points   ------ EXCLUDING notable and arrow points ------ 
  p <- p +
    geom_point(
      data = reduced_df[!reduced_df$is_notable & !reduced_df$is_arrow, ],
      aes(x = .data$plot_x, y = .data$plot_y, 
          color = .data$cluster,
          shape = .data$point_type),
      size = aesthetic_config$point_size,
      alpha = aesthetic_config$point_alpha
    ) +
    
    # Plot notable points with filled shapes and outlines
    geom_point(
      data = reduced_df[reduced_df$is_notable & reduced_df$antigen, ],
      aes(x = .data$plot_x, y = .data$plot_y, 
          fill = .data$cluster),
      shape = 21,  # Filled circle with outline
      color = "black",  # Outline color
      size = aesthetic_config$point_size * 1.2,
      stroke = annotation_config$outline_size,
      alpha = aesthetic_config$point_alpha
    ) +
    # Notable antisera with different filled shape
    geom_point(
      data = reduced_df[reduced_df$is_notable & reduced_df$antiserum, ],
      aes(x = .data$plot_x, y = .data$plot_y, 
          fill = .data$cluster),
      shape = 22,  # Filled square with outline
      color = "black",  # Outline color
      size = aesthetic_config$point_size * 1.2,
      stroke = annotation_config$outline_size,
      alpha = aesthetic_config$point_alpha
    ) +
    # Configure scales
    scale_colour_manual(name = cluster_legend_title, values = colors) +
    scale_fill_manual(name = cluster_legend_title, values = colors, guide = "none") +
    scale_shape_manual(
      name = "Type",
      values = aesthetic_config$point_shapes,
      labels = c(antigen = "Antigen", antiserum = "Antiserum"),
      guide = if(show_shape_legend) "legend" else "none"
    ) +
    guides(colour = if(aesthetic_config$show_legend) 
      guide_legend(ncol = n_legend_cols, title.position = "top", byrow = TRUE)
      else "none") +
    base_theme +
    legend_theme +
    labs(title = if(aesthetic_config$show_title) "Cluster Mapping" else "",
         x = "Dimension 1", 
         y = "Dimension 2",
         colour = "Cluster")
  
  # Add labels to notable points if any exist
  if (nrow(annotation_df) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      if (annotation_config$box) {
        # Use label boxes with background
        p <- p + ggrepel::geom_label_repel(
          data = annotation_df,
          aes(x = .data$plot_x, y = .data$plot_y, label = .data$clean_name),
          size = annotation_config$size / ggplot2::.pt,
          color = annotation_config$color,
          alpha = annotation_config$alpha,
          fontface = annotation_config$fontface,
          segment.size = annotation_config$segment_size,
          segment.alpha = annotation_config$segment_alpha,
          min.segment.length = annotation_config$min_segment_length,
          max.overlaps = annotation_config$max_overlaps,
          box.padding = unit(0.4, "lines"),
          point.padding = unit(0.3, "lines"),
          force = 1,
          direction = "both"
        )
      } else {
        # Use simple text without background
        p <- p + ggrepel::geom_text_repel(
          data = annotation_df,
          aes(x = .data$plot_x, y = .data$plot_y, label = .data$clean_name),
          size = annotation_config$size / ggplot2::.pt,
          color = annotation_config$color,
          alpha = annotation_config$alpha,
          fontface = annotation_config$fontface,
          segment.size = annotation_config$segment_size,
          segment.alpha = annotation_config$segment_alpha,
          min.segment.length = annotation_config$min_segment_length,
          max.overlaps = annotation_config$max_overlaps,
          box.padding = unit(0.4, "lines"),
          point.padding = unit(0.3, "lines"),
          force = 1,
          direction = "both"
        )
      }
    } else {
      # Fallback if ggrepel is not available - basic text labels
      warning("ggrepel package not available - using basic text labels without repulsion")
      p <- p + geom_text(
        data = annotation_df,
        aes(x = .data$plot_x, y = .data$plot_y, label = .data$clean_name),
        size = annotation_config$size / ggplot2::.pt,
        color = annotation_config$color,
        alpha = annotation_config$alpha,
        fontface = annotation_config$fontface,
        nudge_x = 0.1,
        nudge_y = 0.1,
        check_overlap = TRUE
      )
    }
  }
  
  # Initialize velocity data as NULL
  velocity_data <- NULL
  if (draw_arrows) {
    if (!"year" %in% names(reduced_df)) {
      stop("`year` column is required when draw_arrows = TRUE")
    }
    # limit to antigens only:
    positions <- reduced_df[reduced_df$antigen, ]
    positions$V1 <- positions$plot_x
    positions$V2 <- positions$plot_y
    
    if (!is.null(phylo_tree)) {
      if (!requireNamespace("ape", quietly = TRUE)) {
        stop("Package 'ape' is required for this function. Please install it.")
      }
      if (ape::is.rooted(phylo_tree)) {
        phylo_tree <- ape::unroot(phylo_tree)
      }
      #-- identify which tip names actually exist in the tree
      positions$name <- toupper(positions$name)
      all_points      <- positions$name
      tree_tips_up  <- toupper(phylo_tree$tip.label)
      tip_idx <- match(toupper(positions$name), tree_tips_up)
      # tip_idx[i] is the row in D_edge (DN) for positions$name[i], or NA if absent
      
      # compare in uppercase for consistency
      present_mask <- toupper(all_points) %in% tree_tips_up
      tree_present_points <- unique(all_points[present_mask])
      absent_tips <- setdiff(all_points, tree_present_points)
      
      if (length(absent_tips) > 0) {
        cat(
          "\nThe following antigens are not in the provided phylo_tree\n ",
          "Thus, they did not contribute to limiting antigenic velo.\n",
          paste(absent_tips, collapse = ", "), "\n"
          
        )
      }
      
      # phylogenetic clade depth cutoff via leaf-to-backbone (longest-path in terms of edges, aka tree spine)
      # 1) force every edge to length 1
      tree_unit <- phylo_tree
      tree_unit$edge.length <- rep(1, nrow(tree_unit$edge))
      
      # 2) compute node-to-node distances
      DN <- ape::dist.nodes(tree_unit)
      
      # 3) find the two tips with maximum separation
      n_tip <- length(tree_unit$tip.label)
      tip_idx <- seq_len(n_tip)
      tip_dist_mat <- DN[tip_idx, tip_idx, drop = FALSE]
      # pick the first pair achieving the max
      max_pair <- which(tip_dist_mat == max(tip_dist_mat), arr.ind = TRUE)[1, ]
      
      # 4) get the (internal + tip) nodes along that path
      path_nodes <- ape::nodepath(tree_unit, max_pair[1], max_pair[2])
      
      # 5) for each tip, its distance to the nearest node on that path
      leaf_distances <- apply(DN[tip_idx, path_nodes, drop = FALSE], 1, min)
      
      # 6) clade_node_depth = median distance (unless overridden)
      clade_node_depth <- ifelse(
        is.null(clade_node_depth),
        ceiling(median(leaf_distances)),
        clade_node_depth
      )
      
      cat(sprintf("Average leaf-to-backbone distance of the tree (%.3f) was used as clades' depth cutoff.\n", clade_node_depth))
      
      # 7) get the clade nodes for each tip      
      get_clade_node <- function(phy, tip_label, depth) {
        # exact match in uppercase
        ix <- match(toupper(tip_label), tree_tips_up)
        if (is.na(ix)) {
          stop("Internal error: tip '", tip_label, "' lookup failed")
        }
        node <- ix
        for (k in seq_len(depth)) {
          parent <- phy$edge[phy$edge[,2] == node, 1]
          if (length(parent) != 1) break
          node <- parent
        }
        node
      }
      
      # only compute clades for points that are present in the tree
      clade_nodes  <- setNames(
        lapply(tree_present_points, get_clade_node, phy = phylo_tree, depth = clade_node_depth),
        tree_present_points
      )
      clade_members <- lapply(
        clade_nodes,
        function(nd) ape::extract.clade(phylo_tree, nd)$tip.label
      )
    }
    
    # Estimate kernel bandwidths based on Silverman's rule (bw.nrd)
    # Use 'dists' is of class "dist" and length n*(n-1)/2
    distmat_t <- dist(positions$year)
    sigma_t <- ifelse(is.null(sigma_t), bw.nrd(distmat_t), sigma_t)
    
    distmat_x <- dist(positions[, c("V1", "V2")], method = "euclidean")
    sigma_x <- ifelse(is.null(sigma_x),
                      sqrt(0.5*( (bw.nrd(positions$V1))^2 + (bw.nrd(positions$V2))^2 )),
                      sigma_x)
    
    cat(sprintf(
      "Kernel bandwidth for time = %.3f\n", sigma_t
    ))
    cat(sprintf(
      "Kernel bandwidth for antigenic distance = %.3f\n", sigma_x
    ))
    
    
    
    n   <- nrow(positions)
    v1  <- numeric(n)
    v2  <- numeric(n)
    
    for (i in seq_len(n)) {
      this_pt <- positions$name[i]
      # determine past indices, excluding only known "non-clade" tips
      if (!is.null(phylo_tree) && this_pt %in% names(clade_members)) {
        # those present but *not* in the same clade
        bad <- setdiff(tree_present_points, clade_members[[this_pt]])
        past_idx <- which(
          positions$year < positions$year[i] &
            !(positions$name %in% bad)
        )
      } else {
        # either no tree or tip absent from tree -> include *all* past points
        past_idx <- which(positions$year < positions$year[i])
      }
      
      if (length(past_idx)) {
        #  positive dt = current - past
        dt <- positions$year[i] - positions$year[past_idx]
        dx <- positions$V1[i] - positions$V1[past_idx]
        dy <- positions$V2[i] - positions$V2[past_idx]
        w  <- exp(-(dx^2 + dy^2)/(2*sigma_x^2)) * exp(- (dt^2)/(2*sigma_t^2)) 
        v1[i] <- sum(w * (dx / dt)) / sum(w)
        v2[i] <- sum(w * (dy / dt)) / sum(w)
      } else {
        v1[i] <- NA
        v2[i] <- NA
      }
    }
    positions$v1  <- v1
    positions$v2  <- v2
    positions$mag <- sqrt(v1^2 + v2^2)
    
    # Store velocity data for return
    velocity_data <- positions[, c("name", "V1", "V2", "year", "v1", "v2", "mag")]
    
    if (show_one_arrow_per_cluster) {
      top_vel <- positions %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(mag == max(mag, na.rm = TRUE)) %>%
        dplyr::ungroup()
      cat("Showing one longest arrow per cluster\n")
    } else {
      top_vel <- subset(positions, mag >= layout_config$arrow_plot_threshold)
      cat(sprintf("Showing only arrows with magnitude >= %.3f (figure unit)\n", layout_config$arrow_plot_threshold))
    }
    
    top_vel$cluster <- factor(top_vel$cluster, levels = cluster_levels)
    
    # mark `top_vel` rows in reduced_df as arrows
    reduced_df$is_arrow[match(top_vel$name, reduced_df$name)] <- TRUE
    
    # -- overlay top-velocity points with filled shape + black outline --
    p <- p +
      geom_point(
        data        = top_vel,
        inherit.aes = FALSE,
        aes(x = .data$V1, y = .data$V2, fill = .data$cluster),
        shape       = 21,   # same filled-circle shape used for notable antigens
        color       = "black",
        size        = aesthetic_config$point_size * 1.2,
        stroke      = annotation_config$outline_size,
        alpha       = aesthetic_config$point_alpha
      )
    
    # add arrow layer
    p <- p +
      geom_segment(
        data      = top_vel,
        inherit.aes = FALSE,
        aes(x = .data$V1 - .data$v1,
            y = .data$V2 - .data$v2,
            xend = .data$V1,
            yend = .data$V2),
        arrow = arrow(length = unit(aesthetic_config$arrow_head_size, "cm")),
        alpha = aesthetic_config$arrow_alpha
      )
    
    if (annotate_arrows) {
      # Add arrow labels (length in parentheses) at midpoint
      # p <- p + 
      # geom_text(
      #   data = top_vel,
      #   aes(
      #     x = V1 - v1 / 2,
      #     y = V2 - v2 / 2,
      #     label = sprintf("(%.2f)", mag)
      #   ),
      #   size = 0.8*(annotation_config$size / ggplot2::.pt),  # Small font size
      #   hjust = 0.5,
      #   vjust = 0.5,
      #   alpha = 0.8*annotation_config$alpha
      # )
      # Annotate top-velocity points exactly like notable-point labels
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        if (annotation_config$box) {
          p <- p +
            ggrepel::geom_label_repel(
              data        = top_vel,
              inherit.aes = FALSE,
              aes(x = .data$V1, y = .data$V2, label = .data$name),
              size              = annotation_config$size / ggplot2::.pt,
              color             = annotation_config$color,
              alpha             = annotation_config$alpha,
              fontface          = annotation_config$fontface,
              segment.size      = annotation_config$segment_size,
              segment.alpha     = annotation_config$segment_alpha,
              min.segment.length= annotation_config$min_segment_length,
              max.overlaps      = annotation_config$max_overlaps,
              box.padding       = unit(0.4, "lines"),
              point.padding     = unit(0.3, "lines"),
              force             = 1,
              direction         = "both"
            )
        } else {
          p <- p +
            ggrepel::geom_text_repel(
              data        = top_vel,
              inherit.aes = FALSE,
              aes(x = .data$V1, y = .data$V2, label = .data$name),
              size              = annotation_config$size / ggplot2::.pt,
              color             = annotation_config$color,
              alpha             = annotation_config$alpha,
              fontface          = annotation_config$fontface,
              segment.size      = annotation_config$segment_size,
              segment.alpha     = annotation_config$segment_alpha,
              min.segment.length= annotation_config$min_segment_length,
              max.overlaps      = annotation_config$max_overlaps,
              box.padding       = unit(0.4, "lines"),
              point.padding     = unit(0.3, "lines"),
              force             = 1,
              direction         = "both"
            )
        }
      } else {
        warning("ggrepel package not available - using basic text labels without repulsion")
        p <- p + geom_text(
          data = top_vel,
          aes(x = .data$V1, y = .data$V2, label = .data$name),
          size           = annotation_config$size / ggplot2::.pt,
          color          = annotation_config$color,
          alpha          = annotation_config$alpha,
          fontface       = annotation_config$fontface,
          nudge_x        = - 0.15,
          nudge_y        = - 0.15,
          check_overlap  = TRUE
        )
      }
    }
  }
  
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
  
  # Save plot if requested
  if (layout_config$save_plot) {
    if (draw_arrows){
      if (!is.null(phylo_tree)) {
        filename <- sprintf(
          "clustered_map_ndim_%d_s_t_%g_s_x_%g_cladeDepth_%g_arrowthresh_%g.%s",
          ndim,
          sigma_t,
          sigma_x,
          clade_node_depth,
          layout_config$arrow_plot_threshold,
          layout_config$save_format
        )
      } else {
        filename <- sprintf(
          "clustered_map_ndim_%d_s_t_%g_s_x_%g_arrowthresh_%g.%s",
          ndim,
          sigma_t,
          sigma_x,
          layout_config$arrow_plot_threshold,
          layout_config$save_format
        )
      }
    } else {
      filename <- sprintf(
        "clustered_map_ndim_%d.%s",
        ndim,
        layout_config$save_format
      )
    }
    save_plot(p, filename, layout_config, output_dir)
  }
  
  # Return both plot and velocity data
  if (!is.null(velocity_data)) {
    return(list(
      plot = p,
      velocity_data = velocity_data
    ))
  } else {
    return(p)
  }
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
#' @param output_dir Character. Directory for output files. Required if `interactive` is `FALSE`.
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
#' Note: This function requires the rgl package and OpenGL support. If rgl is not 
#' available, the function will return a 2D plot with a message explaining how to
#' enable 3D visualization.
#'
#' @return Invisibly returns the rgl scene ID for further manipulation if rgl is 
#'         available, or a 2D ggplot object as a fallback.
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' data <- data.frame(
#'   V1 = rnorm(100), V2 = rnorm(100), V3 = rnorm(100), V4 = rnorm(100), name = 1:100,
#'   antigen = rep(c(0,1), 50), antiserum = rep(c(1,0), 50),
#'   cluster = rep(1:5, each=20), year = rep(2000:2009, each=10)
#' )
#'
#' # Create a static plot and save to a temporary file
#' # This example requires an interactive session and the 'rgl' package.
#' if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
#'   temp_dir <- tempdir()
#'   # Basic interactive plot (will open a new window)
#'   if(interactive()) {
#'     plot_3d_mapping(data, ndim=4)
#'   }
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
#'   # Create customized static plot and save it
#' plot_3d_mapping(data, ndim=4,
#'   aesthetic_config = aesthetic_config,
#'   layout_config = layout_config,
#'   interactive = FALSE, output_dir = temp_dir
#' )
#'   list.files(temp_dir)
#'   unlink(temp_dir, recursive = TRUE)
#' }
#'
#' @seealso 
#' \code{\link{plot_temporal_mapping}} for 2D temporal visualization
#' \code{\link{plot_cluster_mapping}} for 2D cluster visualization
#' \code{\link{make_interactive}} for converting 2D plots to interactive versions
#' @importFrom grDevices colorRampPalette
#' @export
plot_3d_mapping <- function(df, ndim, 
                            dim_config = new_dim_reduction_config(),
                            aesthetic_config = new_aesthetic_config(),
                            layout_config = new_layout_config(),
                            interactive = TRUE,
                            output_dir) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("rgl package is required for 3D plotting. Install with: install.packages('rgl')")
  }
  # Check if rgl package is available
  has_rgl <- requireNamespace("rgl", quietly = TRUE)
  
  if (!has_rgl) {
    message("3D visualization requires the 'rgl' package with OpenGL support.")
    message("Install rgl with: install.packages('rgl')")
    message("Falling back to 2D visualization.")
    
    # Create a 2D plot as fallback
    return(plot_temporal_mapping(df, ndim, dim_config, aesthetic_config, layout_config, output_dir = if(!missing(output_dir)) output_dir else NULL))
  }

  if (!interactive && missing(output_dir)) {
    stop("An 'output_dir' must be provided when 'interactive' is FALSE.", call. = FALSE)
  }

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
    rgl::open3d(windowRect = c(100, 100, 
                          100 + layout_config$width * 100,
                          100 + layout_config$height * 100))
  }
  
  # Set point colors
  if("cluster" %in% names(reduced_df)) {
    n_clusters <- length(unique(reduced_df$cluster))
    
    # Use the color palette directly from the aesthetic config, truncating if necessary.
    colors <- aesthetic_config$color_palette[1:min(n_clusters, length(aesthetic_config$color_palette))]
    
    # Handle cases with more clusters than colors by recycling
    if (n_clusters > length(colors)) {
        warning("More clusters than available colors. Colors will be recycled.")
        colors <- rep_len(colors, n_clusters)
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
  rgl::open3d(windowRect = c(100, 100, 1200, 1200))
  rgl::plot3d(reduced_df$dim1, reduced_df$dim2, reduced_df$dim3,
         col = point_colors,
         size = aesthetic_config$point_size * 0.8, # Adjust for 3D
         type = "s",  # spheres
         alpha = aesthetic_config$point_alpha)
  
  # Add axes and labels
  if(layout_config$show_axis) {
    rgl::axes3d(edges = "bbox",
           labels = TRUE,
           tick = TRUE,
           box = TRUE)
  }
  
  # Save if not interactive
  if(!interactive) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    filename <- sprintf("3d_mapping_ndim_%d.%s", 
                        ndim, layout_config$save_format)
    full_path <- file.path(output_dir, filename)
    
    rgl::rgl.snapshot(filename = full_path)
    rgl::close3d()
    return(invisible(full_path))
  }
  
  # Return scene ID invisibly
  invisible(rgl::rgl.cur())
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
#' @param output_dir Character. Directory for output files. This argument is required.
#' 
#' @return No return value, called for side effects (saves a plot to a file).
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
#' # The sole purpose of save_plot is to write a file, so its example must demonstrate this. 
#' # For CRAN tests we wrap the example in \donttest{} to avoid writing files.
#' \donttest{
#' # Create a temporary directory for saving all plots
#' temp_dir <- tempdir()
#' 
#' # --- Example 1: Basic ggplot save ---
#' # Create sample data with 3 dimensions to support both 2D and 3D plots
#' data <- data.frame(
#'   V1 = rnorm(10), V2 = rnorm(10), V3 = rnorm(10), name=1:10,
#'   antigen = rep(c(0,1), 5), antiserum = rep(c(1,0), 5),
#'   year = 2000:2009
#' )
#' p <- plot_temporal_mapping(data, ndim=2)
#' save_plot(p, "temporal_plot.png", output_dir = temp_dir)
#'
#' # --- Example 2: Save with custom layout ---
#' layout_config <- new_layout_config(
#'   width = 12,
#'   height = 8,
#'   dpi = 600,
#'   save_format = "pdf"
#' )
#' save_plot(p, "high_res_plot.pdf", layout_config, output_dir = temp_dir)
#' 
#' # --- Verify files and clean up ---
#' list.files(temp_dir)
#' unlink(temp_dir, recursive = TRUE)
#' }
#' 
#' @importFrom tools file_ext
#' @importFrom ggplot2 ggsave
#' @export
save_plot <- function(plot, filename, layout_config = new_layout_config(),
                      output_dir) {
  if (missing(output_dir)) {
    stop("An 'output_dir' must be provided to save the plot.", call. = FALSE)
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
  if(ggplot2::is_ggplot(plot)) {
    ggsave_white_bg(filename = full_path,
           plot = plot,
           width = layout_config$width,
           height = layout_config$height,
           dpi = layout_config$dpi)
  } else if(inherits(plot, c("arrangement", "grob"))) {
    # For arranged plots from gridExtra
    ggsave_white_bg(filename = full_path,
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
#' @return A \code{plotly} object with interactive features.
#'
#' @examples
#' if (interactive() && requireNamespace("plotly", quietly = TRUE)) {
#' # Create sample data and plot
#' data <- data.frame(
#'   V1 = rnorm(100), V2 = rnorm(100), name=1:100,
#'   antigen = rep(c(0,1), 50), antiserum = rep(c(1,0), 50),
#'   year = rep(2000:2009, each=10), cluster = rep(1:5, each=20)
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


# Newed
#' Create Diagnostic Plots for Multiple Sampling Chains
#'
#' @description
#' Creates trace and density plots for multiple sampling or optimization chains to help
#' assess convergence and mixing. It displays parameter trajectories and their
#' distributions across all chains.
#'
#' @param chain_files A character vector of paths to CSV files, where each file contains data for one chain.
#' @param mutual_size Integer. The number of samples to use from the end of each chain for plotting.
#' @param output_file Character. The path for saving the plot. Required if `save_plot` is TRUE.
#' @param output_dir Character. The directory for saving output files. Required if `save_plot` is TRUE.
#' @param save_plot Logical. If TRUE, saves the plot to a file. Default: FALSE.
#' @param width,height,dpi Numeric. The dimensions and resolution for the saved plot.
#' @return A `ggplot` object of the combined plots.
#' @examples
#' # This example uses sample data files that would be included with the package.
#' chain_files <- c(
#'   system.file("extdata", "diag_chain1.csv", package = "topolow"),
#'   system.file("extdata", "diag_chain2.csv", package = "topolow"),
#'   system.file("extdata", "diag_chain3.csv", package = "topolow")
#' )
#'
#' # Only run the example if the files are found
#' if (all(nzchar(chain_files))) {
#'   # Create diagnostic plot without saving to a file
#'   plot_mcmc_diagnostics(chain_files, mutual_size = 50, save_plot = FALSE)
#' }
#'
#' @export
plot_mcmc_diagnostics <- function(chain_files,
                                    mutual_size = 2000,
                                    output_file = "diagnostic_plots.png",
                                    output_dir,
                                    save_plot = FALSE,
                                    width = 3000, height = 3000, dpi = 300) {
                                      # Check if gridextra is available
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package is required for plotting. Please install with install.packages('gridExtra').")
  }
  if (save_plot && (missing(output_dir) || missing(output_file))) {
    stop("`output_dir` and `output_file` must be provided when `save_plot` is TRUE.", call. = FALSE)
  }

  # Read and process chain data
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  chains <- lapply(chain_files, function(file) {
    df <- utils::read.csv(file)
    # Take the last `mutual_size` samples
    tail(df[, par_names], mutual_size)
  })

  n_params <- ncol(chains[[1]])
  n_chains <- length(chains)

  # Create a list to hold all the individual plots
  plot_list <- list()

  for (i in seq_len(n_params)) {
    # Combine data from all chains for the current parameter
    trace_data <- do.call(rbind, lapply(seq_len(n_chains), function(j) {
      data.frame(
        Chain = as.factor(j),
        Iteration = seq_len(nrow(chains[[j]])),
        Value = chains[[j]][, i]
      )
    }))

    # Create Trace Plot
    p_trace <- ggplot2::ggplot(trace_data,
                               ggplot2::aes(x = .data$Iteration, y = .data$Value, color = .data$Chain)) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::labs(title = paste("Trace Plot:", par_names[i]), x = "Iteration", y = "Value") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA)
      )
    plot_list <- c(plot_list, list(p_trace))

    # Create Density Plot
    p_density <- ggplot2::ggplot(trace_data, ggplot2::aes(x = .data$Value, color = .data$Chain)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::labs(title = paste("Density Plot:", par_names[i]), x = "Value", y = "Density") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA)
      )
    plot_list <- c(plot_list, list(p_density))
  }

  # Arrange all plots into a grid
  combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 2)

  # Optionally save the combined plot
  if (save_plot) {
    full_output_path <- file.path(output_dir, output_file)
    ggsave_white_bg(full_output_path, combined_plot,
                    width = width / dpi, height = height / dpi,
                    dpi = dpi, limitsize = FALSE)
  }

  return(combined_plot)
}


# Newed
#' Plot Method for topolow Convergence Diagnostics
#'
#' @description
#' Creates visualizations of convergence diagnostics from a sampling run, including
#' parameter mean trajectories and covariance matrix stability over iterations. This helps
#' assess whether parameter estimation has converged.
#'
#' @details
#' The function generates two types of plots:
#' 1. Parameter mean plots: Shows how the mean value for each parameter changes over iterations.
#'    Stabilization of these plots indicates convergence.
#' 2. Covariance change plot: Shows the relative change in the Frobenius norm of the
#'    covariance matrix. A decreasing trend approaching zero indicates stable relationships
#'    between parameters.
#'
#' @param x A `topolow_convergence` object from `check_gaussian_convergence()`.
#' @param param_names Optional character vector of parameter names for plot titles.
#'   If NULL, names are taken from the input object.
#' @param ... Additional arguments (not currently used).
#' @return A grid of plots showing convergence metrics.
#' @examples
#' # Example with simulated data
#' chain_data <- data.frame(
#'   param1 = rnorm(1000, mean = 1.5, sd = 0.1),
#'   param2 = rnorm(1000, mean = -0.5, sd = 0.2)
#' )
#'
#' # Check convergence
#' results <- check_gaussian_convergence(chain_data)
#'
#' # Plot diagnostics
#' plot(results)
#'
#' # With custom parameter names
#' plot(results, param_names = c("Parameter 1 (log)", "Parameter 2 (log)"))
#'
#' @seealso \code{\link{check_gaussian_convergence}} for generating the convergence object.
#' @method plot topolow_convergence
#' @export
plot.topolow_convergence <- function(x, param_names = NULL, ...) {
  # Check if required packages are available
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package is required for plotting. Please install with install.packages('gridExtra').")
  }
  # Use param_names from the object if not provided by the user
  if (is.null(param_names)) {
    param_names <- x$param_names
  }

  # Create a plot for the mean trajectory of each parameter
  mean_plots <- lapply(seq_along(param_names), function(i) {
    plot_data <- data.frame(
      Iteration = seq_len(nrow(x$mean_history)),
      Value = x$mean_history[, i]
    )
    ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Iteration, y = .data$Value)) +
      ggplot2::geom_line(color = "steelblue") +
      ggplot2::labs(title = paste("Parameter Mean:", param_names[i]), x = "Iteration", y = "Mean Value") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA)
      )
  })

  # Create a plot for the change in the covariance matrix norm
  cov_plot_data <- data.frame(
    Iteration = seq_along(x$cov_changes),
    Change = x$cov_changes
  )
  cov_plot <- ggplot2::ggplot(cov_plot_data, ggplot2::aes(x = .data$Iteration, y = .data$Change)) +
    ggplot2::geom_line(color = "darkred") +
    ggplot2::labs(title = "Covariance Matrix Stability", x = "Iteration", y = "Relative Change in Norm") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )

  # Arrange all plots into a single grid
  gridExtra::grid.arrange(
    grobs = c(mean_plots, list(cov_plot)),
    ncol = 2
  )
}

#' Print Method for topolow Convergence Diagnostics
#'
#' @param x A `topolow_convergence` object.
#' @param ... Additional arguments passed to print.
#' @return No return value; called for its side effect of printing a summary.
#' @method print topolow_convergence
#' @export
print.topolow_convergence <- function(x, ...) {
  cat("--- topolow Convergence Diagnostics ---\n")
  cat(sprintf("Overall Convergence Achieved: %s\n", x$converged))
  cat(sprintf("Mean Vector Converged: %s\n", x$mean_converged))
  cat(sprintf("Covariance Matrix Converged: %s\n", x$cov_converged))
  cat("\nFinal Parameter Means:\n")
  print(setNames(x$final_mean, x$param_names))
}

# Newed
#' Plot Method for topolow parameter estimation Diagnostics
#'
#' @description
#' Creates trace and density plots for multiple chains to assess convergence and mixing.
#' This is an S3 method that dispatches on `topolow_diagnostics` objects.
#'
#' @param x A `topolow_diagnostics` object from `calculate_diagnostics()`.
#' @param output_dir Character. Directory for output files. Required if `save_plot` is TRUE.
#' @param output_file Character path for saving the plot.
#' @param save_plot Logical. Whether to save the plot.
#' @param ... Additional arguments passed to `create_diagnostic_plots`.
#' @return A ggplot object of the combined plots.
#' @method plot topolow_diagnostics
#' @export
plot.topolow_diagnostics <- function(x,
                                        output_dir,
                                        output_file = "topolow_param_diagnostics.png",
                                        save_plot = FALSE,
                                        ...) {
  # This method is a wrapper around the main plotting function
  plot_mcmc_diagnostics(chain_files = x$chains, # Pass the actual chain data
                          mutual_size = x$mutual_size,
                          output_file = output_file,
                          output_dir = if (!missing(output_dir)) output_dir else NULL,
                          save_plot = save_plot,
                          ...)
}

#' Print Method for topolow parameter estimation Diagnostics
#'
#' @param x A `topolow_diagnostics` object.
#' @param ... Additional arguments passed to print.
#' @return No return value; called for its side effect of printing a summary.
#' @method print topolow_diagnostics
#' @export
print.topolow_diagnostics <- function(x, ...) {
  cat("--- topolow Adaptive Sampling Diagnostics ---\n")
  cat("\nR-hat values (should be close to 1 for convergence):\n")
  print(setNames(x$rhat, x$param_names))
  cat("\nEffective Sample Sizes (higher is better):\n")
  print(setNames(x$ess, x$param_names))
}

# Newed

# New:
#' Plot Network Structure
#'
#' @description
#' Creates a visualization of the dissimilarity matrix as a network graph, showing
#' data availability patterns and connectivity between points.
#'
#' @param network_results The list output from `analyze_network_structure()`.
#' @param output_file Character. An optional full path to save the plot. If NULL, the plot is not saved.
#' @param width Numeric. Width in pixels for saved plot (default: 3000).
#' @param height Numeric. Height in pixels for saved plot (default: 3000).
#' @param dpi Numeric. Resolution in dots per inch (default: 300).
#'
#' @return A `ggplot` object representing the network graph.
#'
#' @examples
#' # Create a sample dissimilarity matrix
#' adj_mat <- matrix(runif(25), 5, 5)
#' rownames(adj_mat) <- colnames(adj_mat) <- paste0("Point", 1:5)
#' adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
#' diag(adj_mat) <- 0
#' net_analysis <- analyze_network_structure(adj_mat)
#'
#' # Create and display the plot
#' plot_network_structure(net_analysis)
#'
#' @importFrom ggplot2 ggplot geom_segment geom_point coord_fixed theme_void theme element_text labs ggsave
#' @export
plot_network_structure <- function(network_results, output_file = NULL,
                                   width = 3000, height = 3000, dpi = 300) {
                                    # Check if required packages are available
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph packages is required for this function. Please install it with: install.packages('igraph')")
  }

  if (nrow(network_results$adjacency) == 0 || sum(network_results$adjacency) == 0) {
    # Handle empty network case
    layout_df <- data.frame(
      x = numeric(0),
      y = numeric(0),
      node = character(0)
    )
    edges <- data.frame(
      from = character(0),
      to = character(0),
      x.from = numeric(0),
      y.from = numeric(0),
      x.to = numeric(0),
      y.to = numeric(0)
    )
  } else {
    # Existing logic for non-empty networks
    graph <- igraph::graph_from_adjacency_matrix(
      network_results$adjacency,
      mode = "undirected"
    )

    # ALSO ADD: Check if graph has edges
    if (igraph::ecount(graph) == 0) {
      # No edges - create simple layout
      layout <- matrix(c(seq_along(rownames(network_results$adjacency)),
                         rep(0, nrow(network_results$adjacency))),
                       ncol = 2)
    } else {
      # Use the Fruchterman-Reingold layout algorithm for node positions
      layout <- igraph::layout_with_fr(graph)
    }

    layout_df <- data.frame(
      x = layout[,1],
      y = layout[,2],
      node = rownames(network_results$adjacency)
    )

    # Handle edges
    if (igraph::ecount(graph) > 0) {
      edges <- as.data.frame(igraph::as_edgelist(graph))
      names(edges) <- c("from", "to")
      edges <- merge(edges, layout_df, by.x = "from", by.y = "node", all.x = TRUE)
      edges <- merge(edges, layout_df, by.x = "to", by.y = "node", all.x = TRUE, suffixes = c(".from", ".to"))
    } else {
      edges <- data.frame(
        from = character(0), to = character(0),
        x.from = numeric(0), y.from = numeric(0),
        x.to = numeric(0), y.to = numeric(0)
      )
    }
  }


  # Create the plot with fixed aesthetic values
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x.from, y = .data$y.from, xend = .data$x.to, yend = .data$y.to),
      alpha = 0.3,
      color = "grey50"
    ) +
    ggplot2::geom_point(
      data = layout_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "steelblue",
      size = 4,        # Fixed value instead of aesthetic_config$point_size
      alpha = 0.8      # Fixed value instead of aesthetic_config$point_alpha
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 14,    # Fixed value instead of aesthetic_config$title_size
        hjust = 0.5
      ),
      legend.position = "none"  # Fixed to no legend
    ) +
    ggplot2::labs(
      title = sprintf(
        "Network Structure (%.1f%% Complete)",
        network_results$summary$completeness * 100
      )
    )

  # Save the plot if an output file path is provided
  if (!is.null(output_file)) {
    ggsave_white_bg(
      filename = output_file,
      plot = p,
      width = width / dpi,
      height = height / dpi,
      dpi = dpi
    )
  }

  return(p)
}


#' Plot Fitted vs. True Dissimilarities
#'
#' @description
#' Creates diagnostic plots comparing fitted dissimilarities from a model against the true
#' dissimilarities. It generates both a scatter plot with an identity line and
#' prediction intervals, and a residuals plot.
#'
#' @param dissimilarity_matrix Matrix of true dissimilarities.
#' @param p_dissimilarity_mat Matrix of predicted/fitted dissimilarities.
#' @param scenario_name Character string for output file naming. Used if `save_plot` is TRUE.
#' @param ndim Integer number of dimensions used in the model.
#' @param confidence_level Numeric confidence level for prediction intervals (default: 0.95).
#' @param save_plot Logical. Whether to save plots to files. Default: FALSE.
#' @param output_dir Character. Directory for output files. Required if `save_plot` is TRUE.
#'
#' @return A list containing the `scatter_plot` and `residuals_plot` ggplot objects.
#'
#' @examples
#' # Create sample data
#' true_dist <- matrix(runif(100, 1, 10), 10, 10)
#' pred_dist <- true_dist + rnorm(100)
#'
#' # Create plots without saving
#' plots <- scatterplot_fitted_vs_true(
#'   dissimilarity_matrix = true_dist,
#'   p_dissimilarity_mat = pred_dist,
#'   save_plot = FALSE
#'  )
#'
#' # You can then display a plot, for instance:
#' # plots$scatter_plot
#'
#' @importFrom stats na.omit cor lm coef
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_smooth geom_hline scale_y_continuous
#'   annotate labs theme_classic theme element_text element_blank element_rect scale_x_continuous
#' @export
scatterplot_fitted_vs_true <- function(dissimilarity_matrix, p_dissimilarity_mat,
                                       scenario_name, ndim,
                                       save_plot = FALSE,
                                       output_dir,
                                       confidence_level = 0.95) {

  if (save_plot && (missing(output_dir) || missing(scenario_name) || missing(ndim))) {
    stop("`output_dir`, `scenario_name`, and `ndim` must be provided when `save_plot` is TRUE.", call. = FALSE)
  }

  # Create a data frame for evaluation, removing NA values
  # Extract numeric values from threshold indicators
  true_dissim_vector <- extract_numeric_values(as.vector(dissimilarity_matrix))
  pred_dissim_vector <- as.vector(p_dissimilarity_mat)

  eval_df <- data.frame(
    true_dissimilarity = true_dissim_vector,
    pred_dissimilarity = pred_dissim_vector
  )
  eval_df <- stats::na.omit(eval_df)

  # Calculate performance metrics
  mae <- mean(abs(eval_df$true_dissimilarity - eval_df$pred_dissimilarity))
  pearson_corr <- stats::cor(eval_df$pred_dissimilarity, eval_df$true_dissimilarity, method = "pearson")

  # Calculate prediction interval margin of error
  margin_of_error <- calculate_prediction_interval(dissimilarity_matrix, p_dissimilarity_mat,
                                                   confidence_level = confidence_level)

  # Get regression line parameters
  reg_line <- stats::lm(pred_dissimilarity ~ true_dissimilarity, data = eval_df)
  slope <- stats::coef(reg_line)[2]
  intercept <- stats::coef(reg_line)[1]

  # --- Create Scatter Plot ---
  scatter_plot <- ggplot2::ggplot(eval_df, ggplot2::aes(x = .data$true_dissimilarity, y = .data$pred_dissimilarity)) +
    ggplot2::geom_point(alpha = 0.6, color = "#3366FF", size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed",
                         level = 0.99, fill = "pink", alpha = 0.5, linewidth = 0.5) +
    ggplot2::geom_abline(intercept = intercept + margin_of_error, slope = slope, linetype = "dotted", color = "black", linewidth = 0.7) +
    ggplot2::geom_abline(intercept = intercept - margin_of_error, slope = slope, linetype = "dotted", color = "black", linewidth = 0.7) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 9),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1)
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max(eval_df$true_dissimilarity, na.rm=T))) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(eval_df$pred_dissimilarity, na.rm=T))) +
    ggplot2::labs(
      x = "True Dissimilarity",
      y = "Estimated Dissimilarity",
      title = sprintf("MAE = %.2f, Pearson's r = %.2f, Prediction Interval = %.2f",
                      mae, pearson_corr, margin_of_error)
    )

  # --- Create Residuals Plot ---
  eval_df$residuals <- eval_df$true_dissimilarity - eval_df$pred_dissimilarity
  residuals_plot <- ggplot2::ggplot(eval_df, ggplot2::aes(x = .data$pred_dissimilarity, y = .data$residuals)) +
    ggplot2::geom_point(alpha = 0.6, color = "darkred") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Fitted Values", y = "Residuals", title = "Residuals vs. Fitted Values") +
    ggplot2::annotate("text", x = max(eval_df$pred_dissimilarity) * 0.7, y = Inf, vjust = 1.5,
                      label = sprintf("Mean absolute error = %.2f", mean(abs(eval_df$residuals))),
                      color = "blue", size = 3)

  # --- Save Plots if Requested ---
  if (save_plot) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    # Save scatter plot
    scatter_file <- file.path(output_dir, sprintf("%s_prediction_scatter_dim_%d.png", scenario_name, ndim))
    ggsave_white_bg(scatter_file, scatter_plot, width = 8, height = 8, dpi = 300)

    # Save residuals plot
    residuals_file <- file.path(output_dir, sprintf("%s_residuals_vs_fitted_dim_%d.png", scenario_name, ndim))
    ggsave_white_bg(residuals_file, residuals_plot, width = 8, height = 8, dpi = 300)
  }

  return(list(
    scatter_plot = scatter_plot,
    residuals_plot = residuals_plot
  ))
}
