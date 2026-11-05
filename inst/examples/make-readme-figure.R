# Copyright (c) 2025 Omid Arhami o.arhami@gmail.com
# inst/examples/make-readme-figure.R
#
# Rebuilds man/figures/h3n2-antigenic-map.png, the H3N2 map shown in README.md.
# Run from the package root. The coordinates are the saved 5D Topolow embedding
# of the Smith et al. (2004) H3N2 HI data produced by
# inst/examples/methods-comparison-h3n2-hiv-denv.Rmd.

library(topolow)
library(ggplot2)

pos <- read.csv("inst/examples/data_files/topolow_H3N2_positions.csv",
                stringsAsFactors = FALSE)
# A handful of strains carry no cluster assignment; drop them
pos <- pos[nzchar(trimws(pos$cluster)), ]

# Chronological cluster order, with each cluster's color taken from the data
ord <- c("HK68", "EN72", "VI75", "TX77", "BK79", "SI87",
         "BE89", "BE92", "WU95", "SY97", "FU02")
pal <- vapply(ord, function(cl) unique(pos$color[pos$cluster == cl])[1],
              character(1))

aesthetic_config <- new_aesthetic_config(
  point_shapes = c(antigen = 16, antiserum = 5),
  point_size = 2.2,
  point_alpha = 0.75,
  color_palette = unname(pal),
  show_legend = TRUE,
  legend_position = "right",
  show_title = FALSE,
  show_labels = FALSE,
  title_size = 12,
  axis_title_size = 11,
  axis_text_size = 9
)

# swap_axes puts the dominant (drift) dimension on the horizontal axis
layout_config <- new_layout_config(
  width = 11,
  height = 3.2,
  dpi = 200,
  reverse_x = 1,
  reverse_y = 1,
  swap_axes = TRUE,
  coord_type = "fixed",
  save_plot = FALSE
)

dim_config <- new_dim_reduction_config(method = "pca", n_components = 2,
                                       scale = FALSE)

p <- plot_cluster_mapping(
  pos,
  ndim = 5,
  dim_config = dim_config,
  aesthetic_config = aesthetic_config,
  layout_config = layout_config,
  cluster_legend_title = "Antigenic cluster",
  cluster_legend_order = ord
)

ggsave("man/figures/h3n2-antigenic-map.png", p,
       width = layout_config$width, height = layout_config$height,
       dpi = layout_config$dpi, bg = "white")
