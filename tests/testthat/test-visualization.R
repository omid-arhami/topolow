# tests/testthat/test-visualization.R

# Copyright (c) 2025 Omid Arhami omid.arhami@uga.edu

test_that("scatterplot_fitted_vs_true creates valid plots", {
  true_dist <- matrix(runif(100, 1, 10), 10, 10)
  pred_dist <- true_dist + rnorm(100)

  plots <- scatterplot_fitted_vs_true(
    dissimilarity_matrix = true_dist,
    p_dissimilarity_mat = pred_dist,
    save_plot = FALSE
  )

  expect_true(ggplot2::is_ggplot(plots$scatter_plot))
  expect_true(ggplot2::is_ggplot(plots$residuals_plot))
})

test_that("plot_network_structure creates a valid plot", {
  adj_mat <- matrix(runif(25), 5, 5)
  diag(adj_mat) <- 0
  net_analysis <- analyze_network_structure(adj_mat)

  p <- plot_network_structure(net_analysis)
  expect_true(ggplot2::is_ggplot(p))
})

test_that("create_diagnostic_plots works with temporary files", {
  # Create dummy chain files in a temporary directory
  temp_dir <- tempdir()
  chain_files <- character(2)
  par_names <- c("log_N", "log_k0", "log_cooling_rate", "log_c_repulsion")
  sample_data <- data.frame(
    log_N = rnorm(100), log_k0 = rnorm(100),
    log_cooling_rate = rnorm(100), log_c_repulsion = rnorm(100)
  )
  for (i in 1:2) {
    chain_files[i] <- file.path(temp_dir, paste0("chain", i, ".csv"))
    write.csv(sample_data, chain_files[i], row.names = FALSE)
  }

  p <- create_diagnostic_plots(chain_files, mutual_size = 50, save_plot = FALSE)
  # grid.arrange returns a gtable object
  expect_s3_class(p, "gtable")

  unlink(chain_files)
})

test_that("plot_temporal_mapping creates valid plot", {
  # Create test data
  test_df <- data.frame(
    V1 = rnorm(10),
    V2 = rnorm(10),
    antigen = rep(c(TRUE, FALSE), 5),
    antiserum = rep(c(FALSE, TRUE), 5),
    name = paste0(ifelse(rep(c(TRUE, FALSE), 5), "V/", "S/"), "strain", 1:10),
    year = 2000:2009
  )
  
  plot <- plot_temporal_mapping(test_df, ndim=2)
  expect_true(ggplot2::is_ggplot(plot))
  
  # Test configuration objects
  aesthetic_config <- new_aesthetic_config()
  layout_config <- new_layout_config()
  plot <- plot_temporal_mapping(
    test_df, ndim=2,
    aesthetic_config = aesthetic_config,
    layout_config = layout_config
  )
  expect_true(ggplot2::is_ggplot(plot))
})

test_that("plot_cluster_mapping creates valid plot", {
  test_df <- data.frame(
    V1 = rnorm(10),
    V2 = rnorm(10),
    antigen = rep(c(TRUE, FALSE), 5),
    antiserum = rep(c(FALSE, TRUE), 5),
    cluster = rep(c("A", "B"), 5),
    name = paste0(ifelse(rep(c(TRUE, FALSE), 5), "V/", "S/"), "strain", 1:10)
  )
  
  plot <- plot_cluster_mapping(test_df, ndim=2, )
  expect_true(ggplot2::is_ggplot(plot))
})

test_that("plot configuration objects work correctly", {
  config <- new_aesthetic_config(
    point_size = 2,
    point_alpha = 0.5
  )
  expect_s3_class(config, "aesthetic_config")
  
  config <- new_layout_config(
    width = 8,
    height = 6
  )
  expect_s3_class(config, "layout_config")
})
