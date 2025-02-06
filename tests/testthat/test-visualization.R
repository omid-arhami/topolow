test_that("plot_temporal_mapping creates valid plot", {
  # Create test data
  test_df <- data.frame(
    V1 = rnorm(10),
    V2 = rnorm(10),
    antigen = rep(c(TRUE, FALSE), 5),
    antiserum = rep(c(FALSE, TRUE), 5),
    year = 2000:2009
  )
  
  plot <- plot_temporal_mapping(test_df, ndim=2)
  expect_s3_class(plot, "ggplot")
  
  # Test configuration objects
  aesthetic_config <- new_aesthetic_config()
  layout_config <- new_layout_config()
  plot <- plot_temporal_mapping(
    test_df, ndim=2,
    aesthetic_config = aesthetic_config,
    layout_config = layout_config
  )
  expect_s3_class(plot, "ggplot")
})

test_that("plot_cluster_mapping creates valid plot", {
  test_df <- data.frame(
    V1 = rnorm(10),
    V2 = rnorm(10),
    antigen = rep(c(TRUE, FALSE), 5),
    antiserum = rep(c(FALSE, TRUE), 5),
    cluster = rep(c("A", "B"), 5)
  )
  
  plot <- plot_cluster_mapping(test_df, ndim=2)
  expect_s3_class(plot, "ggplot")
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
