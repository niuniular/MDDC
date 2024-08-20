library(testthat)
library(ggplot2)

test_that("plot_heatmap generates a ggplot object", {
  # Create example data
  data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(data) <- paste0("Row", 1:10)
  colnames(data) <- paste0("Col", 1:10)

  # Generate the plot
  plot <- plot_heatmap(data)

  # Check that the result is a ggplot object
  expect_s3_class(plot, "gg")
})
