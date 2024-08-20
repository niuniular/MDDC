test_that("test MDDC boxplot function", {
  # Generate a random contingency table
  contin_table <- matrix(rpois(30 * 2, 3), nrow = 12)

  # Test different parameter combinations
  test_results <- function(col_specific_cutoff, separate, if_col_cor) {
    mddc_boxplot_res <- mddc_boxplot(contin_table,
      col_specific_cutoff =
        col_specific_cutoff,
      separate = separate,
      if_col_cor = if_col_cor
    )

    expect_type(mddc_boxplot_res, "list")
    expect_equal(
      dim(contin_table),
      dim(mddc_boxplot_res$boxplot_signal)
    )
    expect_equal(
      dim(contin_table),
      dim(mddc_boxplot_res$corr_signal_pval)
    )
    expect_equal(
      dim(contin_table),
      dim(mddc_boxplot_res$corr_signal_adj_pval)
    )
  }

  # Check different parameter combinations
  test_results(TRUE, TRUE, TRUE)
  test_results(FALSE, FALSE, FALSE)
  test_results(TRUE, FALSE, FALSE)
  test_results(FALSE, TRUE, FALSE)

  # Edge Case: Matrix with only zeros
  zero_table <- matrix(0, nrow = 10, ncol = 5)
  mddc_boxplot_res <- mddc_boxplot(zero_table,
    col_specific_cutoff = TRUE,
    separate = TRUE,
    if_col_cor = TRUE
  )
  expect_type(mddc_boxplot_res, "list")
  expect_equal(dim(zero_table), dim(mddc_boxplot_res$boxplot_signal))
  expect_equal(dim(zero_table), dim(mddc_boxplot_res$corr_signal_pval))
  expect_equal(dim(zero_table), dim(mddc_boxplot_res$corr_signal_adj_pval))

  # Edge Case: Matrix with large values
  large_table <- matrix(rpois(30 * 2, 1000), nrow = 12)
  mddc_boxplot_res <- mddc_boxplot(large_table,
    col_specific_cutoff = TRUE,
    separate = TRUE,
    if_col_cor = TRUE
  )
  expect_type(mddc_boxplot_res, "list")
  expect_equal(dim(large_table), dim(mddc_boxplot_res$boxplot_signal))
  expect_equal(dim(large_table), dim(mddc_boxplot_res$corr_signal_pval))
  expect_equal(dim(large_table), dim(mddc_boxplot_res$corr_signal_adj_pval))
})
