test_that("test MDDC MC function", {
  set.seed(42)
  contin_table <- matrix(rpois(30 * 2, 10), nrow = 12)

  test_results <- function(col_specific_cutoff, separate, if_col_cor) {
    mddc_mc_res <- mddc_mc(contin_table,
      col_specific_cutoff = col_specific_cutoff,
      separate = separate,
      if_col_cor = if_col_cor
    )

    expect_type(mddc_mc_res, "list")
    expect_equal(dim(contin_table), dim(mddc_mc_res$mc_pval))
    expect_equal(dim(contin_table), dim(mddc_mc_res$mc_signal))
    expect_equal(dim(contin_table), dim(mddc_mc_res$fisher_signal))
    expect_equal(dim(contin_table), dim(mddc_mc_res$corr_signal_pval))
    expect_equal(dim(contin_table), dim(mddc_mc_res$corr_signal_adj_pval))
  }

  # Check different parameter combinations
  test_results(TRUE, TRUE, TRUE)
  test_results(FALSE, FALSE, FALSE)
  test_results(TRUE, FALSE, FALSE)
  test_results(FALSE, TRUE, FALSE)
})
