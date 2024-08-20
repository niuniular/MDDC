test_that("get_log_bootstrap_cutoff computes cutoffs correctly", {
  # Mocked contingency table
  contin_table <- matrix(c(10, 5, 3, 7), nrow = 2)

  # Test with specific quantile and replication count
  result <- get_log_bootstrap_cutoff(contin_table, 0.95, 100, seed = 123)

  # Check structure of the result
  expect_type(result, "list")
  expect_named(result, c("cutoffs", "null_dist_s"))
  expect_length(result$cutoffs, ncol(contin_table))
  expect_true(is.list(result))
})
