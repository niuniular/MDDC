test_that("validate find optimal coef results", {
  contin_table <- matrix(c(10, 5, 3, 7), nrow = 2)

  results1 <- find_optimal_coef(contin_table, exclude_small_count = TRUE)
  results2 <- find_optimal_coef(contin_table, exclude_small_count = FALSE)

  expect_type(results1, "list")
  expect_type(results2, "list")
  expect_length(results1$coef, ncol(contin_table))
  expect_length(results1$FDR, ncol(contin_table))
  expect_length(results2$coef, ncol(contin_table))
  expect_length(results2$FDR, ncol(contin_table))
})
