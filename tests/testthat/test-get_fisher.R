test_that("get_fisher computes p-value correctly", {
  # Mocked contingency table
  contin_table <- matrix(c(10, 5, 3, 7), nrow = 2)

  # Test Fisher's exact test with specific indices
  p_value1 <- get_fisher(contin_table, 1, 1, FALSE)
  p_value2 <- get_fisher(contin_table, 1, 1, TRUE)

  # Test that p-value is numeric
  expect_type(p_value1, "double")
  expect_type(p_value2, "double")
  expect_length(p_value1, 1)
  expect_length(p_value2, 1)
  expect_true(p_value1 >= 0 && p_value1 <= 1)
  expect_true(p_value2 >= 0 && p_value2 <= 1)
})
