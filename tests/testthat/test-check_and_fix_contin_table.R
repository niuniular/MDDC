library(testthat)

test_that("Valid data matrix with row and column names", {
  # Create a valid matrix with row and column names
  valid_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2)
  rownames(valid_matrix) <- c("AE_1", "AE_2")
  colnames(valid_matrix) <- c("drug_1", "drug_2", "drug_3", "drug_4")

  result <- check_and_fix_contin_table(valid_matrix)
  expect_equal(result, valid_matrix)
  expect_equal(rownames(result), rownames(valid_matrix))
  expect_equal(colnames(result), colnames(valid_matrix))
})

test_that("Data frame input", {
  # Create a data frame and convert it to matrix
  df <- data.frame(Drug_1 = c(1, 2), Drug_2 = c(3, 4))
  rownames(df) <- c("AE_1", "AE_2")
  result <- as.data.frame(check_and_fix_contin_table(df))
  expect_equal(result, df)
})

test_that("Matrix without row and column names", {
  # Create a matrix without row and column names
  no_names_matrix <- matrix(1:6, nrow = 2)

  result <- check_and_fix_contin_table(no_names_matrix)
  expected_row_names <- paste0("AE_", seq_len(nrow(no_names_matrix)))
  expected_col_names <- paste0("drug_", seq_len(ncol(no_names_matrix)))

  expect_equal(rownames(result), expected_row_names)
  expect_equal(colnames(result), expected_col_names)
})

test_that("Matrix with negative values", {
  # Create a matrix with negative values
  negative_matrix <- matrix(c(-1, 2, 3, 4), nrow = 2)

  expect_error(
    check_and_fix_contin_table(negative_matrix),
    "Input is not a data matrix with non-negative integers."
  )
})

test_that("Matrix with non-integer values", {
  # Create a matrix with non-integer values
  non_integer_matrix <- matrix(c(1.5, 2, 3, 4), nrow = 2)

  expect_error(
    check_and_fix_contin_table(non_integer_matrix),
    "Input is not a data matrix with non-negative integers."
  )
})
