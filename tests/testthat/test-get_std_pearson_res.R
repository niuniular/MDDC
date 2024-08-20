library(testthat)

# Define a test for the get_std_pearson_res function
test_that("get_std_pearson_res works correctly", {
  # Create a 6 by 4 data matrix
  set.seed(42)
  dat_mat <- matrix(rpois(6 * 4, 20), nrow = 6)
  dat_mat <- as.data.frame(dat_mat)
  rownames(dat_mat) <- paste0("Row", 1:6)
  colnames(dat_mat) <- paste0("Col", 1:4)

  # Check the format of the data matrix and assign row and column names
  contin_table <- check_and_fix_contin_table(as.matrix(dat_mat))

  # Compute the standardized Pearson residuals
  std_pearson_res <- get_std_pearson_res(contin_table)

  # Test if the result has the same dimensions as the input
  expect_equal(dim(std_pearson_res), dim(contin_table))

  # Test if the row and column names match
  expect_equal(rownames(std_pearson_res), rownames(contin_table))
  expect_equal(colnames(std_pearson_res), colnames(contin_table))

  # Optionally, test with another size of the matrix
  dat_mat2 <- matrix(rpois(3 * 5, 30), nrow = 3)
  dat_mat2 <- as.data.frame(dat_mat2)
  rownames(dat_mat2) <- paste0("Row", 1:3)
  colnames(dat_mat2) <- paste0("Col", 1:5)

  contin_table2 <- check_and_fix_contin_table(as.matrix(dat_mat2))
  std_pearson_res2 <- get_std_pearson_res(contin_table2)

  # Test if the result has the same dimensions as the input
  expect_equal(dim(std_pearson_res2), dim(contin_table2))

  # Test if the row and column names match
  expect_equal(rownames(std_pearson_res2), rownames(contin_table2))
  expect_equal(colnames(std_pearson_res2), colnames(contin_table2))
})
