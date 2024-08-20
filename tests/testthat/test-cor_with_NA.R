library(testthat)

test_that("cor_with_NA computes correlations correctly", {
  # Test matrix with no missing values
  mat1 <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), ncol = 3)
  expected_cor1 <- cor(mat1)
  diag(expected_cor1) <- NA
  result1 <- cor_with_NA(mat1, TRUE)
  expect_equal(result1, expected_cor1)

  # Test matrix where all values are NA
  mat3 <- matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA), ncol = 3)
  expected_cor3 <- matrix(NA, nrow = 3, ncol = 3)
  diag(expected_cor3) <- NA
  row.names(expected_cor3) <- colnames(mat3)
  colnames(expected_cor3) <- colnames(mat3)
  result3 <- cor_with_NA(mat3, TRUE)
  expect_equal(result3, expected_cor3)
})
