test_that("correlation_matrix creates correct matrix", {
  # Test with 3x3 matrix and rho = 0.5
  result <- correlation_matrix(3, 0.5)
  expected <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
  expect_equal(result, expected)
  
  # Test with 4x4 matrix and rho = 0.8
  result <- correlation_matrix(4, 0.8)
  expected <- matrix(c(1, 0.8, 0.8, 0.8, 0.8, 1, 0.8, 0.8, 0.8, 0.8, 1, 0.8, 0.8, 0.8, 0.8, 1), nrow = 4)
  expect_equal(result, expected)
})

