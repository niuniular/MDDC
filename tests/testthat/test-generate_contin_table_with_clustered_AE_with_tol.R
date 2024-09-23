# Load necessary libraries
library(testthat)
library(MASS) # For mvrnorm function

# Define the test cases
test_that("generate_contin_table_with_clustered_AE_with_tol
          works with valid inputs", {
  row_marginal <- c(10, 20, 30)
  column_marginal <- c(15, 25, 20)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)

  AE_idx <- data.frame(idx = c(1, 1, 2), AE = c("AE1", "AE2", "AE3"))
  rho <- 0.5
  n_rep <- 5
  contin_table <- matrix(c(2, 3, 5, 10, 7, 12, 18, 9, 13), nrow = 3)

  result <- generate_contin_table_with_clustered_AE_with_tol(
    row_marginal = row_marginal,
    column_marginal = column_marginal,
    signal_mat = signal_mat,
    contin_table = contin_table,
    AE_idx = AE_idx,
    n_rep = n_rep,
    rho = rho
  )

  expect_length(result, n_rep)
  expect_true(all(vapply(result, is.matrix, logical(1))))
})

test_that("generate_contin_table_with_clustered_AE_with_tol errors if row
          and column marginals don't match", {
  row_marginal <- c(10, 20)
  column_marginal <- c(15, 25, 30)

  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(row_marginal,
      column_marginal,
      signal_mat =
        matrix(1, 2, 3)
    ),
    "The sum of row and column marginals do not match."
  )
})

test_that("generate_contin_table_with_clustered_AE_with_tol
          errors if AE_idx is missing with rho", {
  row_marginal <- c(10, 20, 30)
  column_marginal <- c(15, 25, 20)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)
  rho <- 0.5

  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(row_marginal,
      column_marginal,
      signal_mat,
      rho = rho
    ),
    "User provided `rho` but the `AE_idx` is not provided."
  )
})

test_that("generate_contin_table_with_clustered_AE_with_tol
          errors if rho is out of range", {
  row_marginal <- c(10, 20, 30)
  column_marginal <- c(15, 25, 20)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)
  rho <- -1.5

  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(row_marginal,
      column_marginal,
      signal_mat,
      rho = rho
    ),
    regexp = "The value of `rho` must lie between \\[0,1\\]."
  )
})

test_that("generate_contin_table_with_clustered_AE_with_tol generates
          tables with RTD within tolerance", {
  row_marginal <- c(10, 20, 30)
  column_marginal <- c(15, 25, 20)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)

  AE_idx <- data.frame(idx = c(1, 1, 2), AE = c("AE1", "AE2", "AE3"))
  rho <- 0.5
  n_rep <- 5
  contin_table <- matrix(c(2, 3, 5, 10, 7, 12, 18, 9, 13), nrow = 3)

  tol <- 0.1
  result <- generate_contin_table_with_clustered_AE_with_tol(
    row_marginal = row_marginal,
    column_marginal = column_marginal,
    signal_mat = signal_mat,
    AE_idx = AE_idx,
    n_rep = n_rep,
    rho = rho,
    tol = tol
  )

  simulated_table_sums <- vapply(result, sum, numeric(1))
  original_sum <- sum(row_marginal)

  rtd_values <- (abs(simulated_table_sums - original_sum) / original_sum) * 100
  expect_true(all(rtd_values <= tol))
})

test_that("generate_contin_table_with_clustered_AE_with_tol calculates rho", {
  row_marginal <- c(10, 20, 30)
  column_marginal <- c(15, 25, 20)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)

  AE_idx <- data.frame(idx = c(1, 1, 2), AE = c("AE1", "AE2", "AE3"))
  n_rep <- 5
  contin_table <- matrix(c(2, 3, 5, 10, 7, 12, 18, 9, 13), nrow = 3)

  tol <- 0.1
  result <- generate_contin_table_with_clustered_AE_with_tol(
    contin_table = contin_table,
    signal_mat = signal_mat,
    n_rep = n_rep,
    tol = tol
  )

  simulated_table_sums <- vapply(result, sum, numeric(1))
  original_sum <- sum(contin_table)

  rtd_values <- (abs(simulated_table_sums - original_sum) / original_sum) * 100
  expect_true(all(rtd_values <= tol))
})

test_that("generate_contin_table_with_clustered_AE_with_tol
          accepts custom covariance matrix", {
  row_marginal <- c(10, 20, 30)
  column_marginal <- c(15, 25, 20)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)
  cov_matrix <- diag(1, 3)

  result <- generate_contin_table_with_clustered_AE_with_tol(
    row_marginal = row_marginal,
    column_marginal = column_marginal,
    AE_idx = NULL,
    signal_mat = signal_mat,
    n_rep = 3,
    rho = cov_matrix
  )

  expect_length(result, 3)
  expect_true(all(vapply(result, is.matrix, logical(1))))
})

test_that("generate_contin_table_with_clustered_AE_with_tol errors
          if AE_idx length does not match table", {
  row_marginal <- c(10, 20, 30)
  column_marginal <- c(15, 25, 20)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)

  AE_idx <- data.frame(idx = c(1, 1), AE = c("AE1", "AE2"))
  rho <- 0.5
  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(
      row_marginal = row_marginal,
      column_marginal = column_marginal,
      signal_mat = signal_mat,
      AE_idx = AE_idx,
      rho = rho
    ),
    "The length of `AE_idx` should be same
              as length of `row_marginal`."
  )
})


test_that("generate_contin_table_with_clustered_AE_with_tol errors
          if corr matrix does not match", {
  row_marginal <- c(10, 20, 30)
  column_marginal <- c(15, 25, 20)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)


  AE_idx <- data.frame(idx = c(1, 1), AE = c("AE1", "AE2"))
  rho <- matrix(0.5, 2, 2)
  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(
      row_marginal = row_marginal,
      column_marginal = column_marginal,
      signal_mat = signal_mat,
      AE_idx = AE_idx,
      rho = rho
    ),
    "Please check the shape of the input matrix `rho`.
              It should be an I x I matrix where I is same length
              as `row_marginal`."
  )
})

test_that("generate_contin_table_with_clustered_AE_with_tol errors if
          corr matrix does not match with contin_table", {
  contin_table <- matrix(20, 3, 3)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)
  rho <- matrix(0.5, 2, 2)
  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(
      row_marginal = NULL,
      column_marginal = NULL,
      contin_table = contin_table,
      signal_mat = signal_mat,
      rho = rho
    ),
    "Please check the shape of the input matrix `rho`.
           It should be an I x I matrix where I is the number of rows
           in the contingency table."
  )
})

test_that("generate_contin_table_with_clustered_AE_with_tol errors
          if rho is not provided", {
  contin_table <- matrix(20, 3, 3)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)
  AE_idx <- c(1, 2, 3)
  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(
      row_marginal = NULL,
      column_marginal = NULL,
      contin_table = contin_table,
      signal_mat = signal_mat,
      AE_idx = AE_idx
    ),
    "User provided the `AE` but `rho` has not been provided.
           If user is unable to provide `rho`, then please set `AE_idx`= NULL."
  )
})

test_that("generate_contin_table_with_clustered_AE_with_tol errors
          if contin_table is not provided", {
  signal_mat <- matrix(1, nrow = 3, ncol = 3)
  AE_idx <- c(1, 2, 3)
  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(
      row_marginal = c(20, 20, 20),
      column_marginal = c(20, 20, 20),
      contin_table = NULL,
      signal_mat = signal_mat
    ),
    "`rho` cannot be estimated if no `contin_table` is provided."
  )
})


test_that("generate_contin_table_with_clustered_AE_with_tol errors if
          AE_idx length does not match contin table", {
  contin_table <- matrix(20, 3, 3)
  signal_mat <- matrix(1, nrow = 3, ncol = 3)

  AE_idx <- data.frame(idx = c(1, 1), AE = c("AE1", "AE2"))
  rho <- 0.5
  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(
      row_marginal = NULL,
      column_marginal = NULL,
      contin_table = contin_table,
      signal_mat = signal_mat,
      AE_idx = AE_idx,
      rho = rho
    ),
    "The length of `AE_idx` should be same as
                rows of `contin_table`."
  )
})

test_that("generate_contin_table_with_clustered_AE_with_tol
          errors with all NULL", {
  row_marginal <- NULL
  column_marginal <- NULL
  signal_mat <- matrix(1, nrow = 3, ncol = 3)

  AE_idx <- data.frame(idx = c(1, 1, 1), AE = c("AE1", "AE2", "AE3"))
  rho <- 0.5
  expect_error(
    generate_contin_table_with_clustered_AE_with_tol(
      row_marginal = row_marginal,
      column_marginal = column_marginal,
      signal_mat = signal_mat,
      AE_idx = AE_idx,
      rho = rho
    ),
    "`row_marginal` or `column_marginal` cannot be
            NULL when `contin_table` is also NULL.
      Please provide either `row_marginal` and `column_marginal`
            or `contin_table`."
  )
})



test_that("generate_contin_table_with_clustered_AE_with_tol
          generates correct output", {
  contin_table <- matrix(c(
    10, 20, 30,
    5, 10, 15,
    2, 4, 6
  ), nrow = 3, byrow = TRUE)
  rownames(contin_table) <- c("AE1", "AE2", "AE3")
  colnames(contin_table) <- c("Drug1", "Drug2", "Drug3")

  AE_idx <- data.frame(
    idx = c("Cluster1", "Cluster2", "Cluster3"),
    AE = c("AE1", "AE2", "AE3")
  )

  lambda_matrix <- matrix(1,
    nrow = nrow(contin_table),
    ncol = ncol(contin_table)
  )
  lambda_matrix[2, 1] <- 4 # Assign signal strength

  # Call the function
  simulated_tables1 <- generate_contin_table_with_clustered_AE_with_tol(
    row_marginal = NULL,
    column_marginal = NULL,
    contin_table = contin_table,
    n_rep = 5,
    AE_idx = AE_idx,
    rho = 0.5,
    signal_mat = lambda_matrix,
    tol = 1
  )

  simulated_tables2 <- generate_contin_table_with_clustered_AE_with_tol(
    row_marginal = NULL,
    column_marginal = NULL,
    contin_table = contin_table,
    n_rep = 5,
    AE_idx = NULL,
    rho = diag(ncol(contin_table)),
    signal_mat = lambda_matrix,
    tol = 1
  )

  simulated_tables3 <- generate_contin_table_with_clustered_AE_with_tol(
    row_marginal = rowSums(contin_table),
    column_marginal = colSums(contin_table),
    n_rep = 5,
    AE_idx = AE_idx,
    rho = diag(ncol(contin_table)),
    signal_mat = lambda_matrix,
    tol = 1
  )

  # Check the output
  expect_equal(length(simulated_tables1), 5,
    info = "Should generate 5 simulated tables."
  )

  # Check dimensions of the simulated tables
  for (sim_table in simulated_tables1) {
    expect_equal(dim(sim_table), dim(contin_table),
      info = "Simulated tables should have the same
      dimensions as the input table."
    )
    expect_true(all(rownames(sim_table) == rownames(contin_table)),
      info = "Row names should match the input table."
    )
    expect_true(all(colnames(sim_table) == colnames(contin_table)),
      info = "Column names should match the input table."
    )
  }

  # Check for non-negative values
  for (sim_table in simulated_tables1) {
    expect_true(all(sim_table >= 0),
      info = "Simulated tables should contain
      only non-negative values."
    )
  }

  # Check the output
  expect_equal(length(simulated_tables2), 5,
    info = "Should generate 5 simulated tables."
  )

  # Check dimensions of the simulated tables
  for (sim_table in simulated_tables2) {
    expect_equal(dim(sim_table), dim(contin_table),
      info = "Simulated tables should have the same
      dimensions as the input table."
    )
    expect_true(all(rownames(sim_table) == rownames(contin_table)),
      info = "Row names should match the input table."
    )
    expect_true(all(colnames(sim_table) == colnames(contin_table)),
      info = "Column names should match the input table."
    )
  }

  # Check for non-negative values
  for (sim_table in simulated_tables2) {
    expect_true(all(sim_table >= 0),
      info = "Simulated tables should contain only non-negative values."
    )
  }

  # Check the output
  expect_equal(length(simulated_tables3), 5,
    info = "Should generate 5 simulated tables."
  )

  # Check dimensions of the simulated tables
  for (sim_table in simulated_tables3) {
    expect_equal(dim(sim_table), dim(contin_table),
      info = "Simulated tables should have the same
      dimensions as the input table."
    )
    expect_true(all(rownames(sim_table) == rownames(contin_table)),
      info = "Row names should match the input table."
    )
    expect_true(all(colnames(sim_table) == colnames(contin_table)),
      info = "Column names should match the input table."
    )
  }

  # Check for non-negative values
  for (sim_table in simulated_tables3) {
    expect_true(all(sim_table >= 0),
      info = "Simulated tables should contain only non-negative values."
    )
  }
})
