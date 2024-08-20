# Load necessary libraries
library(testthat)
library(MASS) # For mvrnorm function

# Define the test case
test_that("generate_contin_table_with_clustered_AE generates correct output", {
  # Sample input data
  contin_table <- matrix(c(
    10, 20, 30,
    5, 10, 15,
    2, 4, 6
  ), nrow = 3, byrow = TRUE)
  rownames(contin_table) <- c("AE1", "AE2", "AE3")
  colnames(contin_table) <- c("Drug1", "Drug2", "Drug3")

  AE_idx <- data.frame(
    idx = c("Cluster1", "Cluster1", "Cluster2", "Cluster2"),
    AE = c("AE1", "AE2", "AE2", "AE3")
  )

  lambda_matrix <- matrix(1,
    nrow = nrow(contin_table),
    ncol = ncol(contin_table)
  )
  lambda_matrix[2, 1] <- 4 # Assign signal strength

  # Set seed for reproducibility
  set.seed(123)

  # Call the function
  simulated_tables <- generate_contin_table_with_clustered_AE(
    contin_table = contin_table,
    n_rep = 5,
    AE_idx = AE_idx,
    rho = 0.5,
    signal_mat = lambda_matrix
  )

  # Check the output
  expect_equal(length(simulated_tables), 5,
    info = "Should generate 5 simulated tables."
  )

  # Check dimensions of the simulated tables
  for (sim_table in simulated_tables) {
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
  for (sim_table in simulated_tables) {
    expect_true(all(sim_table >= 0),
      info = "Simulated tables should contain
                  only non-negative values."
    )
  }
})
