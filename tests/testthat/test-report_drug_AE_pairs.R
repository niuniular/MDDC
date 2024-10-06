library(testthat)

# Define test cases for the `report_drug_AE_pairs` function
test_report_drug_AE_pairs <- function() {
  # Sample data
  contin_table1 <- matrix(c(10, 5, 0, 2, 8, 6),
    nrow = 2, ncol = 3,
    dimnames = list(
      c("AE1", "AE2"),
      c("Drug1", "Drug2", "Drug3")
    )
  )

  contin_table2 <- matrix(c(10, 5, 0, 2, 8, 6),
    nrow = 2, ncol = 3
  )

  contin_table_signal1 <- matrix(c(1, 0, 0, 1, 1, 0),
    nrow = 2, ncol = 3,
    dimnames = list(
      c("AE1", "AE2"),
      c("Drug1", "Drug2", "Drug3")
    )
  )

  contin_table_signal2 <- matrix(c(1, 0, 0, 1, 1, 0),
    nrow = 2, ncol = 3
  )

  # Define expected results
  expected_result1 <- data.frame(
    drug = c("Drug1", "Drug2", "Drug3"),
    AE = c("AE1", "AE2", "AE1"),
    observed_num = c("10", "2", "8"),
    expected_num = c("8.7097", "0.8387", "8.129"),
    std_pearson_res = c("0.9398", "1.7205", "-0.0944")
  )

  expected_result2 <- data.frame(
    drug = c("Drug_1", "Drug_2", "Drug_3"),
    AE = c("AE_1", "AE_2", "AE_1"),
    observed_num = c("10", "2", "8"),
    expected_num = c("8.7097", "0.8387", "8.129"),
    std_pearson_res = c("0.9398", "1.7205", "-0.0944")
  )

  # Run tests
  test_that("report_drug_AE_pairs works with valid input1", {
    result1 <- report_drug_AE_pairs(contin_table1, contin_table_signal1)
    expect_equal(nrow(result1), nrow(expected_result1))
    expect_equal(result1, expected_result1, tolerance = 1e-4)
  })

  test_that("report_drug_AE_pairs works with valid input2", {
    result2 <- report_drug_AE_pairs(contin_table2, contin_table_signal2)
    expect_equal(nrow(result2), nrow(expected_result2))
    expect_equal(result2, expected_result2, tolerance = 1e-4)
  })

  test_that("report_drug_AE_pairs checks input types", {
    expect_error(
      report_drug_AE_pairs("invalid", contin_table_signal1),
      "Both inputs must be data matrices."
    )
    expect_error(
      report_drug_AE_pairs(contin_table1, "invalid"),
      "Both inputs must be data matrices."
    )
  })

  test_that("report_drug_AE_pairs checks dimensions", {
    invalid_table <- matrix(c(10, 5, 0, 2),
      nrow = 2, ncol = 2,
      dimnames = list(
        c("AE1", "AE2"),
        c("Drug1", "Drug2")
      )
    )
    expect_error(
      report_drug_AE_pairs(contin_table1, invalid_table),
      "The dimensions of contin_table and contin_table_signal must be the same."
    )
  })

  test_that("report_drug_AE_pairs checks row and column names", {
    invalid_table <- contin_table_signal1
    rownames(invalid_table) <- c("AE1", "AE3") # Different row names
    expect_error(
      report_drug_AE_pairs(contin_table1, invalid_table),
      "The row names of contin_table and contin_table_signal must match."
    )

    invalid_table <- contin_table_signal1
    colnames(invalid_table) <- c("D1", "D3", "D4") # Different row names
    expect_error(
      report_drug_AE_pairs(contin_table1, invalid_table),
      "The column names of contin_table and contin_table_signal must match."
    )
  })
}

# Run the tests
test_report_drug_AE_pairs()
