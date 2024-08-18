#' Verifying and correcting the input \eqn{I} by \eqn{J} contingency table
#' @description Verify and correct the input \eqn{I} by \eqn{J} contingency
#' table to ensure it is properly formatted as a data matrix with row
#' (adverse event) and column (drug) names.
#'
#' @param contin_table A data matrix or a data frame of an \eqn{I} by \eqn{J}
#' contingency table with or without prespecified row and column names.
#'
#' @return A data matrix of an \eqn{I} by \eqn{J} contingency table with row
#' and column names.
#' @export
#'
#' @examples
#' # Create a 6 by 4 data matrix
#' set.seed(42)
#' dat_mat <- matrix(rpois(6 * 4, 20), nrow = 6)
#' dat_mat
#'
#' # Check the format of the data matrix and assign row and column names
#' contin_table <- check_and_fix_contin_table(dat_mat)
#' contin_table
check_and_fix_contin_table <- function(contin_table) {
  # Convert to data matrix if it is a data frame
  if (is.data.frame(contin_table)) {
    contin_table <- as.matrix(contin_table)
  }

  # Check if contin_table is a data matrix with non-negative integers
  if (!is.matrix(contin_table) ||
    !all(contin_table == as.integer(contin_table)) ||
    any(contin_table < 0)) {
    stop(
      "Input is not a data matrix with non-negative integers. Please ",
      "input a correct data matrix."
    )
  }

  # Assign missing rownames or colnames if necessary
  if (is.null(rownames(contin_table))) {
    rownames(contin_table) <- paste0("AE_", seq_len(nrow(contin_table)))
  }
  if (is.null(colnames(contin_table))) {
    colnames(contin_table) <- paste0("drug_", seq_len(ncol(contin_table)))
  }

  return(contin_table)
}
