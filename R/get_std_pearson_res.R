#' Computing the standardized Pearson residuals for a given
#' \eqn{I \times J} contingency table
#'
#' @description Compute the standardized Pearson residuals for a given
#' \eqn{I \times J} contingency table.
#'
#' @usage
#' get_std_pearson_res(contin_table)
#'
#' @param contin_table A matrix of an \eqn{I \times J} contingency table.
#'
#' @return A matrix of the standardized Pearson residuals of the input
#' contingency table.
#' @export
#'
#' @examples
#' # Create a 6 by 4 data matrix
#' set.seed(42)
#' dat_mat <- matrix(rpois(6 * 4, 20), nrow = 6)
#' dat_mat
#'
#' # Check the format of the data matrix and assign row and column names
#' # to the data matrix
#' contin_table <- check_and_fix_contin_table(dat_mat)
#' contin_table
#'
#' # Compute the standardized Pearson residuals
#' get_std_pearson_res(contin_table)
#' @useDynLib MDDC

get_std_pearson_res <- function(contin_table) {
  std_pearson_res <- getZijMat(contin_table, FALSE)
  rownames(std_pearson_res) <- row.names(contin_table)
  colnames(std_pearson_res) <- colnames(contin_table)
  return(std_pearson_res)
}
