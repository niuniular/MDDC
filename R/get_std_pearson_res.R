#' Computing the standardized Pearson residuals for a given \eqn{I} by \eqn{J} contingency table
#'
#' @description Compute the standardized Pearson residuals for a given \eqn{I} by \eqn{J} contingency table.
#'
#' @usage get_std_pearson_res(contin_table)
#'
#' @param contin_table A data matrix of an \eqn{I} by \eqn{J} contingency table
#'
#' @return A data matrix of the standardized Pearson residuals of the input contingency table
#' @export
#'
#' @examples
#' #' # create a 6 by 4 data matrix
#' set.seed(42)
#' dat_mat <- matrix(rpois(6*4, 20), nrow=6)
#' dat_mat
#'
#' # check the format of the data matrix
#' # and assign row and column names to the data matrix
#' contin_table <- check_and_fix_contin_table(dat_mat)
#' contin_table
#'
#' # computed the expected count of the data matrix
#' get_std_pearson_res(contin_table)
get_std_pearson_res <- function(contin_table){

    n_row <- nrow(contin_table)
    n_col <- ncol(contin_table)

    row_names <- row.names(contin_table)
    col_names <- colnames(contin_table)

    n_i_dot <- rowSums(contin_table)
    n_dot_j <- colSums(contin_table)
    n_dot_dot <- sum(contin_table)

    p_i_dot <- n_i_dot / n_dot_dot
    p_dot_j <- n_dot_j / n_dot_dot

    E_ij_mat <- n_i_dot %*% t(n_dot_j) / n_dot_dot
    Z_ij_mat <-
      (contin_table - E_ij_mat) / sqrt(E_ij_mat * ((1 - p_i_dot) %*% t((1 - p_dot_j))))

    row.names(Z_ij_mat) <- row_names
    colnames(Z_ij_mat) <- col_names

    return(Z_ij_mat)

}
