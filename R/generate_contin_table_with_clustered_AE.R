#' Generate simulated contingency tables with option of incorporating
#' adverse event correlation within clusters.
#'
#' @description Generate simulated contingency tables with option of
#' incorporating adverse event correlation within clusters.
#'
#' @param contin_table A data matrix of an \eqn{I} x \eqn{J} contingency
#' table with row (adverse event) and column (drug) names, of which the
#' row and column marginals are used to generate the simulated data.
#' Please first check the input contingency table using the function
#' \code{check_and_fix_contin_table()}.
#' @param n_rep Number of contingency tables to be generated.
#' @param AE_idx A data frame with two variables \code{idx} and \code{AE},
#' where \code{idx} indicates the cluster index (can be either a name or
#' a number), and \code{AE} lists the adverse event names. See the
#' \code{statin49_AE_idx} for \code{statin49} data as an example.
#' @param rho A numeric value indicating the correlation of the AEs within
#' each cluster. Default is 0.5.
#' @param signal_mat A data matrix of the same dimension as the contingency
#' table with entries representing the signal strength. The values should
#' be greater or equal to 1, where 1 indicates no signal, and values
#' greater than 1 indicate signal.
#'
#' @return A list of \code{n_rep} simulated contingency tables.
#'
#' @useDynLib MDDC
#' @importFrom stats rnorm
#' @importFrom MASS mvrnorm
#' @export
#'
#' @examples
#' # using statin49 as an example
#' data(statin49)
#' data(statin49_AE_idx)
#'
#' # Prepare a matrix of signal strength with the same dimension as
#' # statin49, where 1 indicates no signal and values > 1 indicate
#' # signal
#' lambda_matrix <- matrix(1, nrow = nrow(statin49), ncol = ncol(statin49))
#'
#' # Assign the cell (45,1) with signal strength 4
#' lambda_matrix[45, 1] <- 4
#'
#' # Generate 5 simulated tables
#' set.seed(123)
#' simulated_tables <- generate_contin_table_with_clustered_AE(
#'   contin_table = statin49,
#'   n_rep = 5,
#'   AE_idx = statin49_AE_idx,
#'   rho = 0.5,
#'   signal_mat = lambda_matrix
#' )
generate_contin_table_with_clustered_AE <- function(contin_table,
                                                    n_rep = 1,
                                                    AE_idx,
                                                    rho = 0.5,
                                                    signal_mat) {
  set.seed(42)

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

  group_s <- unique(AE_idx$idx)
  n_group_s <- unlist(lapply(
    seq_len(length(group_s)),
    function(a) sum(AE_idx$idx == group_s[a])
  ))

  get_new_contin_table <- function(a) {
    Z_ij_mat <- matrix(NA, nrow = n_row, ncol = n_col)

    for (i in seq_len(length(group_s))) {
      if (n_group_s[i] > 1) {
        Z_ij_mat[which(row.names(contin_table) %in%
          AE_idx$AE[which(AE_idx$idx == group_s[i])]), ] <-
          t(MASS::mvrnorm(
            n = n_col,
            mu = rep(0, n_group_s[i]),
            Sigma = correlation_matrix(n_group_s[i], rho)
          ))
      } else {
        Z_ij_mat[which(row_names %in% AE_idx$AE[which(AE_idx$idx ==
          group_s[i])]), ] <-
          rnorm(n_col)
      }
    }

    new_contin_table <- round(Z_ij_mat * sqrt(E_ij_mat * signal_mat *
      ((1 - p_i_dot) %*%
        t((1 - p_dot_j)))) +
      E_ij_mat * signal_mat)
    new_contin_table[which(new_contin_table < 0)] <- 0

    rownames(new_contin_table) <- row_names
    colnames(new_contin_table) <- col_names

    return(new_contin_table)
  }

  return(lapply(seq_len(n_rep), get_new_contin_table))
}
