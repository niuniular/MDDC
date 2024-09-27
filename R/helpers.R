#' Generate Log Bootstrap Cutoffs for Contingency Table
#'
#' This function computes log-transformed bootstrap cutoffs for each column
#' of a contingency table based on a specified quantile. The cutoffs are
#' derived from a simulated null distribution using multinomial resampling
#' of the contingency table.
#'
#' @param contin_table A matrix representing the contingency table.
#' @param quantile A numeric value between 0 and 1 representing the quantile
#' to compute for the bootstrap cutoff.
#' @param rep An integer specifying the number of bootstrap replications to
#' generate.
#' @param seed An optional integer seed for reproducibility of the random
#' number generator.
#'
#' @return A list containing:
#' \describe{
#'   \item{cutoffs}{A numeric vector with the bootstrap cutoff values for
#'   each column.}
#'   \item{null_dist_s}{A matrix representing the null distributions for
#'   each column across all bootstrap replications.}
#' }
#'
#' @useDynLib MDDC
#' @importFrom stats quantile
#' @importFrom stats rmultinom
#' @noRd
get_log_bootstrap_cutoff <- function(contin_table, quantile, rep, seed) {
  n_col <- ncol(contin_table)
  n_i_dot <- rowSums(contin_table)
  n_dot_j <- colSums(contin_table)
  n_dot_dot <- sum(contin_table)

  p_i_dot <- n_i_dot / n_dot_dot
  p_dot_j <- n_dot_j / n_dot_dot

  p_mat <- p_i_dot %*% t(p_dot_j)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  sim_tab_list <- lapply(seq_len(rep), function(a) {
    matrix(rmultinom(1, n_dot_dot, p_mat), ncol = n_col)
  })

  Z_ij_mat_list <- lapply(sim_tab_list, getZijMat)

  max_list <- do.call(rbind, lapply(Z_ij_mat_list, function(a) {
    apply(log(a), 2, max, na.rm = TRUE)
  }))

  cutoffs <- unlist(lapply(seq_len(n_col), function(a) {
    quantile(max_list[, a], quantile, na.rm = TRUE)
  }))

  names(cutoffs) <- colnames(contin_table)

  rslt <- list(cutoffs, max_list)
  names(rslt) <- c("cutoffs", "null_dist_s")
  return(rslt)
}


#' Calculate Fisher's Exact Test p-value for a Contingency Table
#'
#' This function computes the p-value from Fisher's Exact Test for a 2x2
#' contingency table derived from specified row and column indices. The
#' table is populated based on the values in the given contingency table,
#' and the p-value is calculated using Fisher's Exact Test.
#'
#' @param row_idx An integer specifying the row index to be considered
#' in the contingency table.
#' @param col_idx An integer specifying the column index to be considered
#' in the contingency table.
#'
#' @return A numeric value representing the p-value from Fisher's Exact Test.
#'
#' @details The function constructs a 2x2 contingency table based on the
#' provided indices and computes the p-value using Fisher's Exact Test.
#' The table is filled as follows:
#' \itemize{
#'   \item Row 1, Column 1: Value from the specified cell in the
#'      contingency table.
#'   \item Row 2, Column 1: Sum of the column excluding the specified row.
#'   \item Row 1, Column 2: If \code{exclude_same_drug_class} is \code{TRUE},
#'   the value from the column corresponding to \code{n_col}. Otherwise,
#'   the sum of the row excluding the specified column.
#'   \item Row 2, Column 2: If \code{exclude_same_drug_class} is \code{TRUE},
#'   the sum of the column corresponding to \code{n_col}. Otherwise, the sum
#'   of the matrix excluding the specified row and column.
#' }
#'
#' @seealso \code{\link{fisher.test}} for more details on Fisher's Exact Test.
#'
#' @useDynLib MDDC
#' @importFrom stats fisher.test
#' @noRd
get_fisher <- function(
    contin_table,
    row_idx,
    col_idx,
    exclude_same_drug_class) {
  tabl <- getFisherExactTestTable(
    contin_table,
    row_idx,
    col_idx,
    exclude_same_drug_class
  )
  return(fisher.test(tabl)$p.value)
}

#' Helper function to identify outliers using Boxplot
#'
#' This function identifies the outliers in a dataset by using
#' \eqn{Q_3 + c_j \times IQR}.
#'
#' @param dat A numeric vector containing the data for which outliers
#' are to be identified.
#' @param c_j A numeric value used as a scaling factor to determine
#' outliers.
#' @return A logical vector indicating TRUE for values that are
#' outliers and FALSE otherwise.
#'
#' @useDynLib MDDC
#' @importFrom stats IQR
#' @noRd
get_boxplot_outliers <- function(dat, c_j) {
  outliers <-
    dat > (quantile(dat, 0.75, na.rm = TRUE) +
      c_j * IQR(dat, na.rm = TRUE))
  return(outliers)
}

#' Helper function to compute the FDR
#'
#' This function calculates the average number
#' of outliers for a given column across a number of
#' datasets.
#'
#' @param res_list A list of numeric matrices or data frames. Each element
#' of the list represents a dataset where outliers will be identified.
#' @param c_j A numeric value used as a scaling factor to determine
#' outliers.
#' @param j An integer representing the index of the column to
#' compute the FDR on.
#' @return A numeric value representing the average number of
#' outliers across all datasets for the given column.
#'
#' @useDynLib MDDC
#' @noRd
compute_fdr <- function(res_list, c_j, j) {
  mean(unlist(lapply(res_list, function(a) {
    sum(get_boxplot_outliers(a[, j], c_j), na.rm = TRUE)
  })))
}


#' Helper function to compute the FDR
#'
#' This function calculates the average number
#' of outliers across a number of
#' datasets.
#'
#' @param res_list A list of numeric matrices or data frames. Each element
#' of the list represents a dataset where outliers will be identified.
#' @param c A numeric value used as a scaling factor to determine
#' outliers.
#' @return A numeric value representing the average number of
#' outliers across all datasets.
#'
#' @useDynLib MDDC
#' @noRd
compute_fdr_all <- function(res_list, c) {
  mean(unlist(lapply(res_list, function(a) {
    sum(get_boxplot_outliers(a, c), na.rm = TRUE)
  })))
}


#' Create a Block Diagonal Matrix
#'
#' This function constructs a block diagonal matrix from a list of matrices.
#' Each input matrix is placed on the diagonal of the resulting matrix, with
#' all off-diagonal blocks filled with zeros.
#'
#' @param ... A series of matrices to be placed on the diagonal of the resulting
#' block diagonal matrix. Each matrix should have the same number of columns in
#' order to ensure proper placement.
#'
#' @return A matrix where the input matrices are arranged along the diagonal
#' and all off-diagonal entries are zero.
#'
#' @useDynLib MDDC
#' @noRd
create_block_diagonal_matrix <- function(matrices) {
  total_rows <- sum(vapply(matrices, nrow, integer(1)))
  total_cols <- sum(vapply(matrices, ncol, integer(1)))
  result <- matrix(0, nrow = total_rows, ncol = total_cols)
  current_row <- 1
  current_col <- 1

  for (matrix in matrices) {
    r <- nrow(matrix)
    c <- ncol(matrix)
    result[
      current_row:(current_row + r - 1),
      current_col:(current_col + c - 1)
    ] <- matrix
    current_row <- current_row + r
    current_col <- current_col + c
  }

  return(result)
}


#' Generate a Simulated Contingency Table with specified correlation
#' within cluster
#'
#' This function generates a new contingency table by introducing correlation
#' between AEs within cluster.
#'
#' @param contin_table A data matrix of an \eqn{I} x \eqn{J} contingency
#' table with row (adverse event) and column (drug) names, of which the
#' row and column marginals are used to generate the simulated data.
#' Please first check the input contingency table using the function
#' @param AE_idx A data frame with two variables \code{idx} and \code{AE},
#' where \code{idx} indicates the cluster index (can be either a name or
#' a number), and \code{AE} lists the adverse event names. See the
#' \code{statin49_AE_idx} for \code{statin49} data as an example.
#' @param rho A numeric value indicating the correlation of the
#' AEs within each cluster.
#' @param group_s Unique clusters in the provided \code{AE_idx}.
#' @param n_group_s Count of number of AEs within each cluster.
#' @param E_ij_mat Matrix of expected counts.
#' @param signal_mat Signal matrix.
#' @param p_i_dot Row marginal probabilities.
#' @param p_dot_j Column marginal probabilities.
#'
#' @return Simulated contingency table.
#'
#' @useDynLib MDDC
#' @importFrom MASS mvrnorm
#' @noRd
get_contin_table <- function(a,
                             n_row,
                             n_col,
                             cov_matrix,
                             E_ij_mat,
                             signal_mat,
                             p_i_dot,
                             p_dot_j) {
  Z_ij_mat <- t(mvrnorm(
    n = n_col,
    mu = rep(0, n_row),
    Sigma = cov_matrix
  ))
  new_contin_table <- round(Z_ij_mat * sqrt(E_ij_mat * signal_mat *
    ((1 - p_i_dot) %*%
      t((1 - p_dot_j)))) +
    E_ij_mat * signal_mat)
  new_contin_table[which(new_contin_table < 0)] <- 0
  return(new_contin_table)
}
