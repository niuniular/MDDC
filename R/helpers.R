#' Compute Correlation Matrix with Handling of Missing Values
#'
#' This function calculates a correlation matrix from the input matrix,
#' handling missing values (`NA`) by computing correlations only for
#' non-missing pairs. The function supports computing either row-wise or
#' column-wise correlations, based on the input flag `if_col_cor`.
#'
#' @param mat A numeric matrix or data frame. Each column represents a variable,
#'   and each row represents an observation. The function will compute pairwise
#'   correlations between the columns (or rows if `if_col_cor = FALSE`).
#' @param if_col_cor Logical. If `TRUE`, the function computes correlations
#'   between columns of `mat`. If `FALSE`, correlations are computed between
#'   rows (the matrix is transposed internally).
#'
#' @return A square matrix containing the computed correlation coefficients.
#'   The dimensions of the returned matrix will be equal to the number of
#'   columns (or rows if `if_col_cor = FALSE`) in the input matrix `mat`.
#'   Missing values (`NA`) in the input matrix are handled by only using
#'   observations where both variables have non-missing values. If fewer than
#'   3 valid pairs are available for a correlation, `NA` is returned for
#'   that pair.
#'
#' @details The function iterates through each pair of columns (or rows)
#'   and computes the correlation using Pearson's method, only including
#'   observations where both variables have non-missing values. If the number
#'   of valid pairs is less than 3, the function assigns `NA` to the
#'   corresponding entry in the correlation matrix.
#'
#' @useDynLib MDDC
#' @noRd
cor_with_NA <- function(mat, if_col_cor) {
    if (if_col_cor == FALSE) {
        mat <- t(mat)
    }

    n_col <- ncol(mat)
    cor_mat <- matrix(NA, nrow = n_col, ncol = n_col)
    row.names(cor_mat) <- colnames(mat)
    colnames(cor_mat) <- colnames(mat)

    for (i in seq_len(n_col)) {
        for (j in setdiff(seq_len(n_col), i)) {
            idx <- which((!is.na(mat[, i])) & (!is.na(mat[, j])))
            cor_mat[i, j] <- ifelse((length(idx) >= 3), cor(mat[
                idx,
                i
            ], mat[idx, j]), NA)
        }
    }
    return(cor_mat)
}

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

    max_list <- matrix(unlist(lapply(Z_ij_mat_list, function(a) {
        apply(a, 2, function(b) {
            max(log(b), na.rm = TRUE)
        })
    })), byrow = TRUE, ncol = n_col)

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


#' Create a Correlation Matrix
#'
#' This function generates an n x n matrix where the diagonal elements are 1 
#' and all off-diagonal elements are set to a given correlation coefficient \code{rho}.
#'
#' @param n Integer. The number of rows and columns of the matrix.
#' @param rho Numeric. The correlation coefficient for the off-diagonal elements.
#' @return A numeric matrix of size n x n with 1 on the diagonal and \code{rho} on the off-diagonal.
#' @examples
#' correlation_matrix(3, 0.5)
#' correlation_matrix(4, 0.8)
#' @useDynLib MDDC
#' @noRd
correlation_matrix <- function(n, rho) {
  # Create an n x n matrix filled with rho
  mat <- matrix(rho, nrow = n, ncol = n)
  
  # Set the diagonal elements to 1
  diag(mat) <- 1
  
  return(mat)
}
