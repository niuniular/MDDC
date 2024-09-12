#' Generate simulated contingency tables with the option of incorporating
#' adverse event correlation within clusters.
#'
#' @description Generate simulated contingency tables with the option of
#' incorporating adverse event correlation within clusters.
#'
#' @param contin_table A data matrix of an \eqn{I} x \eqn{J} contingency
#' table with row (adverse event) and column (drug) names, of which the
#' row and column marginals are used to generate the simulated data.
#' Please first check the input contingency table using the function
#' \code{check_and_fix_contin_table()}.
#' @param signal_mat A data matrix of the same dimension as the contingency
#' table with entries representing the signal strength. The values should
#' be greater or equal to 1, where 1 indicates no signal, and values
#' greater than 1 indicate signal.
#' @param AE_idx A data frame with two variables \code{idx} and \code{AE},
#' where \code{idx} indicates the cluster index (can be either a name or
#' a number), and \code{AE} lists the adverse event names. See the
#' \code{statin49_AE_idx} for \code{statin49} data as an example.
#' @param n_rep Number of contingency tables to be generated.
#' @param rho A numeric value, matrix, or NULL indicating the correlation
#' structure.If a numeric value (float or int) is provided, it represents the
#' correlation value \code{rho} to be used between all elements within each
#' cluster specified by \code{AE_idx}. If a matrix is provided, it must be a
#' square matrix with dimensions equal to the number of rows in
#' \code{contin_table}. In this case, \code{rho} defines the correlation
#' structure directly, and \code{AE_idx} is not used. If \code{rho} is NULL,
#' a covariance matrix is generated based on the correlation coefficients of
#' \code{contin_table}.
#' @param seed An optional integer to set the seed for reproducibility.
#' If NULL, no seed is set.
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
                                                    signal_mat,
                                                    AE_idx = NULL,
                                                    n_rep = 1,
                                                    rho = NULL,
                                                    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

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

  if (is.numeric(rho) && length(rho) == 1) {
    if (!(0 <= rho && rho <= 1)) {
      stop("The value of `rho` must lie between [0,1]")
    } else {
      if (is.null(AE_idx)) {
        stop("User provided `rho` but the `AE_idx` is not provided.")
      }

      group_s <- unique(AE_idx$idx)
      n_group_s <- unlist(lapply(
        seq_len(length(group_s)),
        function(a) sum(AE_idx$idx == group_s[a])
      ))
      return(lapply(
        seq_len(n_rep),
        get_contin_table_numeric_rho,
        contin_table,
        AE_idx,
        rho,
        group_s,
        n_group_s,
        E_ij_mat,
        signal_mat,
        p_i_dot,
        p_dot_j
      ))
    }
  } else if (is.matrix(rho)) {
    if (!all(dim(rho) == c(dim(contin_table)[1], dim(contin_table)[1]))) {
      stop("Please check the shape of the input matrix `rho`.
           It should be an I x I matrix where I is the number of rows
           in the contingency table.")
    } else {
      cov_matrix <- rho
      return(lapply(
        seq_len(n_rep),
        get_contin_table_matrix_rho,
        contin_table,
        cov_matrix,
        E_ij_mat,
        signal_mat,
        p_i_dot,
        p_dot_j
      ))
    }
  } else if (is.null(rho)) {
    if (!is.null(AE_idx)) {
      stop("User provided the `AE` but `rho` has not been provided.
           If user is unable to provide `rho`, then please set `AE_idx`= NULL.")
    }
    cov_matrix <- cor(t(contin_table))
    return(lapply(
      seq_len(n_rep),
      get_contin_table_matrix_rho,
      contin_table,
      cov_matrix,
      E_ij_mat,
      signal_mat,
      p_i_dot,
      p_dot_j
    ))
  }
}
