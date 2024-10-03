#' Generate simulated contingency tables with the option of incorporating
#' adverse event correlation within clusters and tolerance for total report
#' count.
#'
#' @description Generate simulated contingency tables with the option of
#' incorporating adverse event correlation within clusters and tolerance
#' for total report count.
#'
#' @param row_marginal Marginal sums for the rows of the contingency table.
#' @param column_marginal Marginal sums for the columns of the contingency
#' table.
#' @param signal_mat A data matrix of the same dimension as the contingency
#' table with entries representing the signal strength. The values should
#' be greater or equal to 1, where 1 indicates no signal, and values
#' greater than 1 indicate signal.
#' @param tol Tolerance for the total report count, expressed in
#' terms of the Relative Total Difference (RTD), defined as:
#'
#' \deqn{RTD = \frac{|n^{orig}_{\cdot \cdot} -
#' n^{sim}_{\cdot \cdot}|}{n^{orig}_{\cdot \cdot}} \times 100}
#'
#' This represents the difference between the total number of reports in the
#' simulated datasets and the original input total number of reports.
#' Sufficiently low tolerance will generate tables with total report counts
#' equal to the actual supplied value. Default is 0.1.
#' @param contin_table A data matrix of an \eqn{I} x \eqn{J} contingency
#' table with row (adverse event) and column (drug or vaccine) names,
#' of which the row and column marginals are used to generate the simulated
#' data. Please first check the input contingency table using the function
#' \code{check_and_fix_contin_table()}. Default is NULL.
#' @param AE_idx A data frame or list.
#' In case of data frame it must contain two variables \code{idx} and \code{AE},
#' where \code{idx} indicates the cluster index (a number),
#' and \code{AE} lists the adverse event names. See the
#' \code{statin49_AE_idx} for \code{statin49} data as an example.
#' In case of a list, make sure the cluster indices are aligned with the
#' corresponding row marginal.
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
#' lambda_matrix[1, 1] <- 4
#'
#' # Generate 5 simulated tables
#' simulated_tables <- generate_contin_table_with_clustered_AE_with_tol(
#'   contin_table = statin49,
#'   signal_mat = lambda_matrix,
#'   n_rep = 5,
#'   AE_idx = statin49_AE_idx,
#'   rho = 0.5,
#'   tol = 0.1
#' )
generate_contin_table_with_clustered_AE_with_tol <-
  function(row_marginal,
           column_marginal,
           signal_mat,
           tol = 0.1,
           contin_table = NULL,
           AE_idx = NULL,
           n_rep = 1,
           rho = NULL) {
    if (is.data.frame(AE_idx)) {
      AE_idx <- as.list(AE_idx$idx)
    }

    if (!is.null(contin_table)) {
      n_row <- nrow(contin_table)
      n_col <- ncol(contin_table)
      row_names <- row.names(contin_table)
      col_names <- colnames(contin_table)

      n_i_dot <- rowSums(contin_table)
      n_dot_j <- colSums(contin_table)
      n_dot_dot <- sum(contin_table)
    } else if (is.null(contin_table) &&
      (!is.null(row_marginal) &&
        !is.null(column_marginal))) {
      if (sum(row_marginal) == sum(column_marginal)) {} else {
        stop("The sum of row and column marginals do not match.")
      }

      n_row <- length(row_marginal)
      n_col <- length(column_marginal)

      n_i_dot <- row_marginal
      n_dot_j <- column_marginal
      n_dot_dot <- sum(row_marginal)
    } else {
      if ((is.null(row_marginal) || is.null(column_marginal)) &&
        is.null(contin_table)) {
        stop("`row_marginal` or `column_marginal` cannot be
            NULL when `contin_table` is also NULL.
      Please provide either `row_marginal` and `column_marginal`
            or `contin_table`.")
      }
    }

    p_i_dot <- n_i_dot / n_dot_dot
    p_dot_j <- n_dot_j / n_dot_dot

    E_ij_mat <- n_i_dot %*% t(n_dot_j) / n_dot_dot

    if (is.numeric(rho) && length(rho) == 1) {
      if (!(0 <= rho && rho <= 1)) {
        stop("The value of `rho` must lie between [0,1].")
      } else {
        if (!is.null(contin_table)) {
          if (length(AE_idx) != dim(contin_table)[1]) {
            stop("The length of `AE_idx` should be same as
                rows of `contin_table`.")
          }
        } else {
          if (is.null(AE_idx)) {
            stop("User provided `rho` but the `AE_idx` is not provided.")
          } else {
            if (length(AE_idx) != length(row_marginal)) {
              stop("The length of `AE_idx` should be same
              as length of `row_marginal`.")
            }
          }
        }
        group_s <- unique(AE_idx)
        n_group_s <- vapply(group_s, function(g) sum(AE_idx == g), numeric(1))
        names(n_group_s) <- group_s

        cov_matrices <- list()

        for (group in group_s) {
          size <- n_group_s[group]
          matrix <- matrix(rho, nrow = size, ncol = size)
          diag(matrix) <- 1
          cov_matrices[[group]] <- matrix
        }
        cov_matrix <- create_block_diagonal_matrix(cov_matrices)
      }
    } else if (is.matrix(rho)) {
      if (!is.null(contin_table)) {
        if (!all(dim(rho) == c(dim(contin_table)[1], dim(contin_table)[1]))) {
          stop("Please check the shape of the input matrix `rho`.
           It should be an I x I matrix where I is the number of rows
           in the contingency table.")
        }
      } else {
        if (!all(dim(rho) == c(length(row_marginal), length(row_marginal)))) {
          stop("Please check the shape of the input matrix `rho`.
              It should be an I x I matrix where I is same length
              as `row_marginal`.")
        }
      }

      cov_matrix <- rho
    } else if (is.null(rho)) {
      if (!is.null(AE_idx)) {
        stop("User provided the `AE` but `rho` has not been provided.
           If user is unable to provide `rho`, then please set `AE_idx`= NULL.")
      }
      if (!is.null(contin_table)) {
        cov_matrix <- cor(t(contin_table))
      } else {
        stop("`rho` cannot be estimated if no `contin_table` is provided.")
      }
    }

    tables <- lapply(seq_len(n_rep), function(i) {
      get_contin_table(
        i, n_row, n_col, cov_matrix, E_ij_mat,
        signal_mat, p_i_dot, p_dot_j
      )
    })

    simulated_table_sums <- vapply(tables, sum, numeric(1))
    rtd_values <- (abs(simulated_table_sums - sum(n_i_dot))
    / sum(n_i_dot)) * 100

    if (max(rtd_values) > tol) {
      indices <- which(rtd_values > tol)
      num_rejected_samples <- length(indices)
      new_samples <- list()

      for (i in seq_len(num_rejected_samples)) {
        current_rtd <- Inf
        while (current_rtd > tol) {
          sample_replacement <- get_contin_table(
            1, n_row, n_col, cov_matrix, E_ij_mat,
            signal_mat, p_i_dot, p_dot_j
          )

          current_rtd <- (abs(sum(n_i_dot) - sum(sample_replacement))
          / sum(n_i_dot)) * 100
        }
        new_samples[[i]] <- sample_replacement
      }
      accepted_samples <- tables[-indices]
      tables <- c(accepted_samples, new_samples)
    }


    if (!is.null(contin_table)) {
      if (all(!is.na(colnames(contin_table))) &&
        all(!is.na(rownames(contin_table)))) {
        tables <- lapply(tables, function(sample) {
          df <- as.matrix(sample)
          colnames(df) <- colnames(contin_table)
          rownames(df) <- rownames(contin_table)
          return(df)
        })
      }
    }

    return(tables)
  }
