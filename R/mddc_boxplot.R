#' Modified Detecting Deviating Cells (MDDC) algorithm for adverse event signal
#' identification with boxplot method for cutoff selection.
#'
#' @description Modified Detecting Deviating Cells (MDDC) algorithm for adverse
#' event signal identification. Boxplot method is used for cutoff selection in
#' step 2 of the algorithm.
#' @param contin_table A data matrix of an \eqn{I} x \eqn{J} contingency table
#' with row (adverse event) and column (drug) names. Please first check the
#' input contingency table using the function
#' \code{check_and_fix_contin_table()}.
#' @param col_specific_cutoff Logical. In the step 2 of the algorithm, whether
#' to apply boxplot method to the standardized Pearson residuals of the entire
#' table, or within each drug column. Default is \code{TRUE}, that is within
#' each drug column (column specific cutoff). \code{FALSE} indicates applying
#' boxplot method on residuals of the entire table.
#' @param separate Logical. In the step 2 of the algorithm, whether to separate
#' the standardized Pearson residuals for the zero cells and non zero cells and
#' apply boxplot method separately or together. Default is \code{TRUE}.
#' @param if_col_cor Logical. In the step 3 of the algorithm, whether to use
#' column (drug) correlation or row (adverse event) correlation. Default is
#' \code{FALSE}, that is using the adverse event correlation. \code{TRUE}
#' indicates using drug correlation.
#' @param cor_lim A numeric value between (0, 1). In the step 3,
#' what correlation threshold should be used to select ``connected''
#' adverse events. Default is 0.8.
#' @param coef A numeric value used for determining the cutoff corresponding
#' to the non-zero cells.
#' @param num_cores Number of cores used to parallelize the MDDC
#' Boxplot algorithm. Default is 2.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{boxplot_signal} returns the signals identified in the step 2.
#' 1 indicates signals, 0 for non signal.
#' \item \code{corr_signal_pval} returns the p values for each cell in the
#' contingency table in the step 5, when the \eqn{r_{ij}} values are mapped
#' back to the standard normal distribution.
#' \item \code{corr_signal_adj_pval} returns the Benjamini-Hochberg adjusted
#' p values for each cell in the step 5. We leave here an option for the user
#' to decide whether to use \code{corr_signal_pval} or
#' \code{corr_signal_adj_pval}, and what threshold for p values should be
#' used (for example, 0.05). Please see the example below.
#' }
#' @export
#'
#' @seealso \code{\link{find_optimal_coef}} for finding an optimal value of
#' \code{coef}.
#'
#' @examples
#' # using statin49 data set as an example
#' data(statin49)
#'
#' # apply the mddc_boxplot
#' boxplot_res <- mddc_boxplot(statin49)
#'
#' # signals identified in step 2 using boxplot method
#' signal_step2 <- boxplot_res$boxplot_signal
#'
#' # signals identified in step 5 by considering AE correlations
#' # In this example, cells with p values less than 0.05 are
#' # identified as signals
#' signal_step5 <- (boxplot_res$corr_signal_pval < 0.05) * 1
#'
#' @useDynLib MDDC
#' @importFrom grDevices boxplot.stats
#' @importFrom stats cor
#' @importFrom stats weighted.mean
#' @importFrom stats sd
#' @importFrom stats pnorm
#' @importFrom stats p.adjust
#' @importFrom stats lm
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster

mddc_boxplot <- function(
    contin_table,
    col_specific_cutoff = TRUE,
    separate = TRUE,
    if_col_cor = FALSE,
    cor_lim = 0.8,
    coef = 1.5,
    num_cores = 2) {
  n_row <- nrow(contin_table)
  n_col <- ncol(contin_table)

  row_names <- row.names(contin_table)
  col_names <- colnames(contin_table)

  Z_ij_mat <- getZijMat(continTable = contin_table, na = FALSE)

  res_all <- as.vector(Z_ij_mat)
  res_nonzero <- as.vector(Z_ij_mat[which(contin_table != 0)])
  res_zero <- as.vector(Z_ij_mat[which(contin_table == 0)])

  if (col_specific_cutoff == TRUE) {
    if (separate == TRUE) {
      c_univ_drug <- unlist(lapply(seq_len(n_col), function(a) {
        boxplot.stats(Z_ij_mat[which(contin_table[
          ,
          a
        ] != 0), a], coef = coef)$stats[[5]]
      }))
      zero_drug_cutoff <- unlist(lapply(seq_len(n_col), function(a) {
        boxplot.stats(Z_ij_mat[which(contin_table[
          ,
          a
        ] == 0), a])$stats[[1]]
      }))
    } else {
      c_univ_drug <-
        apply(Z_ij_mat, 2, function(a) boxplot.stats(a, coef = coef)$stats[[5]])
      zero_drug_cutoff <-
        apply(Z_ij_mat, 2, function(a) boxplot.stats(a)$stats[[1]])
    }
  } else {
    if (separate == TRUE) {
      c_univ_drug <- rep(boxplot.stats(res_nonzero)$stats[5], n_col)
      zero_drug_cutoff <- rep(boxplot.stats(res_zero)$stats[1], n_col)
    } else {
      c_univ_drug <- rep(boxplot.stats(res_all)$stats[5], n_col)
      zero_drug_cutoff <- rep(boxplot.stats(res_all)$stats[1], n_col)
    }
  }


  high_outlier <- matrix(unlist(lapply(seq_len(n_col), function(a) {
    (Z_ij_mat[
      ,
      a
    ] > c_univ_drug[a]) * 1
  })), ncol = n_col)
  colnames(high_outlier) <- colnames(contin_table)
  row.names(high_outlier) <- row.names(contin_table)

  low_outlier <- matrix(unlist(lapply(seq_len(n_col), function(a) {
    (Z_ij_mat[
      ,
      a
    ] < -c_univ_drug[a]) * 1
  })), ncol = n_col)
  colnames(low_outlier) <- colnames(contin_table)
  row.names(low_outlier) <- row.names(contin_table)

  zero_cell_out <- matrix(unlist(lapply(seq_len(n_col), function(a) {
    (Z_ij_mat[
      ,
      a
    ] < zero_drug_cutoff[a])
  })), ncol = n_col)
  zero_cell_outlier <- ((zero_cell_out) & (contin_table == 0)) * 1

  if_outlier_mat <- ((high_outlier + low_outlier + zero_cell_outlier) !=
    0) * 1

  U_ij_mat <- ifelse(1 - if_outlier_mat, Z_ij_mat, NA)

  cor_U <- cor_with_NA(U_ij_mat, if_col_cor)

  if (if_col_cor == TRUE) {
    iter_over <- n_col
  } else {
    iter_over <- n_row
  }

  num_cores <- as.numeric(num_cores)
  registerDoParallel(num_cores)

  results <- foreach(i = seq_len(iter_over), .packages = c("stats")) %dopar% {
    # nocov start
    idx <- which(abs(cor_U[i, ]) >= cor_lim)
    cor_list_i <- idx[!idx %in% i]
    weight_list_i <- abs(cor_U[i, cor_list_i])

    if (length(cor_list_i) == 0) {
      fitted_value_i <- NA
    } else {
      if (if_col_cor == TRUE) {
        mat <- matrix(NA, n_row, length(cor_list_i))
        row.names(mat) <- row_names
        colnames(mat) <- col_names[cor_list_i]
      } else {
        mat <- matrix(NA, length(cor_list_i), n_col)
        row.names(mat) <- row_names[cor_list_i]
        colnames(mat) <- col_names
      }

      coef_list <- list()
      for (j in seq_len(length(cor_list_i))) {
        if (if_col_cor == TRUE) {
          coeff <- lm(U_ij_mat[, i] ~ U_ij_mat[, cor_list_i[j]])$coefficients
          fit_values <- U_ij_mat[, cor_list_i[j]] * coeff[2] +
            coeff[1]
          mat[which(row_names %in% names(fit_values)), j] <- fit_values
        } else {
          coeff <- lm(U_ij_mat[i, ] ~ U_ij_mat[cor_list_i[j], ])$coefficients
          fit_values <- U_ij_mat[cor_list_i[j], ] * coeff[2] +
            coeff[1]
          mat[j, which(col_names %in% names(fit_values))] <- fit_values
        }
        coef_list <- append(coef_list, list(coeff))
      }

      fitted_value_i <- mat
    }

    if (length(weight_list_i) > 0) {
      if (if_col_cor == TRUE) {
        Z_ij_hat_i <- apply(fitted_value_i, 1, function(a) {
          weighted.mean(x = a, w = weight_list_i, na.rm = TRUE)
        })
      } else {
        Z_ij_hat_i <- apply(fitted_value_i, 2, function(a) {
          weighted.mean(x = a, w = weight_list_i, na.rm = TRUE)
        })
      }
    } else {
      Z_ij_hat_i <- NA
    }

    list(Z_ij_hat = Z_ij_hat_i)
  } # nocov end

  stopImplicitCluster()

  Z_ij_hat_mat <- matrix(NA, n_row, n_col)

  for (i in seq_len(iter_over)) {
    result <- results[[i]]
    if (!is.null(result)) {
      if (if_col_cor) {
        Z_ij_hat_mat[, i] <- result$Z_ij_hat
      } else {
        Z_ij_hat_mat[i, ] <- result$Z_ij_hat
      }
    }
  }

  # Step 5: standardize the residuals within each drug column and flag outliers
  R_ij_mat <- Z_ij_mat - Z_ij_hat_mat

  r_ij_mat <- apply(R_ij_mat, 2, function(a) {
    (a - mean(a, na.rm = TRUE)) / sd(a, na.rm = TRUE)
  })

  r_pval <- 1 - pnorm(r_ij_mat)
  colnames(r_pval) <- col_names
  rownames(r_pval) <- row_names

  r_adj_pval <- matrix(p.adjust(r_pval, method = "BH"),
    nrow = n_row,
    ncol = n_col
  )
  colnames(r_adj_pval) <- col_names
  rownames(r_adj_pval) <- row_names

  rslt_list <- list(high_outlier, r_pval, r_adj_pval)
  names(rslt_list) <-
    c("boxplot_signal", "corr_signal_pval", "corr_signal_adj_pval")

  return(rslt_list)
}
