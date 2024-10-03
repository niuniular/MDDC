#' Modified Detecting Deviating Cells (MDDC) algorithm for adverse event signal
#' identification with Monte Carlo (MC) method for cutoff selection.
#'
#' @description Modified Detecting Deviating Cells (MDDC) algorithm for adverse
#' event signal identification. Monte Carlo (MC) method is used for cutoff
#' selection in the second step of the algorithm.
#' @param contin_table A data matrix of an \eqn{I} x \eqn{J} contingency table
#' with row (adverse event) and column (drug or vaccine) names.
#' Please first check the input contingency table using the function
#' \code{check_and_fix_contin_table()}.
#' @param quantile In the second step of the algorithm, the quantile of the null
#' distribution obtained via MC method to use as a threshold for identifying
#' cells with high value of the standardized Pearson residuals. Default is 0.95.
#' @param rep In the second step, the number of Monte Carlo replications in the
#' MC method. Default is 10000.
#' @param exclude_same_drug_class In the second step, when applying Fisher's
#' exact test to cells with a count less than six, a 2 by 2 contingency table
#' needs to be constructed. Does the construction need to exclude other drugs
#' or vaccines in the same class as the drug or vaccine of interest?
#' Default is \code{TRUE}.
#' @param col_specific_cutoff Logical. In the second step of the algorithm,
#' whether to apply MC method to the standardized Pearson residuals
#' of the entire table, or within each drug or vaccine column.
#' Default is \code{TRUE}, that is within each drug or vaccine
#' column (column specific cutoff). \code{FALSE} indicates applying MC method
#' on residuals of the entire table.
#' @param separate Logical. In the second step of the algorithm, whether to
#' separate the standardized Pearson residuals for the zero cells and non zero
#' cells and apply MC method separately or together. Default is \code{TRUE}.
#' @param if_col_cor Logical. In the third step of the algorithm, whether to use
#' column (drug or vaccine) correlation or row (adverse event) correlation.
#' Default is \code{FALSE}, that is using the adverse event correlation.
#' \code{TRUE} indicates using drug or vaccine correlation.
#' @param cor_lim A numeric value between (0, 1).
#' In the third step, what correlation threshold should be used to
#' select ``connected'' adverse events. Default is 0.8.
#' @param num_cores Number of cores used to parallelize the
#' MDDC MC algorithm. Default is 2.
#' @param seed An optional integer to set the seed for reproducibility.
#' If NULL, no seed is set.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{mc_pval} returns the p values for each cell in the second step.
#' For cells with a count greater than five, the p values are obtained
#' via MC method. For cells with a count less than or equal to five,
#' the p values are obtained via Fisher's exact tests.
#' \item \code{mc_signal} returns the signals with a count greater than five and
#' identified in the second step by MC method. 1 indicates signals, 0 for non
#' signal.
#' \item \code{fisher_signal} returns the signals with a count
#' less than or equal to five and identified in the second step by
#' Fisher's exact tests. 1 indicates signals, 0 for non signal.
#' \item \code{corr_signal_pval} returns the p values for each cell in the
#' contingency table in the fifth step, when the \eqn{r_{ij}} values are mapped
#' back to the standard normal distribution.
#' \item \code{corr_signal_adj_pval} returns the Benjamini-Hochberg adjusted p
#' values for each cell in the fifth step. We leave here an option for the user
#' to decide whether to use \code{corr_signal_pval} or
#' \code{corr_signal_adj_pval}, and what threshold for p values should be used
#' (for example, 0.05). Please see the example below.
#' }
#' @export
#'
#' @examples
#' # using statin49 data set as an example
#' data(statin49)
#'
#' # apply the mddc_mc
#' mc_res <- mddc_boxplot(statin49)
#'
#' # signals identified in step 2 using MC method
#' signal_step2 <- mc_res$mc_signal
#'
#' # signals identified in step 5 by considering AE correlations
#' # In this example, cells with p values less than 0.05 are
#' # identified as signals
#' signal_step5 <- (mc_res$corr_signal_pval < 0.05) * 1
#' @useDynLib MDDC
#' @importFrom grDevices boxplot.stats
#' @importFrom stats cor
#' @importFrom stats weighted.mean
#' @importFrom stats sd
#' @importFrom stats pnorm
#' @importFrom stats p.adjust
#' @importFrom stats fisher.test
#' @importFrom stats lm
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster


mddc_mc <- function(
    contin_table,
    quantile = 0.95,
    rep = 10000,
    exclude_same_drug_class = TRUE,
    col_specific_cutoff = TRUE,
    separate = TRUE,
    if_col_cor = FALSE,
    cor_lim = 0.8,
    num_cores = 2,
    seed = NULL) {
  V_get_pval <- Vectorize(getPVal, "obs")

  cutoffs_and_dist_s <- suppressWarnings(get_log_bootstrap_cutoff(
    contin_table,
    quantile, rep, seed
  ))
  c_univ_drug <- cutoffs_and_dist_s[[1]]
  null_dist_s <- cutoffs_and_dist_s[[2]]

  n_row <- nrow(contin_table)
  n_col <- ncol(contin_table)

  row_names <- row.names(contin_table)
  col_names <- colnames(contin_table)

  Z_ij_mat <- getZijMat(continTable = contin_table, na = FALSE)
  p_val_mat <- suppressWarnings(matrix(unlist(lapply(
    seq_len(n_col),
    function(a) {
      V_get_pval(as.vector(log(Z_ij_mat)[, a]), as.vector(null_dist_s[
        ,
        a
      ]))
    }
  )), ncol = n_col, byrow = FALSE))


  for (j in seq_len(n_col)) {
    for (i in which((contin_table[, j] < 6) & (contin_table[, j] >
      0))) {
      p_val_mat[i, j] <- get_fisher(
        contin_table,
        i - 1,
        j - 1,
        exclude_same_drug_class
      )
    }
  }

  p_val_mat[is.na(p_val_mat)] <- 1
  row.names(p_val_mat) <- row_names
  colnames(p_val_mat) <- col_names


  signal_mat <- ifelse((p_val_mat < (1 - quantile)) & (contin_table >
    5), 1, 0)
  second_signal_mat <- ifelse((p_val_mat < (1 - quantile)) & (contin_table <
    6), 1, 0)

  res_all <- as.vector(Z_ij_mat)
  res_nonzero <- as.vector(Z_ij_mat[which(contin_table != 0)])
  res_zero <- as.vector(Z_ij_mat[which(contin_table == 0)])


  if (col_specific_cutoff == TRUE) {
    if (separate == TRUE) {
      zero_drug_cutoff <- unlist(lapply(seq_len(n_col), function(a) {
        boxplot.stats(Z_ij_mat[which(contin_table[
          ,
          a
        ] == 0), a])$stats[[1]]
      }))
    } else {
      zero_drug_cutoff <-
        apply(Z_ij_mat, 2, function(a) boxplot.stats(a)$stats[[1]])
    }
  } else {
    if (separate == TRUE) {
      zero_drug_cutoff <- rep(boxplot.stats(res_zero)$stats[1], n_col)
    } else {
      zero_drug_cutoff <- rep(boxplot.stats(res_all)$stats[1], n_col)
    }
  }


  high_outlier <- matrix(unlist(lapply(seq_len(n_col), function(a) {
    (Z_ij_mat[
      ,
      a
    ] > exp(c_univ_drug)[a]) * 1
  })), ncol = n_col)
  colnames(high_outlier) <- colnames(contin_table)
  row.names(high_outlier) <- row.names(contin_table)

  low_outlier <- matrix(unlist(lapply(seq_len(n_col), function(a) {
    (Z_ij_mat[
      ,
      a
    ] < -exp(c_univ_drug)[a]) * 1
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

  cor_U <- pearsonCorWithNA(U_ij_mat, if_col_cor)
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

  R_ij_mat <- Z_ij_mat - Z_ij_hat_mat

  r_ij_mat <- apply(R_ij_mat, 2, function(a) {
    (a - mean(a, na.rm = TRUE)) / sd(a, na.rm = TRUE)
  })

  r_pval <- 1 - pnorm(r_ij_mat)
  r_adj_pval <- matrix(p.adjust(r_pval, method = "BH"),
    nrow = n_row,
    ncol = n_col
  )

  colnames(r_pval) <- col_names
  rownames(r_pval) <- row_names

  colnames(r_adj_pval) <- col_names
  rownames(r_adj_pval) <- row_names

  list_mat <- list(
    p_val_mat, signal_mat, second_signal_mat, r_pval,
    r_adj_pval
  )

  names(list_mat) <- c(
    "mc_pval", "mc_signal", "fisher_signal", "corr_signal_pval",
    "corr_signal_adj_pval"
  )

  return(list_mat)
}
