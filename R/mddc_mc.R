#' Modified Detecting Deviating Cells (MDDC) algorithm for adverse event signal identification with Monte Carlo (MC) method for cutoff selection.
#'
#' @description Modified Detecting Deviating Cells (MDDC) algorithm for adverse event signal identification. Monte Carlo (MC) method is used for cutoff selection in step 2 of the algorithm.
#' @param contin_table A data matrix of an \eqn{I} x \eqn{J} contingency table with row (adverse event) and column (drug) names. Please first check the input contingency table using the function \code{check_and_fix_contin_table()}.
#' @param quantile_mc In the step 2 of the algorithm, the quantile of the null distribution obtained via MC method to use as a threshold for identify cells with high value of the standardized Pearson residuals. Default is 0.95.
#' @param mc_num In the step 2, the number of Monte Carlo replications in the MC method. Default is 10000.
#' @param exclude_same_drug_class In the step 2, when applying Fisher's exact test to cells with a count less than six, a 2 by 2 contingency table is needed to be constructed. Does the construction exclude other drugs in the same class as the drug of interest? Default is \code{TRUE}.
#' @param col_specific_cutoff Logical. In the step 2 of the algorithm, whether to apply MC method to the standardized Pearson residuals of the entire table, or within each drug column. Default is \code{TRUE}, that is within each drug column (column specific cutoff). \code{FALSE} indicates applying MC method on residuals of the entire table.
#' @param separate Logical. In the step 2 of the algorithm, whether to separate the standardized Pearson residuals for the zero cells and non zero cells and apply MC method separately or together. Default is \code{TRUE}.
#' @param if_col_cor Logical. In the step 3 of the algorithm, whether to use column (drug) correlation or row (adverse event) correlation. Default is \code{FALSE}, that is using the adverse event correlation. \code{TRUE} indicate using drug correlation.
#' @param cor_lim A numeric value between (0, 1). In the step 3, what correlation threshold should be used to select ``connected'' adverse events. Default is 0.8.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{mc_pval} returns the p values for each cell in the step 2. For cells with count greater than five, the p values are obtained via MC method. For cells with count less than or equal to five, the p values are obtained via Fisher's exact tests.
#' \item \code{mc_signal} returns the signals with a count greater than five and identified in the step 2 by MC method. 1 indicates signals, 0 for non signal.
#' \item \code{fisher_signal} returns the signals with a count less than or equal to five and identified in the step 2 by Fisher's exact tests. 1 indicates signals, 0 for non signal.
#' \item \code{corr_signal_pval} returns the p values for each cell in the contingency table in the step 5, when the \eqn{r_{ij}} values are mapped back to the standard normal distribution.
#' \item \code{corr_signal_adj_pval} returns the Benjamini-Hochberg adjusted p values for each cell in the step 5. We leave here an option for the user to decide whether to use \code{corr_signal_pval} or \code{corr_signal_adj_pval}, and what threshold for p values should be used (for example, 0.05). Please see the example below.
#' }
#' @export
#'
#' @examples
#' # using statin49 data set as an example
#' data(statin49)
#'
#' # apply the mddc_mc
#' eg2 <- mddc_boxplot(statin49)
#'
#' # signals identified in step 2 using MC method
#' signal_step2 <- eg2$mc_signal
#'
#' # signals identified in step 5 by considering AE correlations
#' # In this example, cells with p values less than 0.05 are identified as signals
#' signal_step5 <- (eg2$corr_signal_pval < 0.05) * 1
#' @importFrom grDevices boxplot.stats
#' @importFrom stats cor
#' @importFrom stats weighted.mean
#' @importFrom stats sd
#' @importFrom stats pnorm
#' @importFrom stats p.adjust
#' @importFrom stats rmultinom
#' @importFrom stats fisher.test
#' @importFrom stats quantile
#' @importFrom stats lm

mddc_mc <- function(contin_table,
                    quantile_mc = 0.95,
                    mc_num = 1e4,
                    exclude_same_drug_class = TRUE,
                    col_specific_cutoff = TRUE,
                    separate = TRUE,
                    if_col_cor = FALSE, # if consider column correlations? set it to F since we are using row correlations (AE correlations)
                    cor_lim = 0.8 # c_corr in step 3
) {
  # get p-value

  get_pval <- function(obs, dist) {
    dist <- as.vector(dist)
    dist <- dist[!is.na(dist)]
    pval <- (1 + sum(obs < dist)) / (1 + length(dist))
    return(pval)
  }

  V_get_pval <- Vectorize(get_pval, "obs")

  ### get bootstrap cutoff
  get_Z_ij_mat <- function(contin_table) {
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

    ##    !!!!!! keep only the residual with a cell count > 5 ########
    Z_ij_mat[which(contin_table < 6)] <- NA

    return(Z_ij_mat)
  }

  # get log scale cutoff
  get_log_bootstrap_cutoff <- function(contin_table, quantile_mc) {
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
    p_mat <- p_i_dot %*% t(p_dot_j)

    set.seed(42)

    sim_tab_list <- lapply(
      seq_len(mc_num),
      function(a) {
        matrix(rmultinom(1, n_dot_dot, p_mat), ncol = n_col)
      }
    )

    Z_ij_mat_list <- lapply(sim_tab_list, get_Z_ij_mat)

    max_list <-
      matrix(unlist(lapply(Z_ij_mat_list, function(a) {
        apply(a, 2, function(b) {
          max(log(b), na.rm = TRUE)
        })
      })), byrow = TRUE, ncol = n_col)

    cutoffs <-
      unlist(lapply(seq_len(n_col), function(a) {
        quantile(max_list[, a], quantile_mc, na.rm = TRUE)
      }))
    #  cutoffs <- unlist(lapply(1:n_col, function(a) quantile(max_list[,a], 0.75, na.rm = TRUE)))

    names(cutoffs) <- colnames(contin_table)

    rslt <- list(cutoffs, max_list)
    names(rslt) <- c("cutoffs", "null_dist_s")
    return(rslt)
  }


  cutoffs_and_dist_s <- suppressWarnings(get_log_bootstrap_cutoff(contin_table, quantile_mc))
  c_univ_drug <- cutoffs_and_dist_s[[1]]
  null_dist_s <- cutoffs_and_dist_s[[2]]


  n_row <- nrow(contin_table)
  n_col <- ncol(contin_table)

  row_names <- row.names(contin_table)
  col_names <- colnames(contin_table)

  n_i_dot <- rowSums(contin_table)
  n_dot_j <- colSums(contin_table)
  n_dot_dot <- sum(contin_table)

  p_i_dot <- n_i_dot / n_dot_dot
  p_dot_j <- n_dot_j / n_dot_dot


  # Step 1: for each cell compute the standardized Pearson residual

  E_ij_mat <- n_i_dot %*% t(n_dot_j) / n_dot_dot
  Z_ij_mat <-
    (contin_table - E_ij_mat) / sqrt(E_ij_mat * ((1 - p_i_dot) %*% t((1 - p_dot_j))))


  p_val_mat <- suppressWarnings(matrix(unlist(lapply(seq_len(n_col), function(a) V_get_pval(as.vector(log(Z_ij_mat)[, a]), as.vector(null_dist_s[, a])))),
    ncol = n_col, byrow = FALSE
  ))


  ###### Fisher's exact test for cell <= 5 #########

  get_fisher <- function(row_idx, col_idx) {
    tabl <- matrix(NA, 2, 2)
    tabl[1, 1] <- contin_table[row_idx, col_idx]
    tabl[2, 1] <- sum(contin_table[-row_idx, col_idx])

    if (exclude_same_drug_class == TRUE) {
      tabl[1, 2] <- contin_table[row_idx, n_col]
      tabl[2, 2] <- sum(contin_table[-row_idx, n_col])
    } else {
      tabl[1, 2] <- sum(contin_table[row_idx, -col_idx])
      tabl[2, 2] <- sum(contin_table[-row_idx, -col_idx])
    }
    return(fisher.test(tabl)$p.value)
  }


  for (j in seq_len(n_col)) {
    for (i in which((contin_table[, j] < 6) & (contin_table[, j] > 0))) {
      p_val_mat[i, j] <- get_fisher(i, j)
    }
  }

  ##################

  p_val_mat[is.na(p_val_mat)] <- 1
  row.names(p_val_mat) <- row_names
  colnames(p_val_mat) <- col_names


  signal_mat <- ifelse((p_val_mat < (1 - quantile_mc)) & (contin_table > 5), 1, 0)
  second_signal_mat <- ifelse((p_val_mat < (1 - quantile_mc)) & (contin_table < 6), 1, 0)



  ###################

  res_all <- as.vector(Z_ij_mat)
  res_nonzero <- as.vector(Z_ij_mat[which(contin_table != 0)])
  res_zero <- as.vector(Z_ij_mat[which(contin_table == 0)])


  if (col_specific_cutoff == TRUE) {
    if (separate == TRUE) {
      zero_drug_cutoff <- unlist(lapply(seq_len(n_col), function(a) boxplot.stats(Z_ij_mat[which(contin_table[, a] == 0), a])$stats[[1]]))
    } else {
      zero_drug_cutoff <- apply(Z_ij_mat, 2, function(a) boxplot.stats(a)$stats[[1]])
    }
  } else {
    if (separate == TRUE) {
      zero_drug_cutoff <- rep(boxplot.stats(res_zero)$stats[1], n_col)
    } else {
      zero_drug_cutoff <- rep(boxplot.stats(res_all)$stats[1], n_col)
    }
  }

  # Step 2: apply univariate outlier detection to all the cells

  high_outlier <- matrix(unlist(lapply(seq_len(n_col), function(a) (Z_ij_mat[, a] > exp(c_univ_drug)[a]) * 1)), ncol = n_col)
  colnames(high_outlier) <- colnames(contin_table)
  row.names(high_outlier) <- row.names(contin_table)

  low_outlier <- matrix(unlist(lapply(seq_len(n_col), function(a) (Z_ij_mat[, a] < -exp(c_univ_drug)[a]) * 1)), ncol = n_col)
  colnames(low_outlier) <- colnames(contin_table)
  row.names(low_outlier) <- row.names(contin_table)

  zero_cell_out <- matrix(unlist(lapply(seq_len(n_col), function(a) (Z_ij_mat[, a] < zero_drug_cutoff[a]))), ncol = n_col)
  zero_cell_outlier <- ((zero_cell_out) & (contin_table == 0)) * 1




  if_outlier_mat <-
    ((high_outlier + low_outlier + zero_cell_outlier) != 0) * 1

  U_ij_mat <- ifelse(1 - if_outlier_mat, Z_ij_mat, NA)

  # Step 3 & 4: consider the bivariate relations between AEs and predict values based on the connected AEs

  cor_with_NA <- function(mat, if_col_corr) {
    # this function can be replaced by cor(mat, use = "pairwise.complete.obs")
    if (if_col_corr == FALSE) {
      mat <- t(mat)
    }

    n_col <- ncol(mat)
    cor_mat <- matrix(NA, nrow = n_col, ncol = n_col)
    row.names(cor_mat) <- colnames(mat)
    colnames(cor_mat) <- colnames(mat)

    for (i in seq_len(n_col)) {
      for (j in setdiff(seq_len(n_col), i)) {
        idx <- which((!is.na(mat[, i])) & (!is.na(mat[, j])))
        # cor_mat[i,j] <- cor(mat[idx,i], mat[idx,j])
        cor_mat[i, j] <-
          ifelse((length(idx) >= 3), cor(mat[idx, i], mat[idx, j]), NA)
      }
    }
    return(cor_mat)
  }

  cor_orig <-
    cor(t(contin_table)) # correlation between original cell counts
  cor_Z <-
    cor(t(Z_ij_mat)) # correlation between the standardized Pearson residuals
  cor_U <-
    cor_with_NA(U_ij_mat, if_col_cor) # correlation between the standardized Pearson residuals without the outlying cells
  # cor_U <- cor(t(U_ij_mat), use = "pairwise.complete.obs")

  if (if_col_cor == TRUE) {
    cor_list <- list()
    weight_list <- list()
    fitted_value_list <- list()
    Z_ij_hat_mat <- matrix(NA, n_row, n_col)
    coef_list <- list()

    for (i in seq_len(n_col)) {
      idx <-
        which(abs(cor_U[i, ]) >= cor_lim)
      cor_list[[i]] <- idx[!idx %in% i]
      weight_list[[i]] <- abs(cor_U[i, cor_list[[i]]])

      if (length(cor_list[[i]]) == 0) {
        fitted_value_list[[i]] <- NA
      } else {
        mat <- matrix(NA, n_row, length(cor_list[[i]]))
        row.names(mat) <- row_names
        colnames(mat) <- col_names[cor_list[[i]]]

        for (j in seq_len(length(cor_list[[i]]))) {
          coeff <-
            lm(U_ij_mat[, i] ~ U_ij_mat[, cor_list[[i]][j]])$coefficient
          fit_values <-
            U_ij_mat[, cor_list[[i]][j]] * coeff[2] + coeff[1]
          mat[which(row_names %in% names(fit_values)), j] <-
            fit_values
          coef_list <- append(coef_list, coeff)
        }

        fitted_value_list <- append(fitted_value_list, list(mat))
      }

      if (length(weight_list[[i]]) == 0) {

      } else {
        Z_ij_hat_mat[, i] <- apply(
          fitted_value_list[[i]], 1,
          function(a) {
            weighted.mean(
              x = a,
              w = weight_list[[i]],
              na.rm = TRUE
            )
          }
        )
      }
    }
  } else {
    cor_list <- list()
    weight_list <- list()
    fitted_value_list <- list()
    Z_ij_hat_mat <- matrix(NA, n_row, n_col)
    coef_list <- list()

    for (i in seq_len(n_row)) {
      idx <- which(abs(cor_U[i, ]) >= cor_lim)
      cor_list[[i]] <- idx[!idx %in% i]
      weight_list[[i]] <- abs(cor_U[i, cor_list[[i]]])

      if (length(cor_list[[i]]) == 0) {
        fitted_value_list[[i]] <- NA
      } else {
        mat <- matrix(NA, length(cor_list[[i]]), n_col)
        row.names(mat) <- row_names[cor_list[[i]]]
        colnames(mat) <- col_names

        for (j in seq_len(length(cor_list[[i]]))) {
          coeff <-
            lm(U_ij_mat[i, ] ~ U_ij_mat[cor_list[[i]][j], ])$coefficient
          fit_values <-
            U_ij_mat[cor_list[[i]][j], ] * coeff[2] + coeff[1]
          mat[j, which(col_names %in% names(fit_values))] <-
            fit_values
          coef_list <- append(coef_list, coeff)
        }

        fitted_value_list <- append(fitted_value_list, list(mat))
      }

      if (length(weight_list[[i]]) == 0) {

      } else {
        Z_ij_hat_mat[i, ] <- apply(
          fitted_value_list[[i]], 2,
          function(a) {
            weighted.mean(
              x = a,
              w = weight_list[[i]],
              na.rm = TRUE
            )
          }
        )
      }
    }
  }



  # Step 5: standardize the residuals within each drug column and flag outliers

  R_ij_mat <- Z_ij_mat - Z_ij_hat_mat

  r_ij_mat <-
    apply(R_ij_mat, 2, function(a) {
      (a - mean(a, na.rm = TRUE)) / sd(a, na.rm = TRUE)
    })

  r_pval <- 1 - pnorm(r_ij_mat)
  r_adj_pval <- matrix(p.adjust(r_pval, method = "BH"),
    nrow = n_row,
    ncol = n_col
  )

  list_mat <- list(
    p_val_mat,
    signal_mat,
    second_signal_mat,
    r_pval,
    r_adj_pval
  )

  names(list_mat) <- c(
    "mc_pval",
    "mc_signal",
    "fisher_signal",
    "corr_signal_pval",
    "corr_signal_adj_pval"
  )

  return(list_mat)
}
