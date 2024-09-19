#' Find Adaptive Boxplot Coefficient `coef` via Grid Search
#'
#' This function performs a grid search to determine the optimal
#' adaptive boxplot coefficient `coef` for each column of a contingency
#' table, ensuring the target false discovery rate (FDR) is met.
#'
#' @param contin_table A matrix representing the \eqn{I \times J}
#' contingency table.
#' @param n_sim An integer specifying the number of simulated tables under
#' the assumption of independence between rows and columns. Default is 1000.
#' @param target_fdr A numeric value specifying the desired level of false
#' discovery rate (FDR). Default is 0.05.
#' @param grid A numeric value representing the size of the grid added to
#' the default value of \code{coef = 1.5} as suggested by Tukey. Default is 0.1.
#' @param col_specific_cutoff Logical. If TRUE, then a single value of the
#' coefficient is returned for the entire dataset, else when FALSE specific
#' values corresponding to each of the columns are returned.
#' @param exclude_small_count A logical indicating whether to exclude cells
#' with counts smaller than or equal to five when computing boxplot statistics.
#' Default is \code{TRUE}.
#'
#' @return A list with the following components:
#' \describe{
#'  \code{coef}: A numeric vector containing the optimal coefficient
#'              `coef` for each column of the input contingency table.
#'
#'  \code{FDR}: A numeric vector with the corresponding false discovery
#'              rate (FDR) for each column.
#' }
#' @export
#'
#' @examples
#' \donttest{
#' # This example uses the statin49 data
#' data(statin49)
#' find_optimal_coef(statin49)
#' }
#' @useDynLib MDDC
#' @importFrom stats rmultinom
find_optimal_coef <- function(contin_table,
                              n_sim = 1000,
                              target_fdr = 0.05,
                              grid = 0.1,
                              col_specific_cutoff = TRUE,
                              exclude_small_count = TRUE) {
  I <- nrow(contin_table)
  J <- ncol(contin_table)

  expected <- get_expected_counts(contin_table)

  tab_list <- lapply(
    1:n_sim,
    function(a) {
      matrix(
        rmultinom(1,
          sum(contin_table),
          prob = as.vector(expected)
        ),
        nrow = I
      )
    }
  )

  if (exclude_small_count == FALSE) {
    res_list <- lapply(
      1:n_sim,
      function(a) {
        temp_tab <- tab_list[[a]]
        res <- get_std_pearson_res(temp_tab)
        res[which(temp_tab == 0)] <- NA
        return(res)
      }
    )
  } else {
    res_list <- lapply(
      1:n_sim,
      function(a) {
        temp_tab <- tab_list[[a]]
        res <- get_std_pearson_res(temp_tab)
        res[which(temp_tab <= 5)] <- NA
        return(res)
      }
    )
  }

  c_vec <- rep(1.5, J)
  c <- 1.5

  if (col_specific_cutoff == TRUE) {
    # Optimize c_j for each column j
    fdr_vec <- rep(1, J)

    for (j in 1:J) {
      while (fdr_vec[j] > target_fdr) {
        c_vec[j] <- c_vec[j] + grid
        fdr_vec[j] <- compute_fdr(res_list, c_vec[j], j)
      }
    }
    rslt <- list(c_vec, fdr_vec)
    names(rslt) <- c("coef", "FDR")
  } else {
    fdr <- 1
    while (fdr > target_fdr) {
      c <- c + grid
      fdr <- compute_fdr_all(res_list, c)
    }
    rslt <- list(c, fdr)
    names(rslt) <- c("coef", "FDR")
  }
  return(rslt)
}
