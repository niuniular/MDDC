#' Report the potential adverse events for drugs from contingency table
#'
#' @param contin_table A data matrix of an \eqn{I} x \eqn{J} contingency
#' table with row (adverse event) and column (drug or vaccine) names.
#' Please first check the input contingency table using the function
#' \code{check_and_fix_contin_table()}.
#' @param contin_table_signal A data matrix with the same dimension and row
#' and column names as \code{contin_table}, with entries either 1
#' (indicating signal) or 0 (indicating non-signal). This data matrix can
#' be obtained via applying the function \code{mddc_boxplot} or
#' \code{mddc_mc}.
#' @param along_rows Specifies the content along the rows of the
#' \code{contin_table} (e.g. AE or Drug).
#' @param along_cols Specifies the content along the columns of the
#' \code{contin_table} (e.g. AE or Drug).
#' @return A data frame with five variables:
#' \itemize{
#' \item \code{drug} the drug name.
#' \item \code{AE} the potential adverse event for the drug or vaccine.
#' \item \code{observed_num} the observed count of the (drug or vaccine, AE)
#' pair.
#' \item \code{expected_num} the expected count of the (drug or vaccine , AE)
#' pair.
#' \item \code{std_pearson_res} the value of the standardized Pearson residual.
#' }
#' @export
#'
#' @examples
#' # load statin49 data
#' data(statin49)
#'
#' # run mddc boxplot method
#' test1 <- mddc_boxplot(statin49)
#'
#' # get the signals from step 2
#' contin_table_signal <- test1$boxplot_signal
#'
#' # get the signals from step 5
#' contin_table_signal_corr <- test1$corr_signal_pval < 0.05
#'
#' # identify the (drug, AE) signals for step 2
#' result_1 <- report_drug_AE_pairs(
#'   contin_table = statin49,
#'   contin_table_signal = contin_table_signal
#' )
#' result_1
#'
#' # identify the (drug, AE) signals for step 5
#' result_2 <- report_drug_AE_pairs(
#'   contin_table = statin49,
#'   contin_table_signal = contin_table_signal_corr
#' )
#' result_2
report_drug_AE_pairs <- function(contin_table, contin_table_signal,
                                 along_rows = "AE", along_cols = "Drug") {
  # Check if the inputs are data matrix
  if (!is.matrix(contin_table) || !is.matrix(contin_table_signal)) {
    stop("Both inputs must be data matrices.")
  }

  # Check if the dimensions match
  if (!all(dim(contin_table) == dim(contin_table_signal))) {
    stop(paste(
      "The dimensions of contin_table and contin_table_signal",
      "must be the same."
    ))
  }

  # Check if the row and column names match
  if (!all(rownames(contin_table) == rownames(contin_table_signal))) {
    stop("The row names of contin_table and contin_table_signal must match.")
  }

  if (!all(colnames(contin_table) == colnames(contin_table_signal))) {
    stop("The column names of contin_table and contin_table_signal must match.")
  }

  if (is.null(rownames(contin_table)) && is.null(colnames(contin_table))) {
    rownames(contin_table) <- paste(along_rows, seq_len(nrow(contin_table)),
      sep = "_"
    )
    colnames(contin_table) <- paste(along_cols, seq_len(ncol(contin_table)),
      sep = "_"
    )
  }

  row_names <- rownames(contin_table)
  col_names <- colnames(contin_table)

  # Replace NA entries in contin_table_signal with 0
  contin_table_signal[is.na(contin_table_signal)] <- 0

  mat_expected_count <- round(get_expected_counts(contin_table), 4)
  mat_std_res <- round(get_std_pearson_res(contin_table), 4)

  pairs <- list()

  for (j in seq_len(ncol(contin_table_signal))) {
    for (i in seq_len(nrow(contin_table_signal))) {
      if (contin_table_signal[i, j] == 1 && contin_table[i, j] !=
        0) {
        pairs <- c(pairs, list(c(col_names[j], row_names[i], contin_table[
          i,
          j
        ], mat_expected_count[i, j], mat_std_res[i, j])))
      }
    }
  }

  if (length(pairs) > 0) {
    pairs_df <- do.call(rbind, pairs)
    colnames(pairs_df) <- c(
      "drug", "AE", "observed_num", "expected_num",
      "std_pearson_res"
    )
    return(as.data.frame(pairs_df))
  } else {
    return(data.frame(Column = character(0), Row = character(0)))
  }
}
