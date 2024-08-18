#' FDA Anaphylaxis Dataset
#'
#' A 130 by 4 data matrix from the FDA Adverse Event Reporting System (FAERS)
#' database, processed for 2016. It contains drugs and associated adverse
#' events, focusing on events related to anaphylaxis.
#'
#' @format A data matrix with 130 rows and 4 columns:
#' \describe{
#'   \item{drugs}{130 drugs represented by rows.}
#'   \item{adverse_events}{4 adverse events related to anaphylaxis represented
#'   by columns.}
#' }
#'
#' @details
#' The 130 rows represent different drugs, while the 4 columns correspond to
#' adverse events associated with anaphylaxis. Data is from the FDA Adverse
#' Event Reporting System (FAERS) database, specifically from 2016 reports.
#'
#' The dataset helps identify drugs associated with the 4 anaphylaxis-related
#' adverse events, rather than focusing on identifying adverse events
#' associated with a drug class.
#'
#' @name anaphylaxis
#' @docType data
#' @usage data(anaphylaxis)
#' @examples
#' head(anaphylaxis)
#'
"anaphylaxis"
