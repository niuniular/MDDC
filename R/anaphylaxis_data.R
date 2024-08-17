#' FDA Anaphylaxis Dataset
#'
#' A 130 by 4 data matrix of a contingency table processed from the FDA Adverse Event Reporting System (FAERS) database from 2016.
#' The dataset consists of drugs and associated adverse events, particularly focusing on events related to anaphylaxis.
#'
#' @format A data matrix with 130 rows and 4 columns:
#' \describe{
#'   \item{drugs}{130 drugs represented by rows.}
#'   \item{adverse_events}{4 adverse events related to anaphylaxis represented by columns.}
#' }
#'
#' @details
#' The 130 rows represent different drugs, while the 4 columns correspond to adverse events associated with anaphylaxis.
#' The data is derived from the FDA Adverse Event Reporting System (FAERS) database, specifically from reports collected in 2016.
#'
#' The dataset aims to help identify the drugs that are associated with the 4 adverse events related to anaphylaxis, rather than focusing
#' on identifying adverse events associated with a drug class.
#'
#' @name anaphylaxis
#' @docType data
#' @usage data(anaphylaxis)
#' @examples
#' head(anaphylaxis)
#' 
"anaphylaxis"
