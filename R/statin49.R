#' FDA statin dataset with 49 adverse events
#'
#' @description A 50 by 7 data matrix of a contingency table processed from
#' the FDA Adverse Event Reporting System (FAERS) database from Q3 2014 to
#' Q4 2020.
#'
#' @format A data matrix with 50 rows and 7 columns.
#'
#' @details A 50 by 7 data matrix of a contingency table from the FDA Adverse
#' Event Reporting System (FAERS) database, covering Q3 2014 to Q4 2020.
#'
#' The 49 rows represent important adverse events associated with statins,
#' with the final row aggregating the remaining 5,990 events. The 49 AEs are
#' classified into three clusters (see \code{statin49_AE_idx}
#' for cluster indices):
#'
#' 1) AEs associated with muscle injury signs and symptoms, 2) AEs associated
#' with muscle injury lab tests, and 3) AEs associated with kidney injury and
#' its diagnosis and treatment.
#'
#' The 7 columns include six statin medications and an "other drugs" column.
#' Marginal totals for each drug: 197,390 for Atorvastatin,
#' 5,742 for Fluvastatin, 3,230 for Lovastatin, 22,486 for
#' Pravastatin, 122,450 for Rosuvastatin, 85,445 for Simvastatin, and
#' 63,539,867 for Other drugs.
#'
#' @examples
#' data(statin49)
#' head(statin49)
#'
"statin49"
