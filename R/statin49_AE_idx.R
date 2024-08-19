#' Cluster index of the FDA statin dataset with 49 adverse events
#'
#' @description A 50 by 2 data frame showing the cluster index of the
#' \code{statin49} dataset's adverse events.
#'
#' @format A data frame with 50 rows and 2 columns: \code{idx} and
#' \code{AE}.
#'
#' @details A data frame with 50 rows and 2 columns: \code{idx} and
#' \code{AE}. \code{AE} lists the adverse event names in the \code{statin49}
#' dataset, while \code{idx} lists the cluster index of each adverse event.
#'
#' The 49 AEs are classified into three clusters: 1) AEs associated with
#' signs and symptoms of muscle injury, 2) AEs associated with laboratory
#' tests for muscle injury, and 3) AEs associated with kidney injury and
#' its laboratory diagnosis and treatment.
#' @examples
#' data(statin49_AE_idx)
#' head(statin49_AE_idx)
"statin49_AE_idx"
