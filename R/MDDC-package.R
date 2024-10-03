#' @name MDDC-package
#' @title Modified Detecting Deviating Cells Algorithm in Pharmacovigilance
#' @description
#' In this package, we present the Modified Detecting Deviating Cells (MDDC)
#' algorithm for adverse event identification.
#' For a certain time period, the spontaneous reports can be extracted from
#' the safety database and depicted as an
#' \eqn{I \times J} contingency table, where:
#' \itemize{
#'   \item \eqn{I} denotes the total number of adverse events (AEs).
#'   \item \eqn{J} denotes the total number of drugs or vaccines.
#' }
#' The cell counts \eqn{n_{ij}} represent the total number of reported cases
#' corresponding to the \eqn{j}-th drug/vaccine and
#' the \eqn{i}-th AE. We are interested in identifying which
#' (AE, drug or vaccine) pairs are signals.
#' Signals refer to potential adverse events that may be caused by a
#' drug or vaccine.
#' In the contingency table setting, signals refer to the cells where
#' \eqn{n_{ij}} is abnormally higher than the expected values.
#'
#' The Detecting Deviating Cells (DDC) algorithm, originally proposed by
#' Rousseeuw and Bossche (2018), was designed for outlier identification in
#' multivariate datasets. However, the original DDC algorithm assumes
#' multivariate normality of the data, with cutoffs based on this assumption.
#' In contrast, the MDDC algorithm is designed for the discrete nature
#' of adverse event data in pharmacovigilance, which clearly do not follow
#' a multivariate normal distribution.
#'
#' Our Modified Detecting Deviating Cells (MDDC) algorithm has the following
#' characteristics:
#' \enumerate{
#'   \item It is easy to compute.
#'   \item It considers AE relationships.
#'   \item It uses data-driven cutoffs.
#'   \item It is independent of the use of ontologies.
#' }
#' The MDDC algorithm consists of five steps. The first two steps identify
#' univariate outliers via cutoffs, and the next three steps
#' evaluate the signals using AE correlations. More details can be found in the
#' \href{https://mddc.readthedocs.io/en/latest/user_guide/mddc_algorithm.html}{MDDC algorithm documentation}.
#'
#' For an introduction to the `MDDC` package, see the vignette:
#' \href{https://niuniular.github.io/MDDC/articles/Introduction_to_MDDC.html}{Usage Examples for MDDC in R}.
#'
#' @details This work has been supported by the Food and Drug Administration
#' and the Kaleida Health Foundation.
#'
#' @author
#' Anran Liu, Raktim Mukhopadhyay, and Marianthi Markatou
#'
#' Maintainer: Anran Liu \email{anranliu@buffalo.edu}
#'
#' @references
#' Liu, A., Mukhopadhyay, R., and Markatou, M. (2024). MDDC:
#' An R and Python package for adverse event identification in
#' pharmacovigilance data. arXiv preprint. arXiv:2410.01168.
#'
#' Liu, A., Markatou, M., Dang, O., and Ball, R. (2024).
#' Pattern discovery in pharmacovigilance through the Modified Detecting
#' Deviating Cells (MDDC) algorithm. Technical Report, Department of
#' Biostatistics, University at Buffalo.
#'
#' Rousseeuw, P. J., and Van den Bossche, W. (2018).
#' Detecting deviating data cells. Technometrics, 60(2), 135-145.
#'
"_PACKAGE"
