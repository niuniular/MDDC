#include <Rcpp.h>
#include <RcppEigen.h>
#include <tuple>
#include <limits>
#include <algorithm>
#include <Eigen/Dense>

// [[Rcpp::depends(RcppEigen)]]
//'
//' Compute the Expected Count Matrix from a Contingency Table
//'
//' This function computes the expected counts matrix \eqn{E_{ij}} from a
//' given \eqn{I \times J} contingency table using the formula:
//' \deqn{E_{ij} = \frac{n_{i.} n_{.j}}{n_{..}}}
//' where \eqn{n_{i.}} is the sum of the \eqn{i}-th row, \eqn{n_{.j}} is
//' the sum of the \eqn{j}-th column, and \eqn{n_{..}} is the total sum
//' of the table.
//'
//' @param continTable A numeric matrix representing the \eqn{I \times J}
//' contingency table.
//'
//' @return A numeric matrix of the same dimension as \code{continTable},
//' containing the expected counts for each cell \eqn{(i, j)} of the
//' contingency table. The expected counts are based on the row and column
//' marginal sums of the contingency table.
//'
//' @examples
//' # Create a 6 by 4 contingency table
//' set.seed(42)
//' dat_mat <- matrix(rpois(6*4, 20), nrow=6)
//' dat_mat
//'
//' # Assign row and column names
//' contin_table <- check_and_fix_contin_table(dat_mat)
//' contin_table
//'
//' # Compute the expected counts of the contingency table
//' get_expected_counts(contin_table)
//'
//' @useDynLib MDDC
//' @importFrom Rcpp evalCpp
//' @rdname get_expected_counts
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_expected_counts(const Eigen::MatrixXd &continTable)
{
    // This calculates the expected count matrix Eij

    int nRow = continTable.rows();
    int nCol = continTable.cols();
    Eigen::VectorXd niDot = continTable.rowwise().sum();
    Eigen::VectorXd nDotj = continTable.colwise().sum();
    double nDotDot = continTable.sum();

    Eigen::VectorXd piDot = niDot / nDotDot;
    Eigen::VectorXd pDotj = nDotj / nDotDot;

    Eigen::MatrixXd EijMat = (niDot * nDotj.transpose()) / nDotDot;

    return EijMat;
}

// [[Rcpp::depends(RcppEigen)]]
//'
//' Standardized Pearson Residuals for Contingency Tables
//'
//' Compute the standardized Pearson residuals \eqn{Z_{ij}} for a given
//' \eqn{I \times J} contingency table. Standardized Pearson residuals
//' are calculated as \eqn{Z_{ij} = \frac{(O_{ij} - E_{ij})}{\sqrt{E_{ij}
//' (1 - p_{i.}) (1 - p_{.j})}}} where \eqn{O_{ij}} is the observed count,
//' \eqn{E_{ij}} is the expected count, and \eqn{p_{i.}} and \eqn{p_{.j}}
//' are the row and column marginal probabilities, respectively.
//'
//' @param continTable A numeric matrix representing the \eqn{I \times J}
//' contingency table.
//' @param na A boolean flag (default is \code{TRUE}) to replace cells with
//' observed counts less than 6 with \code{NA} values.
//'
//' @return A matrix containing:
//' \itemize{
//'   \item A matrix of standardized Pearson residuals (\code{ZijMat}).
//'   \item The total sum of the contingency table (\code{nDotDot}).
//'   \item A vector of row marginal proportions (\code{piDot}).
//'   \item A vector of column marginal proportions (\code{pDotj}).
//' }
//'
//' @useDynLib MDDC
//' @rdname StandardizedPearsonResiduals
//' @keywords internal
//'
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd getZijMat(const Eigen::MatrixXd &continTable, bool na = true)
{
    // This calculates the standardized Pearson residuals Zij
    int nRow = continTable.rows();
    int nCol = continTable.cols();
    Eigen::VectorXd niDot = continTable.rowwise().sum();
    Eigen::VectorXd nDotj = continTable.colwise().sum();
    double nDotDot = continTable.sum();

    Eigen::VectorXd piDot = niDot / nDotDot;
    Eigen::VectorXd pDotj = nDotj / nDotDot;

    Eigen::MatrixXd EijMat = (niDot * nDotj.transpose()) / nDotDot;

    Eigen::MatrixXd temp = (1 - piDot.array()).matrix() * (1 - pDotj.array()).transpose().matrix();
    Eigen::MatrixXd sqrtEijMat = EijMat.cwiseProduct(temp).array().sqrt();

    Eigen::MatrixXd ZijMat = (continTable - EijMat).array() / sqrtEijMat.array();

    if (na)
    {
        ZijMat = (continTable.array() < 6).select(std::numeric_limits<double>::quiet_NaN(), ZijMat);
    }

    return ZijMat;
}

// [[Rcpp::depends(RcppEigen)]]
//'
//' Compute P-values Based on an Empirical Distribution
//'
//' This function calculates the p-values for each observation in a vector
//' \code{obs} based on a sorted empirical distribution vector \code{dist}.
//' The p-value is calculated as the proportion of values in \code{dist} that
//' are greater than the observation, adjusted for ranking. If an observation
//' is \code{NaN}, the corresponding p-value will be \code{NaN}.
//'
//' @param obs A numeric vector (Eigen::VectorXd) containing the observations
//' for which p-values are to be computed.
//' @param dist A numeric vector (Eigen::VectorXd) representing the empirical
//' distribution from which the p-values will be calculated.
//'
//' @return A numeric vector (Eigen::VectorXd) of p-values corresponding to
//' each observation in \code{obs}. The p-values are based on the number of
//' elements in \code{dist} greater than each observation.
//'
//' @useDynLib MDDC
//' @rdname getPVal
//' @keywords internal
//'
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd getPVal(const Eigen::VectorXd &obs, const Eigen::VectorXd &dist)
{
    // This calculates the p values

    // Create a sorted copy of dist
    Eigen::VectorXd sortedDist = dist;
    std::sort(sortedDist.data(), sortedDist.data() + sortedDist.size());

    // Initialize the p-values vector
    Eigen::VectorXd pVals = Eigen::VectorXd::Constant(obs.size(), std::numeric_limits<double>::quiet_NaN());

    // Convert Eigen vectors to std::vectors for use with STL algorithms
    std::vector<double> obsVec(obs.data(), obs.data() + obs.size());
    std::vector<double> sortedDistVec(sortedDist.data(), sortedDist.data() + sortedDist.size());

    // Compute p-values for each observation
    std::transform(obsVec.begin(), obsVec.end(), pVals.data(),
                   [&sortedDistVec, distSize = dist.size()](double obs_i)
                   {
                       if (std::isnan(obs_i))
                           return std::numeric_limits<double>::quiet_NaN();
                       auto it = std::upper_bound(sortedDistVec.begin(), sortedDistVec.end(), obs_i);
                       int count = std::distance(it, sortedDistVec.end());
                       return static_cast<double>(1 + count) / (1 + distSize);
                   });

    return pVals;
}

// [[Rcpp::depends(RcppEigen)]]
//'
//' Generate Contingency Table for Fisher's Exact Test
//'
//' This function constructs a 2x2 contingency table used for performing
//' Fisher's Exact Test. The table is derived from the given contingency
//' matrix \code{continTable} based on the specified row and column indices.
//' The user has the option to exclude data from the same drug class by
//' adjusting the calculations for the contingency table.
//'
//' @param continTable A numeric matrix (Eigen::MatrixXd) representing the
//' contingency table from which the 2x2 table will be constructed.
//' @param rowIdx An integer representing the row index of the element to be
//' used in the contingency table.
//' @param colIdx An integer representing the column index of the element to
//' be used in the contingency table.
//' @param excludeSameDrugClass A boolean value indicating whether to exclude
//' the same drug class in the calculations. If \code{true}, the function
//' will adjust the calculations by excluding the same drug class.
//'
//' @return A 2x2 numeric matrix (Eigen::MatrixXd) representing the contingency
//' table for Fisher's Exact Test.
//'
//' @useDynLib MDDC
//' @rdname getFisherExactTestTable
//' @keywords internal
//'
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd getFisherExactTestTable(Eigen::MatrixXd continTable, int rowIdx, int colIdx, bool excludeSameDrugClass)
{
    // This generates the tables used for Fisher Exact Test
    Eigen::MatrixXd tabl(2, 2);

    // Set the values of tabl
    tabl(0, 0) = continTable(rowIdx, colIdx);
    tabl(1, 0) = continTable.col(colIdx).sum() - continTable(rowIdx, colIdx);

    if (excludeSameDrugClass)
    {
        int n_col = continTable.cols() - 1;
        tabl(0, 1) = continTable(rowIdx, n_col);
        tabl(1, 1) = continTable.col(n_col).sum() - continTable(rowIdx, n_col);
    }
    else
    {
        tabl(0, 1) = continTable.row(rowIdx).sum() - continTable(rowIdx, colIdx);
        tabl(1, 1) = continTable.sum() - continTable.row(rowIdx).sum() - continTable.col(colIdx).sum() + continTable(rowIdx, colIdx);
    }

    return tabl;
}

// [[Rcpp::depends(RcppEigen)]]
//' Compute Correlation Matrix with Handling of Missing Values
//'
//' This function calculates a correlation matrix from the input matrix,
//' handling missing values (`NA`) by computing correlations only for
//' non-missing pairs. The function supports computing either row-wise or
//' column-wise correlations, based on the input flag `if_col_cor`.
//'
//' @param mat A numeric matrix or data frame. Each column represents a variable,
//'   and each row represents an observation. The function will compute pairwise
//'   correlations between the columns (or rows if `if_col_cor = FALSE`).
//' @param if_col_cor Logical. If `TRUE`, the function computes correlations
//'   between columns of `mat`. If `FALSE`, correlations are computed between
//'   rows (the matrix is transposed internally).
//'
//' @return A square matrix containing the computed correlation coefficients.
//'   The dimensions of the returned matrix will be equal to the number of
//'   columns (or rows if `if_col_cor = FALSE`) in the input matrix `mat`.
//'   Missing values (`NA`) in the input matrix are handled by only using
//'   observations where both variables have non-missing values. If fewer than
//'   3 valid pairs are available for a correlation, `NA` is returned for
//'   that pair.
//'
//' @details The function iterates through each pair of columns (or rows)
//'   and computes the correlation using Pearson's method, only including
//'   observations where both variables have non-missing values. If the number
//'   of valid pairs is less than 3, the function assigns `NA` to the
//'   corresponding entry in the correlation matrix.
//'
//' @useDynLib MDDC
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd pearsonCorWithNA(const Eigen::MatrixXd &mat, bool ifColCorr = true)
{
    // Computes the Pearson Correlation in the presence of NA values
    Eigen::MatrixXd matrix = ifColCorr ? mat : mat.transpose();
    int nCol = matrix.cols();
    Eigen::MatrixXd corMat(nCol, nCol);
    corMat.fill(std::numeric_limits<double>::quiet_NaN());

    for (int i = 0; i < nCol; ++i)
    {
        for (int j = i + 1; j < nCol; ++j)
        {
            std::vector<double> x, y;
            for (int k = 0; k < matrix.rows(); ++k)
            {
                if (!std::isnan(matrix(k, i)) && !std::isnan(matrix(k, j)))
                {
                    x.push_back(matrix(k, i));
                    y.push_back(matrix(k, j));
                }
            }
            if (x.size() >= 3)
            {
                double meanX = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
                double meanY = std::accumulate(y.begin(), y.end(), 0.0) / y.size();

                double numerator = 0.0;
                double denomX = 0.0;
                double denomY = 0.0;

                for (size_t k = 0; k < x.size(); ++k)
                {
                    double diffX = x[k] - meanX;
                    double diffY = y[k] - meanY;
                    numerator += diffX * diffY;
                    denomX += diffX * diffX;
                    denomY += diffY * diffY;
                }

                double denominator = std::sqrt(denomX * denomY);
                double correlation = numerator / denominator;

                {
                    corMat(i, j) = correlation;
                    corMat(j, i) = correlation;
                }
            }
        }
    }

    return corMat;
}
