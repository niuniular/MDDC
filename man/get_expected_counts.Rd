% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_expected_counts}
\alias{get_expected_counts}
\title{Compute the Expected Count Matrix from a Contingency Table}
\usage{
get_expected_counts(continTable)
}
\arguments{
\item{continTable}{A numeric matrix representing the \eqn{I \times J}
contingency table.}
}
\value{
A numeric matrix of the same dimension as \code{continTable},
containing the expected counts for each cell \eqn{(i, j)} of the
contingency table. The expected counts are based on the row and column
marginal sums of the contingency table.
}
\description{
This function computes the expected counts matrix \eqn{E_{ij}} from a
given \eqn{I \times J} contingency table using the formula:
\deqn{E_{ij} = \frac{n_{i.} n_{.j}}{n_{..}}}
where \eqn{n_{i.}} is the sum of the \eqn{i}-th row, \eqn{n_{.j}} is
the sum of the \eqn{j}-th column, and \eqn{n_{..}} is the total sum
of the table.
}
\examples{
# Create a 6 by 4 contingency table
set.seed(42)
dat_mat <- matrix(rpois(6*4, 20), nrow=6)
dat_mat

# Assign row and column names
contin_table <- check_and_fix_contin_table(dat_mat)
contin_table

# Compute the expected counts of the contingency table
get_expected_counts(contin_table)

}
