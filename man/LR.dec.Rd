% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BenFcns.R
\name{LR.dec}
\alias{LR.dec}
\title{Likelihood ratio test: Benford vs. Decreasing probabilities}
\usage{
LR.dec(frequencies, digits = 1:9)
}
\arguments{
\item{frequencies}{Vector of multinomial frequencies.}

\item{digits}{A significant digits vector. If unspecified, default is 1:9.  \code{digits}
must match the frequencies \code{frequencies}.}
}
\value{
The output is a list of 3 items: the likelihood ratio test statistic
  \code{$LR}, p-value \code{$SL}, and a convergence number
  \code{$conv} (0 meaning successful convergence).

  The convergence number indicates convergence of the numerical method
  used to obtain the MLEs of the cell probabilities.
}
\description{
\code{LR.dec} tests whether the first significant digits frequencies of a
sample are compatible with Benford's Law using a likelihood ratio test with
test statistic \eqn{-2 log\lambda}, where \eqn{latexascii}{H_0:} cell
probabilities follow Benford's, vs. \eqn{H_1}: cell probabilities are
decreasing, i.e. \eqn{p_1 \ge p_2 \dots \ge p_9}.
}
\examples{
set.seed(123)
LR.dec((firstdigitsfreq(rnorm(100), 1)))

}
\references{
Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
  (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
  simultaneous confidence intervals. PLoS one; 11(3).
  Wong, S. (2010)  Testing Benford's Law with the first two significant
  digits.  University of Victoria, Master's thesis.
}
