% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BenExtraFcns.R
\name{firstdigitsfreq}
\alias{firstdigitsfreq}
\title{Count first digit frequencies directly from data}
\usage{
firstdigitsfreq(x, numdigits = 1)
}
\arguments{
\item{x}{Numeric matrix or vector.}

\item{numdigits}{A positive integer, default is 1.}
}
\value{
The output will always be a vector of length
  \eqn{9*(10^(numdigits-1))}.
}
\description{
\code{firstdigitsfreq} returns a vector containing the \code{numdigits} first
significant digits frequencies of data \code{x}.
}
\details{
This function takes the user's data \code{x} as an argument, and outputs the
\code{numdigits} first digits frequencies of the data. For example, inputting
\code{x}=\{1,2,30,400\} and \code{numdigits}=1 will return the vector
\{1,1,1,1,0,0,0,0,0\}.

Note that running \code{firstdigitsfreq{x, numdigits}} is equivalent to running:

\code{\link{freqdigits}(\link{extractdigits}(x, numdigits), numdigits)}

Thus zeroes are excluded. For more details on this process,
see \code{\link{extractdigits}} and \code{\link{freqdigits}}.
}
\examples{
firstdigitsfreq(1:9, 1)

}
\references{
Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
  (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
  simultaneous confidence intervals. PLoS one; 11(3).
}
\seealso{
\code{\link{extractdigits}}

\code{\link{freqdigits}}
}
