#' Extract first digits from user data
#'
#' \code{extractdigits} extracts the first nonzero \code{numdigits}
#' significant digits of the user's data, \code{x}.
#'
#' "0" entries are omitted from the data.
#'
#' The first n significant digits of a non-zero number are the first n non-zero
#' digit of the number. For example, the first two significant digits of "012"
#' are "12".
#'
#' The first n significant digits of a number that has less than m < n digits is
#' given by the m-digit number concatenated with (m-n) zeroes. For example, the
#' first four significant digits of "1.2" are "1200".
#'
#' @return The output is a vector of first \code{numdigits} for nonzero entries
#'   in \code{x} (taken in column-major order if \code{x} is a matrix).
#'
#' @param x Numeric matrix or vector.
#' @param numdigits A positive integer, default is 1.
#'
#' @examples
#' extractdigits(1:10, 1)
#' extractdigits(c(012, 1.2, 1, 0),2)
#'
#' @seealso \code{\link{freqdigits}}
#' @seealso \code{\link{firstdigitsfreq}}
#'
#' @export
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).

extractdigits=function(x,numdigits=1){
  xx=abs(as.numeric(x[x!=0]))
  numdigitsx=abs(as.integer(numdigits))
  b <- floor(log10(xx))
  a <- (xx*10^(numdigitsx-1))/(10^b)
  a <- trunc(a)
  return(a)
}

#' Count frequencies of a significant digits vector
#'
#' \code{freqdigits} computes the frequencies of numbers with \code{numdigits}
#' first significant digits in a vector \code{x} of first digits.
#'
#' \code{freqdigits} may return different vectors for a single \code{x}
#' depending upon the \code{numdigits} chosen. For example,
#' \code{freqdigits(1:9, 1)} will return a vector of ones of length 9, but
#' \code{freqdigits(1:9, 2)} will return a vector of zeroes of length 90
#' because the values in 1:9 all have only 1 digit.
#'
#' If the vector \code{x} contains first digits of varying lengths, then \code{freqdigits}
#' will only count the first digits containing
#' \code{numdigits} digits. For example, inputting \code{x}=\{1,10,100\} with
#' \code{numdigits}=1 will return a vector of length 9 with a 1 in the first
#' position and zeroes otherwise.
#'
#' @return The output will always be a vector of length
#'   \eqn{9*(10^(numdigits-1))}.
#'
#' @param x A significant digits vector.
#' @param numdigits A positive integer, default is 1.
#'
#' @examples
#' freqdigits(1:9, 1)
#' freqdigits(c(456,290,180), 3)
#' freqdigits(extractdigits(c(0.41, 1.25, 0.21, 0.54, 0.19),1),1)
#'
#' @seealso \code{\link{extractdigits}}
#' @seealso \code{\link{firstdigitsfreq}}
#'
#' @export
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
freqdigits=function(x,numdigits=1){
  numdigitsx=abs(as.integer(numdigits))
  digits=seq(from=1*10^(numdigitsx-1),length=9*(10^(numdigitsx-1)))
  freq=tabulate(x,10^numdigitsx-1)
  freq=freq[digits]
  return(freq)
}

#' Count first digit frequencies directly from data
#'
#' \code{firstdigitsfreq} returns a vector containing the \code{numdigits} first
#' significant digits frequencies of data \code{x}.
#'
#' This function takes the user's data \code{x} as an argument, and outputs the
#' \code{numdigits} first digits frequencies of the data. For example, inputting
#' \code{x}=\{1,2,30,400\} and \code{numdigits}=1 will return the vector
#' \{1,1,1,1,0,0,0,0,0\}.
#'
#' Note that running \code{firstdigitsfreq{x, numdigits}} is equivalent to running:
#'
#' \code{\link{freqdigits}(\link{extractdigits}(x, numdigits), numdigits)}
#'
#' Thus zeroes are excluded. For more details on this process,
#' see \code{\link{extractdigits}} and \code{\link{freqdigits}}.
#'
#'
#' @return The output will always be a vector of length
#'   \eqn{9*(10^(numdigits-1))}.
#'
#' @param x Numeric matrix or vector.
#' @param numdigits A positive integer, default is 1.
#'
#' @examples
#' firstdigitsfreq(1:9, 1)
#'
#' @seealso \code{\link{extractdigits}}
#' @seealso \code{\link{freqdigits}}
#'
#' @export
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
firstdigitsfreq=function(x,numdigits=1){
  xx=abs(as.numeric(x[x!=0]))
  numdigitsx=abs(as.integer(numdigits))
  b <- floor(log10(xx))
  a <- (xx*10^(numdigitsx-1))/(10^b)
  a <- trunc(a)
  digits=seq(from=1*10^(numdigitsx-1),length=9*(10^(numdigitsx-1)))
  freq=tabulate(a,10^numdigitsx-1)
  freq=freq[digits]
  return(freq)
}

