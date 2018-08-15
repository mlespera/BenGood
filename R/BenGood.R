#' Tools for Assessing Goodness-of-Fit of Data to Benford's Law
#'
#' Welcome to \code{BenGood}! This package provides easy-to-use tools for
#' generating simultaneous confidence intervals, running Cramer-von Mises
#' and likelihood ratio tests
#' and running power tests to aid in assessing the conformance of your
#' dataset to Benford's Law.
#'
#' @section The main functions are:
#' \itemize{
#'     \item \code{\link{firstdigitsfreq}}
#'     \item \code{\link{Ben.ps}}
#'     \item \code{\link{LR.mat.mult}}
#'     \item \code{\link{LR.dec}}
#'     \item \code{\link{LR.genben}}
#'     \item \code{\link{LR.rod}}
#'     \item \code{\link{SimultConf}}
#'     \item \code{\link{SimultBenp}}
#'     \item \code{\link{plotCI}}
#'     \item \code{\link{CVMStats}}
#'     }
#'
#' @section References:
#'   Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.

#' @docType package
#' @name bengood-package
#' @importFrom stats qchisq pchisq dpois ppois integrate nlminb qnorm
#' @importFrom graphics legend locator plot points title lines

NULL
