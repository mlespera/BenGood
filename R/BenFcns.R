#' Generate Benford probabilities
#'
#' \code{Ben.ps} returns a vector of Benford probabilities corresponding to the
#' digits of the argument vector \code{digits}.
#'
#' @param digits A significant digits vector, e.g. 1:9, 10:99
#'
#' @examples
#' Ben.ps(1:9)
#'
#' @export
Ben.ps <- function(digits)
{
	return(log10(1+1/digits))
}

#' aLoga
#'
#' \code{aLoga} is used in \code{\link{LR.mat.mult}}.
#'
#' @param a A non-negative number.
#'
#' @keywords internal
aLoga <- function(a)
{
	return(a*log(a+(a==0)))
}

#' bLoga
#'
#' \code{bLoga} is used in \code{\link{LR.genben}}.
#'
#' @param a A non-negative number.
#' @param b A non-negative number.
#'
#' @keywords internal
bLoga <- function(b, a)
{
	return(b*log(a+(a==0 & b==0)))
}


#' Div
#'
#' \code{Div} is used to prevent division by zero in \code{\link{CVMEigen}}.
#'
#' @param a A non-negative number.
#' @param b A non-negative number.
#'
#' @keywords internal
Div <- function(a, b)
{
 	return((b!=0)*a/(b+(b==0)))
}


#' Likelihood ratio test: Benford vs. General Multinomial
#'
#' \code{LR.mat.mult} is a likelihood ratio test for whether the first
#' significant digits of a set of samples are compatible with Benford's Law,
#' where \eqn{H_0:} cell probabilities follow Benford's, vs. \eqn{H_1}: the most
#' general multinomial, i.e. all probabilities non-negative,
#' sum to 1.
#'
#' The LR statistic \eqn{latexascii}{\lambda} for testing \eqn{H_0} vs.
#' \eqn{H_1} is:
#'
#' \eqn{latexascii}{-2log\lambda = -2 \sum f_i{log(p_i) - log(f_i/n)}}
#'
#' where \eqn{f_i} is the frequency of the ith first digit, and \eqn{p_i} is the
#' cell probability of the ith first digit.
#'
#' @return The output is a vector of p-values, one per sample (row).
#'
#' @param freq Vector or matrix with multinomial samples in the rows.
#' @param digits A significant digits vector. If unspecified, default is 1:9. \code{digits}
#' must match the frequencies in the rows of \code{freq}.
#'
#' @examples
#' LR.mat.mult(firstdigitsfreq(c(0.41, 1.25, 0.21, 0.54, 0.19), 1))
#' set.seed(123)
#' LR.mat.mult(firstdigitsfreq(rnorm(100), 1))
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @export
LR.mat.mult <- function(freq, digits=1:9)
{
  if(is.vector(freq)) freq <- matrix(freq,nrow=1)
	log.l1 <- apply(aLoga(freq), 1, sum) - aLoga(apply(freq, 1, sum))
	log.l0 <- apply(t(freq)*log(Ben.ps(digits)), 2, sum)
	LR <- 2*(log.l1-log.l0)
	return( 1 - pchisq(LR, length(digits)-1))
}


#' Likelihood ratio test: Benford vs. Decreasing probabilities
#'
#' \code{LR.dec} tests whether the first significant digits frequencies of a
#' sample are compatible with Benford's Law using a likelihood ratio test with
#' test statistic \eqn{-2 log\lambda}, where \eqn{latexascii}{H_0:} cell
#' probabilities follow Benford's, vs. \eqn{H_1}: cell probabilities are
#' decreasing, i.e. \eqn{p_1 \ge p_2 \dots \ge p_9}.
#'
#' @return The output is a list of 3 items: the likelihood ratio test statistic
#'   \code{$LR}, p-value \code{$SL}, and a convergence number
#'   \code{$conv} (0 meaning successful convergence).
#'
#'   The convergence number indicates convergence of the numerical method
#'   used to obtain the MLEs of the cell probabilities.
#'
#' @param frequencies Vector of multinomial frequencies.
#' @param digits A significant digits vector. If unspecified, default is 1:9.  \code{digits}
#' must match the frequencies \code{frequencies}.
#'
#' @examples
#' set.seed(123)
#' LR.dec((firstdigitsfreq(rnorm(100), 1)))
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @export
LR.dec <- function(frequencies, digits=1:9)
{
  freq <- as.vector(frequencies)
  k <- length(digits)
  log.l1 <-  nlminb(objective=dec.objf, start = rep(1, k), lower = rep(0, k),
                    freq = freq)
  log.l0 <- sum(freq * log(Ben.ps(digits)))
  LR <- 2 * (-log.l1$objective - log.l0)
  return(list(LR=LR, SL= 1 - pchisq(LR, k-1),
              conv=log.l1$convergence))
}


#' Objective function for numerically finding MLEs under LR.dec, given z's
#'
#' \code{dec.objf} contains the (negative) objective function of the nonlinear
#' programming problem used to compute the MLEs in \code{LR.dec}.
#'
#' \code{LR.dec} computes the  MLEs of the first significant digit probabilities
#' numerically under the alternative \eqn{H_1: p_1 \ge p_2 \dots \ge p_9} by
#' solving the nonlinear programming problem:
#'
#' \eqn{ maximize \sum f_i*log(\sum z_j)}
#'
#' subject to constraints, where \eqn{f_i} is the frequency of the ith
#' significant digit, and \eqn{z_i = p_i - p_{i+1}}. \code{dec.objf}
#' contains the negative of this expression in order to adhere to standard
#' minimization optimization procedure. See Lesperance et. al for
#' further details.
#'
#' @param zs Numeric vector, defined as the difference in probabilities:
#'   \eqn{z_i = p_i - p_{i+1}}
#' @param freq Vector of multinomial frequencies.
#'
#' @seealso \code{\link{LR.dec}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
dec.objf <- function(zs, freq)
{
	ps <- rev.cumsum(zs)
	ps <- ps/sum(ps)
 	return(- sum(bLoga(freq, ps)))
}


#' Compute reverse cumulative sum
#'
#' \code{rev.cumsum} is used in \code{dec.objf}.
#'
#' @param v Numeric vector.
#'
#' @seealso \code{\link{dec.objf}}
#'
#' @keywords internal
rev.cumsum <- function(v)
{
	return(rev(cumsum(rev(v))))
}

#' Likelihood ratio test: Benford vs. Generalized Benford
#'
#' \code{LR.genben} tests whether the first significant digits frequencies of a
#' sample are compatible with Benford's Law using a likelihood ratio test, where
#' \eqn{latexascii}{H_0:} cell probabilities follow Benford's, vs. \eqn{H_1}:
#' cell probabilities follow generalized Benford.
#'
#' The generalized Benford's law depends upon parameter \eqn{\alpha}. Because
#' the generalized Benford's law reduces to Benford's law when \eqn{\alpha = 0},
#' equivalent hypotheses for this test are \eqn{latexascii}{H_0: \alpha = 0} vs.
#' \eqn{H_1: \alpha =/= 0}.
#'
#' @param frequencies Vector of multinomial frequencies.
#' @param digits A significant digits vector. If unspecified, default is 1:9.  \code{digits}
#' must match the frequencies \code{frequencies}.
#'
#' @return The output is a list containing a p-value \code{$SL} and a
#'   convergence number \code{$conv} (0 meaning successful convergence)
#'   respectively.
#'
#'   The convergence number indicates convergence of the numerical method
#'   used to obtain the MLE of \eqn{\alpha}.
#'
#' @examples
#' set.seed(123)
#' LR.genben((firstdigitsfreq(rnorm(100), 1)))
#'
#' @seealso \code{\link{genben.ps}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @export
LR.genben <- function(frequencies, digits=1:9)
{
  freq <- as.vector(frequencies)
  digs <- as.vector(digits)
	log.l1 <-  nlminb(objective=genben.ml.objf, start = 0.2, freq = freq, digits=digs )
	log.l0 <- sum(freq * log(Ben.ps(digits)))
	LR <- 2 * (-log.l1$objective - log.l0)
	return(list(SL=1 - pchisq(LR, 1),conv=log.l1$conv))
}


#' Objective function for numerically finding MLE under LR.genben
#'
#' \code{genben.ml.objf} contains the (negative) objective function of the
#' nonlinear programming problem used to compute the MLE in \code{LR.genben}.
#'
#' \code{LR.genben} computes the  MLE for \eqn{\alpha} numerically under the
#' generalized Benford alternative by maximizing a log-likelihood objective
#' function, or equivalently minimizing the negative of said function. See
#' Lesperance et. al for further details.
#'
#' @param alpha A real number.
#' @param freq Vector of multinomial frequencies.
#' @param digits A significant digits vector.
#'
#' @seealso \code{\link{LR.genben}}
#' @seealso \code{\link{genben.ps}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
genben.ml.objf <- function(alpha, freq, digits)
{
 	return(- sum(bLoga(freq,genben.ps(alpha,digits))))
}


#' Generate Generalized Benford probabilities
#'
#' \code{genben.ps} returns a vector of generalized Benford probabilities
#' corresponding to the digits of the argument vector \code{digits}, and
#' the chosen \eqn{\alpha}.
#'
#' @param alpha A real number.
#' @param digits A significant digits vector, default is 1:9.
#'
#' @examples
#' genben.ps(1, 1:9)
#'
#' @seealso \code{\link{Ben.ps}}
#' @seealso \code{\link{LR.genben}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @export
genben.ps <- function(alpha, digits=1:9)
{
	return((digits^(-alpha)-(digits+1)^(-alpha))/(min(digits)^(-alpha)-max(digits + 1)^(-alpha)))
}


#' Likelihood ratio test: Benford vs. Rodriguez family
#'
#' \code{LR.genben} tests whether the first significant digits frequencies of a
#' sample are compatible with Benford's Law using a likelihood ratio test, where
#' \eqn{H_0:} cell probabilities follow Benford's, vs. \eqn{H_1}:
#' cell probabilities follow a Rodriguez's generalization of Benford.
#'
#' The Rodriguez formula depends upon parameter \eqn{\beta}. Because the
#' Rodriguez formula reduces to Benford's law when \eqn{\beta = -1}, equivalent
#' hypotheses for this test are \eqn{latexascii}{H_0: \beta = -1} vs.
#' \eqn{H_1: \beta =/= -1}.
#'
#' @param frequencies Vector of multinomial frequencies.
#' @param digits A significant digits vector. If unspecified, default is 1:9. \code{digits}
#' must match the frequencies \code{frequencies}.
#'
#' @examples
#' set.seed(123)
#' LR.rod((firstdigitsfreq(rnorm(100), 1)))
#'
#' @return The output is a list containing a p-value \code{$SL} and a
#'   convergence number \code{$conv} (0 meaning successful convergence)
#'   respectively.
#'
#'   The convergence number indicates convergence of the numerical method
#'   used to obtain the MLE of \eqn{\alpha}.
#'
#' @seealso \code{\link{rodben.ps}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @export
LR.rod <- function(frequencies, digits=1:9)
{
  freq <- as.vector(frequencies)
  digs <- as.vector(digits)
	log.l1 <-  nlminb(objective=rodben.ml.objf, start = -0.5, freq = freq,
              	digits=digs, control=list(trace=0))

	if(log.l1$par>40)
	{
		log.l1$objective <- -(sum(aLoga(freq))-aLoga(sum(freq)))
    		log.l1$convergence<-2
	}
	log.l0 <- sum(freq * log(Ben.ps(digits)))
	LR <- 2 * (-log.l1$objective - log.l0)
	return(list(SL=1 - pchisq(LR, 1),conv=log.l1$convergence))
}

#' Objective function for numerically finding MLE under LR.rod
#'
#' \code{rodben.ml.objf} contains the objective function of the nonlinear programming
#' problem used to compute the MLE in \code{LR.rod}.
#'
#' \code{LR.rod} computes the  MLE for \eqn{\beta} numerically under the
#' Rodriguez generalization of Benford by maximizing a log-likelihood objective
#' function, or equivalently minimizing the negative of said function. See
#' Lesperance et. al for further details.
#'
#' @param beta A real number.
#' @param freq Vector of multinomial frequencies.
#' @param digits A significant digits vector.
#'
#' @seealso \code{\link{LR.rod}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
rodben.ml.objf <- function(beta, freq, digits)
{
	return( - sum(bLoga(freq, rodben.ps(beta, digits))))
}


#' Generate Rodriguez probabilities
#'
#' \code{rodben.ps} returns a vector of Rodriguez family probabilities
#' corresponding to the digits of the argument vector \code{digits}, and
#' the chosen \eqn{\beta}.
#'
#' @param beta A real number.
#' @param digits A significant digits vector.
#'
#' @seealso \code{\link{LR.rod}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @export
rodben.ps <- function(beta, digits)
{
	return((beta + 1)/(length(digits) * beta) - ((digits + 1)^(beta + 1) - digits^(beta + 1))/
	    (beta * (max(digits + 1)^(beta + 1) - min(digits)^(beta + 1))))
}


#' Compute Cramer-von Mises Statistics
#'
#' \code{CVMStats} computes p-values for the CVM statistics for testing
#' \eqn{H_0:} probabilities follow Benford vs. the most general alternative,
#' \eqn{H_1:} probabilities are multinomial. The four statistics are \eqn{W^2}
#' (Cramer-von Mises), \eqn{U^2} (Watson), \eqn{A^2} (Anderson-Darling) and
#' \eqn{X^2} (Pearson's chi-square).
#'
#' There are three CVM type statistics: Cramer-von Mises as
#' \code{W}, Watson as \code{U}, Anderson-Darling as \code{A}, as
#' well as Pearson's chi-square as \code{X}.
#'
#' The three CVM type statistics can be computed using one of two
#' methods: Imhof's numerical method, or a chi-square approximation. The
#' argument \code{method} can be used to indicate which statistics to compute
#' with which method. Note that the chi-square approximation is faster to
#' compute than the Imhof numerical method. See Lesperance et. al (2016) for further
#' details.
#'
#' @return The output is a vector containing the user specified statistics by
#'   method.
#'
#' @param freq Vector of multinomial probabilities.
#' @param method Parameter to specify the statistic and the method by which to
#'   compute the statistic. This parameter can take on any of the following
#'   values, or a vector containing any combination of these values separated
#'   using commas:
#'
#' \enumerate{
#'     \item 'W'
#'     \item 'Wap'
#'     \item 'U'
#'     \item 'Uap'
#'     \item 'A'
#'     \item 'Aap'
#'     \item 'X'
#'     }
#'
#'     These parameter values correspond to the statistics and their methods
#'     respectively:
#'
#'      \enumerate{
#'      \item Cramer-von Mises, Imhof
#'      \item Cramer-von Mises, approximation
#'      \item Watson, Imhof
#'      \item Watson, approximation
#'      \item Anderson-Darling, Imhof
#'      \item Anderson-Darling, approximation
#'      \item Pearson
#'      }
#'
#'      Specifying no method will return all statistics and methods.
#'
#' @param digits A significant digits vector. If unspecified, default is 1:9. \code{digits}
#' must match the frequencies \code{frequencies}.
#'
#'
#' @examples
#' set.seed(123)
#' CVMStats((firstdigitsfreq(rnorm(100), 1)))

#' @seealso \code{\link{CVMEigen}}
#'
#' @export
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
CVMStats <- function(freq, method = c('W','Wap','U','Uap','A','Aap','X'), digits=1:9)
{
	methods = c('W','Wap','U','Uap','A','Aap','X')
	meth = match.arg(method, methods, several.ok=TRUE)
#	CVME <- CVMEigen(n=100, digits)
	CVME <- CVMEigen(digits)
	CR.eig.W <- CVME$eig.W
	CR.eig.U <- CVME$eig.U
	CR.eig.A <- CVME$eig.A

	k <- length(freq)
	n <- sum(freq)
	ps <- Ben.ps(digits)
	e <- n*ps
	S <- cumsum(freq)
	Tt <- cumsum(e)
	Z <- S-Tt
	H <- Tt/n
	tj <- .5*(ps + c(ps[-1],ps[1]))
	Zbar <- sum(tj*Z)
	W2 <- sum(tj*Z^2)/n
	W <- Imhof(W2, CR.eig.W)
	if (W<0)  {W <- Imhof(W2, CR.eig.W, UPPER=Inf, subdiv=500)}
	Wap <- Chi.approx(W2, CR.eig.W)
	U2 <- sum(tj*(Z-Zbar)^2)/n
	U <- Imhof(U2, CR.eig.U)
	if (U<0)  {U <- Imhof(U2, CR.eig.U, UPPER=Inf, subdiv=500)}
	Uap <- Chi.approx(U2, CR.eig.U)
	A2 <- sum(Div(tj*Z^2,H*(1-H)))/n
	A <- Imhof(A2, CR.eig.A)
	if (A<0)  {A <- Imhof(A2, CR.eig.A, UPPER=Inf, subdiv=500)}
	Aap <- Chi.approx(A2, CR.eig.A)
	X2 <- sum(Div((freq-e)^2, e))
	X <- 1-pchisq(X2, k-1)
  Sigs <- c(W=W, Wap=Wap, U=U, Uap=Uap, A=A, Aap=Aap, X=X)
	m = is.element(methods,meth)
	return(Sigs[m])
}


#' Compute eigenvalues required for function CVMStats
#'
#' \code{CVMEigen} computes eigenvalues needed for computation of significance
#' levels of Cramer-von Mises type statistics for testing \eqn{H_0:} Benford vs
#' \eqn{H_1:} multinomial p's.
#'
#' @return The output is a list of three component vectors of eigenvalues:
#'   \code{$eig.W} as Cramer-von Mises, \code{$eig.U} as Watson, and
#'   \code{$eig.A} as Anderson-Darling respectively.
#'
#' @param digits A significant digits vector.
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @seealso \code{\link{CVMStats}}
CVMEigen <- function(digits)
{
	k <- length(digits)
	ps <- Ben.ps(digits)
	# e <- n*ps
	# Tt <- cumsum(e)
	# H <- Tt/n
	H <- cumsum(ps)
	tj <- .5*(ps + c(ps[-1],ps[1]))
	Amat <- matrix(0,nrow=k,ncol=k)
	Amat[row(Amat)>=col(Amat)] <- 1
	Emat <- diag(tj)
	Dmat <- diag(ps)
	Kmat <- diag(Div(1,H*(1-H)))
	Kmat[k,k] <- 0
	Sigma.0 <- (Dmat-(ps%o%ps))
	Sigma.y <- Amat%*%Sigma.0%*%t(Amat)
      Sigma.y[,k] <- 0
      Sigma.y[k,] <- 0
      B <- mgcv::mroot(Sigma.y,rank=k-1,method='svd')
  	eig.W <- eigen(t(B)%*%Emat%*%B)$values
	Utem <- diag(1,k)-(Emat%*%(rep(1,k)%o%rep(1,k)))
	eig.U <- eigen(t(B)%*%Utem%*%Emat%*%t(Utem)%*%B)$values
	eig.A <- eigen(t(B)%*%Emat%*%Kmat%*%B)$values
	return(list(eig.W=eig.W, eig.U=eig.U, eig.A=eig.A))
}


#' Imhof numerical method
#'
#' \code{Imhof} computes upper tail probabilities of asymptotic distributions
#' for the CVM type statistics using Imhof's numerical method.
#'
#' @param x MLDEFINITION
#' @param eigvals MLDEFINITION
#' @param UPPER MLDEFINITION. If unspecified, default is NA.
#' @param subdiv MLDEFINITION. If unspecified, default is 100.
#' @param eps1 MLDEFINITION. If unspecified, default is 0.0001.
#'
#' @seealso \code{\link{Chi.approx}}
#'
#' @references Imhof, J. P. (1961). Computing the Distribution of Quadratic Forms in Normal Variables.
#' Biometrika 48, 419-426.  Imhof, J. P. (1962). Corrigenda: Computing the Distribution of
#' Quadratic Forms in Normal Variables. Biometrika 49, p. 284.
#' Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @export
Imhof <- function(x, eigvals, UPPER=NA, subdiv=100, eps1=.0001)
{
	klam <- length(eigvals)*.5
	if (is.na(UPPER))
	{
		UB <- exp(-(log(eps1*pi*klam)+.5*sum(log(abs(eigvals))))/klam)
	} else {
		UB <- UPPER
	}
   	res <- integrate(dImhof, lower=0, upper=UB, xcrit=x, lambda=eigvals, subdivisions=subdiv,
 	stop.on.error=FALSE)
 	if(res$message!="OK") print(res$message)
   	return(.5+res$value/pi)
}


#' Integrand for the \code{Imhof} function
#'
#' \code{dImhof} contains the integrand used in the Imhof numerical method function.
#'
#' @param u MLDEFINITION evaluate the integrand at u
#' @param xcrit MLDEFINITION x critical P (Q>x)
#' @param lambda MLDEFINITION eigenvalues of A\%*\%Sigma
#'
#' @seealso \code{\link{Imhof}}
#'
#' @references Imhof, J. P. (1961). Computing the Distribution of Quadratic Forms in Normal Variables.
#' Biometrika 48, 419-426.  Imhof, J. P. (1962). Corrigenda: Computing the Distribution of
#' Quadratic Forms in Normal Variables. Biometrika 49, p. 284.
#' Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
dImhof <- function(u, xcrit, lambda)
{
	ulambda <- u%o%lambda
	theta <- .5*as.vector(apply(atan(ulambda), 1, sum))- .5*xcrit*u
	rho <- exp(.25*as.vector(apply(log(1+ulambda^2), 1, sum)))
	return(sin(theta)/u/rho)
}


#' Pearson's three-moment Chi-square approximation
#'
#' \code{Chi.approx} computes upper tail probabilities of asymptotic
#' distributions for the CVM type statistics using the Pearson chi-square
#' approximation.
#'
#' The Chi-square approximation only requires the first three cumulants of the
#' statistic of interest, and is therefore faster than using the Imhof method.
#' See Lesperance et. al for further details.
#'
#' @param x MLDEFINITION
#' @param eigvals MLDEFINITION
#'
#' @seealso \code{\link{Imhof}}
#'
#' @export
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
Chi.approx <- function(x,eigvals)
{
	k1 <- sum(eigvals)
	k2 <- 2*sum(eigvals^2)
	k3 <- 8*sum(eigvals^3)
	b <- k3/4/k2
	p <- 8*k2^3/k3^2
	a <- k1-b*p
	return(1-pchisq((x-a)/b,p))
}


#' Compare Simultaneous Confidence Intervals against Benford
#'
#' \code{SimultBenp} indicates whether the simultaneous confidence intervals
#' (generated by the user's significant digit frequencies in
#' \code{\link{SimultConf}}) cover the Benford probabilities.
#'
#' @return The output is a vector containing 0's and 1's, where each 0 or 1
#'   corresponds to one of the seven simultaneous confidence intervals
#'   (Quesenberry and Hurst, Goodman, etc.). A zero indicates that at least
#'   Benford probability falls outside a given CI. A one indicates that all
#'   Benford probabilities fall within a given CI.
#'
#' @param freq Vector of multinomial frequencies.
#' @param alpha Choose \code{alpha} to generate \eqn{(1-\alpha)}% CIs. If
#'   unspecified, default is 0.05.
#' @param digits A significant digits vector which must match the frequecies \code{freq}.
#' If unspecified, default is 1:9.
#'
#' #' @examples
#' set.seed(123)
#' SimultConf((firstdigitsfreq(rnorm(100), 1)), alpha=.05)
#' SimultBenp((firstdigitsfreq(rnorm(100), 1)), alpha=.05)
#'
#' @seealso \code{\link{SimultConf}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @export
SimultBenp <- function(freq, alpha=0.05, digits=1:9)
{
	n <- sum(freq)
	ps <- matrix(Ben.ps(digits), nrow=length(digits), ncol=7)
	CI <- SimultConf(freq, alpha=alpha)
	Ind <- (CI$Lower<= ps) & (CI$Upper >= ps)
	return(apply(Ind, 2, prod))
}
