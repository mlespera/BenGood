#' Generate Simultaneous Confidence Intervals
#'
#' \code{SimultConf} takes a vector of observed cell frequencies from a
#' multinomial distribution, \code{N}, and given significance level \eqn{\alpha},
#' returns simultaneous confidence intervals for proportions
#' using seven different methods:
#' \enumerate{
#'     \item Quesenberry and Hurst
#'     \item Goodman
#'     \item Bailey angular transformation
#'     \item Bailey square root transformation
#'     \item Fitzpatrick and Scott
#'     \item Sison and Glaz
#'     \item Univariate approximate Binomial confidence intervals.
#' }
#'
#' An easy way to check whether your simultaneous confidence intervals generated
#' by \code{SimultConf} cover the Benford probabilities is to use
#' \code{\link{SimultBenp}}.
#'
#' See Lesperance et. al (2016) for the formulae of each method.
#'
#' @return The output is two matrices: a Lower matrix containing the lower
#'   bounds, and an Upper matrix containing the upper bounds of the
#'   simultaneous confidence intervals using the methods specified by the user.
#'   The univariate approximate binomial CIs are always included.
#'
#' @param N Vector of observed frequencies.
#' @param method The methods for simultaneous confidence intervals for
#'  multinomial proportions. This parameter can take on the following values, or
#'  a vector containing any combination of these values separated using commas:
#'
#'  \enumerate{
#'       \item 'Quesenberry'
#'       \item 'Goodman'
#'       \item 'BaileyAng'
#'       \item 'Baileysqrt'
#'       \item 'Fitzpatrick'
#'       \item 'Sison'
#'       \item 'binomial'}
#'
#'  Specifying no method will return all confidence intervals.
#' @param alpha Choose \code{alpha} to generate \eqn{(1-\alpha)}\% CIs. If
#'   unspecified, default is 0.05.
#'
#' @seealso \code{\link{SimultBenp}}
#'
#' @examples
#' myCI = SimultConf(1:9, c('Good', 'Quesenberry'))
#' lower = myCI$Lower
#'
#' @export
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
SimultConf <- function(N, method=c('Quesenberry','Goodman','BaileyAng','Baileysqrt','Fitzpatrick','Sison',
                             'binomial'), alpha=.05)
{
   methods = c('Quesenberry','Goodman','BaileyAng','Baileysqrt','Fitzpatrick','Sison',
              'binomial')
	 meth = match.arg(method, methods, several.ok=TRUE)

	n <- sum(N)
	k <- length(N)
	A <- qchisq(1-alpha, k-1)
	B <- qchisq(1-(alpha/k), 1)
	C <- B/(4*n)
	D <- qnorm(1-alpha/4)
	E <- qnorm(1-alpha/2)
	phat <- N/n

	### Quesenberry and Hurst's method ###
	Lower1 <- pmax(0, round((A + 2*N - sqrt(A*(A+4*(N*(n-N)/n))))/(2*(n+A)), 4))
	Upper1 <- pmin(1, round((A + 2*N + sqrt(A*(A+4*(N*(n-N)/n))))/(2*(n+A)), 4))

	### Goodman's method ###
	Lower2 <- pmax(0, round((B + 2*N - sqrt(B*(B+4*(N*(n-N)/n))))/(2*(n+B)), 4))
	Upper2 <- pmin(1, round((B + 2*N + sqrt(B*(B+4*(N*(n-N)/n))))/(2*(n+B)), 4))

	### Bailey's angular transformation method ###
	Lower3 <- pmax(0, round((sin(asin(sqrt((N+3/8)/(n+3/4))) - sqrt(B/(4*n+2))))^2, 4))
	Upper3 <- pmin(1, round((sin(asin(sqrt((N+3/8)/(n+3/4))) + sqrt(B/(4*n+2))))^2, 4))

	### Bailey's square root transformation method ###
	Lower4 <- pmax(0, round((sqrt((N+3/8)/(n+1/8)) - sqrt(C*(C+1-(N+3/8)/(n+1/8))))^2 / (C+1)^2, 4))
	Upper4 <- pmin(1, round((sqrt((N+3/8)/(n+1/8)) + sqrt(C*(C+1-(N+3/8)/(n+1/8))))^2 / (C+1)^2, 4))

	### Fiztpatrick and Scott's method ###
	Lower5 <- pmax(0, round((phat) - D/(2*sqrt(n)), 4))
	Upper5 <- pmin(1, round((phat) + D/(2*sqrt(n)), 4))

	### Sison and Glaz's method ###
	tau <- 1
	vtau <- v(N, tau)
	vtaupl <- v(N, tau+1)

	while(vtaupl < (1-alpha)){
		tau <- tau+1
		vtau <- vtaupl
		vtaupl <- v(N, tau+1)
	}
	gamma <- ((1-alpha)-vtau)/(vtaupl-vtau)

	Lower6 <- pmax(0, round((phat)-(tau/n), 4))
	Upper6 <- pmin(1, round((phat)+(tau+2*gamma)/n, 4))

	### Simple univariate Normal approx to Binomial intervals ###
	Lower7 <- pmax(0, round(phat - E*sqrt(phat*(1-phat)/n), 4))
	Upper7 <- pmin(1, round(phat + E*sqrt(phat*(1-phat)/n), 4))

	tnames <- c("Ques", "Good", "Bang", "Bsqrt", "Fitz", "Sison", "UniBin")
  tnames = paste(tnames,1-alpha)
	Lower <- cbind(Lower1, Lower2, Lower3, Lower4, Lower5, Lower6, Lower7)
	Upper <- cbind(Upper1, Upper2, Upper3, Upper4, Upper5, Upper6, Upper7)
	dimnames(Lower)[[2]] <- tnames
	dimnames(Upper)[[2]] <- tnames
	m = is.element(methods,meth)
	return(list(Lower=Lower[,m, drop=FALSE], Upper=Upper[,m, drop=FALSE]))
}


#' Function for Sison and Glaz's method
#'
#' \code{v} is a function required by SimultConf to compute Sison and Glaz's
#' method.
#'
#' @param N MLDEFINITION
#' @param tau MLDEFINITION An integer.
#'
#' @seealso \code{\link{SimultConf}}
#' @seealso \code{\link{f}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @keywords internal
	v <- function(N,tau)
{
	n <- sum(N)
	k <- length(N)
	lambda <- N
	a <- N+tau
	b <- N-tau

	mu1 <- lambda^1*(1+(((ppois(b-1,N)-ppois(b-2,N))-(ppois(a,N)-ppois(a-1,N)))/(ppois(a,N)-ppois(b-1,N))))
	mu2 <- lambda^2*(1+(((ppois(b-1,N)-ppois(b-3,N))-(ppois(a,N)-ppois(a-2,N)))/(ppois(a,N)-ppois(b-1,N))))
	mu3 <- lambda^3*(1+(((ppois(b-1,N)-ppois(b-4,N))-(ppois(a,N)-ppois(a-3,N)))/(ppois(a,N)-ppois(b-1,N))))
	mu4 <- lambda^4*(1+(((ppois(b-1,N)-ppois(b-5,N))-(ppois(a,N)-ppois(a-4,N)))/(ppois(a,N)-ppois(b-1,N))))

	sigma <- sqrt(mu2+mu1-mu1^2)
	mu3i <- mu3+3*mu2-3*mu2*mu1+mu1-3*mu1^2+2*mu1^3
	mu4i <- mu4+6*mu3+7*mu2+mu1-4*mu1*mu3-12*mu1*mu2-4*mu1^2+6*mu1^2*mu2+6*mu1^3-3*mu1^4

	gamma1 <- sum(mu3i)/sum(sigma^2)^(3/2)
	gamma2 <- sum(mu4i-3*sigma^4)/(sum(sigma^2)^2)

	A <- 1/(dpois(n,n))
	S <- (n-sum(mu1))/(sqrt(sum(sigma^2)))
	vtau <- A*prod(ppois(a,N)-ppois(b-1,N))*f(S,gamma1,gamma2)/sqrt(sum(sigma^2))
	return(vtau)
}

#' Function for Sison and Glaz's \code{v}
#'
#' \code{f} is a function required by \code{v}.
#'
#' @param x MLDEFINITION
#' @param gamma1 MLDEFINITION
#' @param gamma2 MLDEFINITION
#'
#' @seealso \code{\link{SimultConf}}
#' @seealso \code{\link{f}}
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#'   (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#'   simultaneous confidence intervals. PLoS one; 11(3).
#'   Wong, S. (2010)  Testing Benford's Law with the first two significant
#'   digits.  University of Victoria, Master's thesis.
#'
#' @keywords internal
f <- function(x, gamma1, gamma2)
{
	return((exp(-x^2/2)/sqrt(2*pi)) * (1+(gamma1/6 *(x^3-3*x)) +
		( gamma2/24*(x^4-6*x^2+3)) + ( gamma1^2/72*(x^6-15*x^4+45*x^2-15))))
}


#' Plot Simultaneous Confidence Intervals
#'
#' \code{plotCI} plots the simultaneous confidence intervals from the Lower and
#' Upper bound matrices generated by \code{\link{SimultConf}}.
#'
#' The Lower and Upper matrices are matrices of lower and upper limits, where
#' the columns represent the different interval types. Both
#' matrices must be of the same size.
#'
#' @param Lower Matrix of lower interval limits.
#' @param Upper Matrix of upper interval limits.
#' @param digits A significant digits vector. If unspecified, default is 1:9.
#' @param benps \code{TRUE} will plot Benford probabilities on the graphs produced.
#' @param legloc \code{TRUE} will allow user to click to locate the legend.
#' If unspecified, default is FALSE. Note that this feature is only supported on
#' machines where \code{locator} is supported.
#'
#' @seealso \code{\link{SimultConf}}
#' @seealso \code{\link{SimultBenp}}
#' @seealso \code{\link{Ben.ps}}
#'
#' @export
#'
#' @examples
#' myN = c(28, 19, 18,  7,  9,  6,  3,  4,  6)
#' myCI = SimultConf(myN, c('Good','Fitz'))
#' plotCI(myCI$Lower, myCI$Upper)
#'
#' @references Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
#' (2016) Assessing conformance with Benford's Law: goodness-of-fit tests and
#' simultaneous confidence intervals. PLoS one; 11(3).
#' Wong, S. (2010)  Testing Benford's Law with the first two significant
#' digits.  University of Victoria, Master's thesis.
plotCI=function(Lower, Upper, digits=1:9, benps=TRUE, legloc=FALSE){
  Ben <- Ben.ps(digits)
  methods = dimnames(Lower)[[2]]
  for (i in 1:dim(Lower)[2]){
    plot(digits, Ben, ylim=c(0,1), xlab="digits", ylab="proportions", cex=1.25)
    title(paste(methods[i],' Simultaneous CI'))
    points(digits, Lower[,i], pch="-", cex=2)
    points(digits, Upper[,i], pch="-", cex=2)
    x <- rep(digits,each=3)
    y <- as.vector(t(cbind(Lower[,i],Upper[,i],NA)))
    lines(x,y)
    if(legloc){  cat('click on plot for legend location \n')
      legend(locator(1), legend="Benford's p", pch=1)
    }
  }
}
