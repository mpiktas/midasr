\name{midas.u}
\alias{midas.u}
\title{Unrestricted MIDAS regression}
\usage{
  midas.u(y, x, k)
}
\arguments{
  \item{y}{the response variable}

  \item{x}{the predictor variable}

  \item{k}{the number of lags to include in MIDAS
  regression}
}
\value{
  \code{lm} object
}
\description{
  Estimate unrestricted MIDAS regression using OLS. This
  function is basically a wrapper for \code{lm}.
}
\details{
  MIDAS regression has the following form:

  \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i}
  x_{(t-j)m-i}+u_t}

  or alternatively

  \deqn{y_t=\sum_{h=0}^{(k+1)m}\theta_hx_{tm-h}+u_t,} where
  \eqn{m} is the frequency of high-frequency data \eqn{x}
  and \eqn{k} is the number of lags included in the
  regression.

  Given certain assumptions the coefficients can be
  estimated using usual OLS and they have the familiar
  properties associated with simple linear regression.
}
\examples{
##The parameter function
theta.h0 <- function(p, dk) {
   i <- (1:dk-1)/100
   pol <- p[3]*i + p[4]*i^2
   (p[1] + p[2]*i)*exp(pol)
}

##Generate coefficients
theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)

##Plot the coefficients
plot(theta0)

##Generate the predictor variable
x <- simplearma.sim(list(ar=0.6),1500*12,1,12)

##Simulate the response variable
y <- midas.sim(500,theta0,x,1)

##Fit unrestricted model
midas.u(y,x,3)
}
\author{
  Vaidotas Zemlys (zemlys
}
\references{
  Kvedaras V., Zemlys, V. \emph{Testing the functional
  constraints on parameters in regressions with variables
  of different frequency}
  \url{http://dx.doi.org/10.1016/j.econlet.2012.03.009}
}

