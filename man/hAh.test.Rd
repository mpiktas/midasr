\name{hAh.test}
\alias{hAh.test}
\title{Test restrictions on coefficients of MIDAS regression}
\usage{
  hAh.test(restricted, gr = NULL, ...)
}
\arguments{
  \item{restricted}{the restricted model, estimated with
  \code{\link{midas.r}}}

  \item{gr}{the gradient of the restriction function. Must
  return the matrix with dimensions \eqn{d_k \times q},
  where \eqn{d_k} and \eqn{q} are the numbers of
  coefficients in unrestricted and restricted regressions
  correspondingly. Default value is \code{NULL}, which
  means that the numeric approximation of gradient is
  calculated.}

  \item{...}{the parameters supplied to gradient function,
  if \code{gr} is not \code{NULL}}
}
\value{
  a \code{htest} object
}
\description{
  Perform a test whether the restriction on MIDAS
  regression coefficients holds.
}
\details{
  Given MIDAS regression:

  \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i}
  x_{(t-j)m-i}+u_t}

  test the null hypothesis that the following restriction
  holds:

  \deqn{\theta_h=g(h,\lambda),} where \eqn{h=0,...,(k+1)m}.
}
\examples{
##The parameter function
theta.h0 <- function(p, dk) {
   i <- (1:dk-1)
   (p[1] + p[2]*i)*exp(p[3]*i + p[4]*i^2)
}

##Generate coefficients
theta0 <- theta.h0(c(-0.1,0.1,-0.1,-0.001),4*12)

##Plot the coefficients
plot(theta0)

##Generate the predictor variable
set.seed(13)
x <- simplearma.sim(list(ar=0.6),1500*12,1,12)

##Simulate the response variable
y <- midas.sim(500,theta0,x,1)

##Fit restricted model
mr <- midas.r(y,x,resfun=theta.h0,start=list(resfun=c(-0.1,0.1,-0.1,-0.001)),dk=4*12)
mu <- midas.u(y,x,k=3)

##The gradient function
grad.h0<-function(p, dk) {
   i <- (1:dk-1)
   a <- exp(p[3]*i + p[4]*i^2)
   cbind(a, a*i, a*i*(p[1]+p[2]*i), a*i^2*(p[1]+p[2]*i))
}

##Perform test (the expected result should be the acceptance of null)

hAh.test(mr,grad.h0,dk=4*12)

##Use numerical gradient instead of supplied one
hAh.test(mr)
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}
\references{
  Kvedaras V., Zemlys, V. \emph{Testing the functional
  constraints on parameters in regressions with variables
  of different frequency} Economics Letters 116 (2012)
  250-254
}

