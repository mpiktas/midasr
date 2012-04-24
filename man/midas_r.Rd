\name{midas_r}
\alias{midas_r}
\title{Restricted MIDAS regression}
\usage{
  midas_r(formula, ldata, hdata, start, method = "BFGS",
    control.optim = list(), ...)
}
\arguments{
  \item{formula}{}

  \item{ldata}{}

  \item{hdata}{}

  \item{start}{the starting values for optimisation. Must
  be a list with named elements \code{resfun} containing
  vector of starting values for restriction function and
  \code{exo} containing vector of starting values for
  exogenous variables.}

  \item{method}{the method used for optimisation, see
  \code{\link{optim}} documentation. All methods are
  suported except "L-BFGS-B" and "Brent". Default method is
  "BFGS".}

  \item{control.optim}{a list of control parameters for
  \code{\link{optim}}.}

  \item{...}{additional arguments supplied for
  \code{resfun}}
}
\value{
  a list with the following elements:

  \item{coefficients}{the estimates of parameters of
  restrictions} \item{midas.coefficientas}{the estimates of
  restricted coefficients of MIDAS regression}
  \item{model}{model data} \item{opt}{the output of call to
  \code{\link{optim}}} \item{restr.fun}{the restriction
  function used in estimation.}
  \item{unrestricted}{unrestricted regression estimated
  using \code{\link{midas.u}} \item{restr.no}{the number of
  parameters used in restriction function} }
}
\description{
  Estimate restricted MIDAS regression using non-linear
  least squares. Uses \code{optim} for optimisation.
  Currently only the estimates of the parameters are given,
  without their standard errors.
}
\details{
  Given MIDAS regression:

  \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i}
  x_{(t-j)m-i}+u_t}

  estimate the parameters of the restriction

  \deqn{\theta_h=g(h,\lambda),} where \eqn{h=0,...,(k+1)m}.

  MIDAS regression involves times series with different
  frequencies. In R the frequency property is set when
  creating time series objects \code{\link{ts}}. Hence the
  frequency ratio \eqn{m} which figures in MIDAS regression
  is calculated from frequency property of time series
  objects supplied.

  The restriction function must return the restricted
  coefficients of the MIDAS regression. If specific lag
  structure is needed (see more in
  \code{\link{mmatrix.midas}}) make sure that the number of
  coefficients returned by restriction function coincides
  with number of columns of the lagged high frequency
  predictor variable matrix.
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

##Remove unnecessary history of x
x <- window(x,start=c(1500-500+1,1))

##Fit restricted model
mr <- midas_r(y~mdslag(x,4,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)

##Include intercept and trend in regression

mr.it <- midas_r(y~mdslag(x,4,theta0)+trend,data.frame(y=y,trend=1:500),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

