\name{midas.r}
\alias{midas.r}
\title{Restricted MIDAS regression}
\usage{
  midas.r(y, x, exo = NULL, k = NULL, resfun,
    start = list(resfun = NULL, exo = NULL),
    method = "BFGS", control.optim = list(), ...)
}
\arguments{
  \item{y}{the response variable, \code{\link{ts}} object.}

  \item{x}{the predictor variable, \code{\link{ts}} object,
  such that \code{frequency(x) / frequency(y)} is an
  integer.}

  \item{exo}{exogenous variables, \code{\link{ts}} object,
  having the same properties as \code{y}. The default is
  NULL, meaning that no exogenous variables are used in the
  regressidon.}

  \item{k}{the number of lags to include in MIDAS
  regression. The default is \code{NULL} meaning that the
  number of lags is calculated depending on output of
  \code{resfun}.}

  \item{resfun}{the function which returns restricted
  parameters given the restriction function. The parameters
  for the restriction function must be the supplied as
  numeric vector in the first argument of the function.
  Number of lags of the regression is calculated from the
  output of this function.}

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

  \code{coefficients} - the estimates of restricted
  coefficients

  \code{parameters} - the estimates of parameters of the
  restriction function \code{data} - output of
  \code{\link{mmatrix.midas}}, the data matrix used for
  fitting.

  \code{opt} - the output of call to \code{\link{optim}}

  \code{exo.coef} - the estimates of the coefficients of
  exogenous variables

  \code{restr.fun} - the restriction function used in
  estimation
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

##Fit restricted model
mr <- midas.r(y,x,resfun=theta.h0,start=list(resfun=c(-0.1,10,-10,-10)),dk=4*12)

##Include intercept and trend in regression

mr.it <- midas.r(y,x,exo=cbind(1,1:500),resfun=theta.h0,start=list(resfun=c(-0.1,10,-10,-10),exo=c(0,0)),dk=4*12)
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

