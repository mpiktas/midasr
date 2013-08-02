\name{midas_u}
\alias{midas_u}
\title{Estimate unrestricted MIDAS regression}
\usage{
  midas_u(formula, data, ...)
}
\arguments{
  \item{formula}{MIDAS regression model formula}

  \item{data}{a named list containing data with mixed
  frequencies}

  \item{...}{further arguments, which could be passed to
  \code{\link{lm}} function.}
}
\value{
  \code{\link{lm}} object.
}
\description{
  Estimate unrestricted MIDAS regression using OLS. This
  function is a wrapper for \code{lm}.
}
\details{
  MIDAS regression has the following form:

  \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i}
  x_{(t-j)m-i}+\mathbf{z_t}\mathbf{\beta}+u_t}

  or alternatively

  \deqn{y_t=\sum_{h=0}^{(k+1)m}\theta_hx_{tm-h}+\mathbf{z_t}\mathbf{\beta}+u_t,}
  where \eqn{m} is the frequency ratio and \eqn{k} is the
  number of lags included in the regression.

  Given certain assumptions the coefficients can be
  estimated using usual OLS and they have the familiar
  properties associated with simple linear regression.

  MIDAS regression involves times series with different
  frequencies.
}
\examples{
##The parameter function
theta.h0 <- function(p, dk, ...) {
   i <- (1:dk-1)/100
   pol <- p[3]*i + p[4]*i^2
   (p[1] + p[2]*i)*exp(pol)
}

##Generate coefficients
theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)

##Plot the coefficients
##Do not run
#plot(theta0)

##Generate the predictor variable
x <- simplearma.sim(list(ar=0.6),1500*12,1,12)

##Simulate the response variable
y <- midas.sim(500,theta0,x,1)

##Create low frequency data.frame
ldt <- data.frame(y=y,trend=1:length(y))

##Create high frequency data.frame

hdt <- data.frame(x=window(x, start=start(y)))

##Fit unrestricted model
mu <- midas_u(y~fmls(x,2,12)-1, list(ldt, hdt))

##Include intercept and trend in regression

mu.it <- midas_u(y~fmls(x,2,12)+trend, list(ldt, hdt))

##Pass data as partialy named list

mu.it <- midas_u(y~fmls(x,2,12)+trend, list(ldt, x=hdt$x))
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

