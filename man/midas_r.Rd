\name{midas_r}
\alias{midas_r}
\alias{midas_r.formula}
\title{Restricted MIDAS regression}
\usage{
  midas_r(x, ...)

  \method{midas_r}{formula} (formula, ldata, hdata, start,
    optim = list(func = "optim", method = "BFGS"), ...)
}
\arguments{
  \item{x}{either formula for restricted MIDAS regression
  or \code{midas__r} object.}

  \item{formula}{formula for restricted MIDAS regression}

  \item{ldata}{low frequency data, a \code{data.frame}
  object}

  \item{hdata}{high frequency data, a \code{data.frame}
  object}

  \item{start}{the starting values for optimisation. Must
  be a list with named elements.}

  \item{optim}{the list containing the name of optimisation
  function and its arguments. The default is to use
  \code{\link{optim}} with \code{method="BFGS"}}

  \item{...}{additional arguments supplied for
  \code{resfun}}
}
\value{
  a \code{midas_r} object which is the list with the
  following elements:

  \item{coefficients}{the estimates of parameters of
  restrictions} \item{midas.coefficientas}{the estimates of
  restricted coefficients of MIDAS regression}
  \item{model}{model data} \item{restrictions}{the
  restriction function(s) used in estimation.}
  \item{unrestricted}{unrestricted regression estimated
  using \code{\link{midas_u}}} \item{param.map}{parameter
  map for optimisation function} \item{fn0}{optimisation
  function for non-linear least squares problem solved in
  restricted MIDAS regression} \item{rhs}{the function
  which evaluates the right-hand side of the MIDAS
  regression} \item{allcoef}{the function which evaluates
  the restricted coefficientsof MIDAS regression}
  \item{opt}{the output of optimisation procedure}
  \item{argmap.opt}{the list containing the name of
  optimisation function together with arguments for
  optimisation function} \item{start.opt}{the starting
  values used in optimisation}
}
\description{
  Estimate restricted MIDAS regression using non-linear
  least squares. Currently only the estimates of the
  parameters are given, without their standard errors.
}
\details{
  Given MIDAS regression:

  \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i}
  x_{(t-j)m-i}+\mathbf{z_t}\beta+u_t}

  estimate the parameters of the restriction

  \deqn{\theta_h=g(h,\lambda),} where \eqn{h=0,...,(k+1)m},
  together with coefficients \eqn{\beta} corresponding to
  additional low frequency regressors.

  MIDAS regression involves times series with different
  frequencies.

  The restriction function must return the restricted
  coefficients of the MIDAS regression.
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
x <- window(x,start=start(y))

##Fit restricted model
mr <- midas_r(y~mdslag(x,3,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)

##Include intercept and trend in regression

mr.it <- midas_r(y~mdslag(x,3,theta.h0)+trend,data.frame(y=y,trend=1:500),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}
\seealso{
  midas_r.midas_r
}

