\name{midas.auto.sim}
\alias{midas.auto.sim}
\title{Simulate autoregressive MIDAS model}
\usage{
midas.auto.sim(n, theta, alpha, x, eps.sd, n.start = NA)
}
\arguments{
  \item{n}{sample size}

  \item{theta}{a vector with MIDAS weights for predictor
  variable}

  \item{alpha}{autoregressive coefficients}

  \item{x}{a high frequency predictor variable}

  \item{eps.sd}{the standard error of the regression
  disturbances, which are assumed to be independent normal
  zero mean random variables}

  \item{n.start}{number of observations to omit for the
  burn.in}
}
\value{
a \code{ts} object
}
\description{
Given the predictor variable, the weights and
autoregressive coefficients calculate MIDAS regression
response variable.
}
\examples{
theta.h0 <- function(p, dk) {
  i <- (1:dk-1)/100
  pol <- p[3]*i + p[4]*i^2
  (p[1] + p[2]*i)*exp(pol)
}

##Generate coefficients
theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)

##Generate the predictor variable
xx <- simplearma.sim(list(ar=0.6),3000*12,1,12)

y <- midas.auto.sim(500,theta0,c(0.5),xx,1,n.start=100)
x <- window(xx,start=start(y))
midas_r(y~mls(y,1,1)+fmls(x,4*12-1,12,theta.h0),start=list(x=c(-0.1,10,-10,-10)))
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}

