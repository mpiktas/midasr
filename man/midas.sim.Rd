\name{midas.sim}
\alias{midas.sim}
\title{Simulate MIDAS regresion response variable}
\usage{
  midas.sim(n, theta, x, eps.sd)
}
\arguments{
  \item{n}{The sample size}

  \item{theta}{a vector with MIDAS regression coefficients}

  \item{x}{a \code{ts} object with MIDAS regression
  predictor variable}

  \item{eps.sd}{the standard error of the regression
  disturbances, which are assumed to be independent normal
  zero mean random variables}
}
\value{
  a \code{ts} object
}
\description{
  Given the predictor variable and the coefficients
  calculate MIDAS regression response variable. It is
  assumed that predictor variable has long enough history.
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
}
\author{
  Vaidotas Zemlys
}

