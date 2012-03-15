\name{midas.sim}
\alias{midas.sim}
\title{Simulate MIDAS regresion response variable}
\usage{
  midas.sim(n, theta, x, eps.sd)
}
\arguments{
  \item{n}{The sample size}

  \item{theta}{a vector with MIDAS regression parameters}

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
\author{
  Vaidotas Zemlys
}

