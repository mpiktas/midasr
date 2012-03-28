\name{hAh.test}
\alias{hAh.test}
\title{Test restrictions on coefficients of MIDAS regression}
\usage{
  hAh.test(unrestricted, restricted, gr, ...)
}
\arguments{
  \item{unrestricted}{the unrestricted model}

  \item{restricted}{the restricted model}

  \item{gr}{the gradient of the restriction function. Must
  return the matrix with dimensions \code{l x p}, where
  \code{l} is the number of unrestricted coefficients of
  the regression and \code{p} is the number of parameters
  in the restriction function.}

  \item{...}{the parameters supplied to gradient function}
}
\value{
  a \code{htest} object
}
\description{
  Perform a test whether the restriction on MIDAS
  regression coefficients holds.
}
\author{
  Vaidotas Zemlys
}

