\name{hessian}
\alias{hessian}
\alias{hessian.midas_r}
\title{Calculate numerical approximation of gradient of RSS of
MIDAS regression}
\usage{
  hessian(x, ...)

  \method{hessian}{midas_r} (x, ...)
}
\arguments{
  \item{x}{a \code{\link{midas_r}} object}

  \item{...}{further arguments to function
  \code{\link{hessian}} from package \code{numDeriv}}
}
\value{
  numeric vector
}
\description{
  Calculate nummerical approximation of gradient of
  residual sum of squares function of estimated MIDAS
  regression
}
\author{
  Vaidotas Zemlys
}
\seealso{
  midas_r.formula
}

