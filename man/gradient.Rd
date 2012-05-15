\name{gradient}
\alias{gradient}
\alias{gradient.midas_r}
\title{Calculate numerical approximation of gradient of RSS of
MIDAS regression}
\usage{
  gradient(x, ...)

  \method{gradient}{midas_r} (x, ...)
}
\arguments{
  \item{x}{a \code{\link{midas_r}} object}

  \item{...}{further arguments to function
  \code{\link{grad}} from package \code{numDeriv}}
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

