\name{mmatrix.midas}
\alias{mmatrix.midas}
\title{Combine one low and one high frequency time series objects}
\usage{
  mmatrix.midas(y, x, k = 0)
}
\arguments{
  \item{y}{the response variable, a \code{ts} object}

  \item{x}{the predictor variable, a \code{ts} object}

  \item{k}{the number of lags of predictor variable to
  include}
}
\value{
  a matrix
}
\description{
  Given low and high frequency time series, combine them
  into matrix suitable for MIDAS regresion estimation.
}
\author{
  Vaidotas Zemlys
}

