\name{mmatrix.midas}
\alias{mmatrix.midas}
\title{Combine one low and one high frequency time series objects}
\usage{
  mmatrix.midas(y, x, exo = NULL, k = 0)
}
\arguments{
  \item{y}{the response low frequency variable, a \code{ts}
  object}

  \item{x}{the predictor high frequency variable, a
  \code{ts} object}

  \item{exo}{the exogenous low frequency variables, a
  \code{ts} object.}

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
  Virmantas Kvedaras, Vaidotas Zemlys
}

