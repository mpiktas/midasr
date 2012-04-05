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
  \code{ts} object, vector or matrix. It is forcibly
  converted to \code{ts} object with exact properties of
  \code{y}.}

  \item{k}{the number of lags of predictor variable to
  include. This can either be a number indicating the
  highest lag, or the vector with lags which should be
  included. Note that the here the lag is low frequency
  lag. Hence effectively high frequency variable is lagged
  \eqn{k\times m} times, where \eqn{m} is frequency ratio,
  which is calculated as \code{frequency(x) %/%
  frequency(y)}}
}
\value{
  a matrix
}
\description{
  Given low and high frequency time series, combine them
  into matrix suitable for MIDAS regresion estimation.
}
\details{
  Note that this function does not check whether data is
  conformable. If \code{exo} variable does not have the
  same \code{start}, \code{end} and \code{frequency} as
  response variable \code{y}, error is produced.
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

