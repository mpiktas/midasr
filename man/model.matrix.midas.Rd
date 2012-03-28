\name{model.matrix.midas}
\alias{model.matrix.midas}
\title{Arange data for fitting}
\usage{
  model.matrix.midas(y, x, k = 0)
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
  Given response and predictor variables arrange them into
  form for easy fitting.
}
\author{
  Vaidotas Zemlys (zemlys
}

