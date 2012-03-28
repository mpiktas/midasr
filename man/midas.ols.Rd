\name{midas.ols}
\alias{midas.ols}
\title{Unrestricted MIDAS regression}
\usage{
  midas.ols(y, x, k)
}
\arguments{
  \item{y}{the response variable}

  \item{x}{the predictor variable}

  \item{k}{the number of lags to include in MIDAS
  regression}
}
\value{
  \code{lm} object
}
\description{
  Calculate unrestricted MIDAS regression using OLS. This
  function is basically a wrapper for \code{lm}
}
\author{
  Vaidotas Zemlys
}

