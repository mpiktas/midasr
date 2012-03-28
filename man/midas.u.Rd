\name{midas.u}
\alias{midas.u}
\title{Unrestricted MIDAS regression}
\usage{
  midas.u(y, x, k)
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
  Estimate unrestricted MIDAS regression using OLS. This
  function is basically a wrapper for \code{lm}.
}
\author{
  Vaidotas Zemlys
}

