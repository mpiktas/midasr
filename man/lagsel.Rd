\name{lagsel}
\alias{lagsel}
\title{Select lags of MIDAS regression using information criteria}
\usage{
  lagsel(x, ldata = NULL, hdata = NULL, start, kmin = NULL,
    kmax = NULL, IC = list(AIC, BIC), Ofunction = "optim",
    user.gradient = FALSE, ...)
}
\arguments{
  \item{x}{the formula for MIDAS regression, the lag
  selection is performed for the first mls term in the
  formula}

  \item{ldata}{low frequency data, a \code{data.frame}
  object}

  \item{hdata}{high frequency data, a \code{data.frame}
  object}

  \item{start}{the starting values for optimisation}

  \item{kmin}{the minimum high frequency lag, defaults to
  zero.}

  \item{kmax}{the highest high frequency lag, defaults to
  square root of number of low fequency observations}

  \item{IC}{the information criterias on which to base
  selection}

  \item{Ofunction}{see \link{midasr}}

  \item{user.gradient}{see \link{midas_r}}

  \item{...}{additional parameters to optimisation
  function, see \link{midas_r}}
}
\value{
  a table with
}
\description{
  Select lags of MIDAS regression using information
  criteria
}
\author{
  Vaidotas Zemlys
}

