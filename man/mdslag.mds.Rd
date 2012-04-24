\name{mdslag.mds}
\alias{mdslag}
\alias{mdslag.default}
\alias{mdslag.mds}
\title{Lag a mixed data sampling time series}
\usage{
  \method{mdslag}{mds} (x, k, ...)
}
\arguments{
  \item{x}{an object of class \code{\link{mds}}}

  \item{k}{the number of lags of low frequency. This can
  either be a number indicating the highest lag, or the
  vector with lags which should be included.}

  \item{...}{}
}
\value{
  a matrix containing lags
}
\description{
  Compute a lagged version of mixed data sampling time
  series
}
\author{
  Vaidotas Zemlys
}

