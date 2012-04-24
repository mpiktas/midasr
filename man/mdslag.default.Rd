\name{mdslag.default}
\alias{mdslag}
\alias{mdslag.default}
\alias{mdslag.mds}
\title{Lag a mixed data sampling time series}
\usage{
  mdslag.default(x, k, m, ...)
}
\arguments{
  \item{x}{a vector for which mixed data sampling lag has
  to be computed}

  \item{k}{the number of lags of low frequency. This can
  either be a number indicating the highest lag, or the
  vector with lags which should be included.}

  \item{m}{frequency ratio}

  \item{...}{}
}
\value{
  a matrix containing the lags
}
\description{
  Compute a lagged version of mixed data sampling time
  series
}
\author{
  Vaidotas Zemlys
}

