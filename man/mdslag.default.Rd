\name{mdslag.default}
\alias{mdslag.default}
\title{Lag a mixed data sampling time series}
\usage{
  \method{mdslag}{default} (x, k, m, ...)
}
\arguments{
  \item{x}{a vector for which mixed data sampling lag has
  to be computed}

  \item{k}{the number of lags of low frequency. This can
  either be a number indicating the highest lag, or the
  vector with lags which should be included.}

  \item{m}{frequency ratio}

  \item{...}{further arguments}
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
\seealso{
  mdslag mdslag.mds
}

