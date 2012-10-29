\name{embedlf.mds}
\alias{embedlf.mds}
\title{Lag a mixed data sampling time series}
\usage{
  \method{embedlf}{mds} (x, k, ...)
}
\arguments{
  \item{x}{an object of class \code{\link{mds}}}

  \item{k}{the number of lags of low frequency. This can
  either be a number indicating the highest lag, or the
  vector with lags which should be included.}

  \item{...}{further arguments}
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
\seealso{
  embedlf embedlf.default
}

