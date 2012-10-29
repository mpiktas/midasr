\name{embedlf.default}
\alias{embedlf.default}
\title{Lag a mixed data sampling time series}
\usage{
  \method{embedlf}{default} (x, k, m, ...)
}
\arguments{
  \item{x}{a vector for which mixed data sampling lag has
  to be computed}

  \item{k}{the number of high frequency lags plus 1.}

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
  embedlf embedlf.mds
}

