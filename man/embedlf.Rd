\name{embedlf}
\alias{embedlf}
\alias{embedlf.default}
\title{Lag a mixed data sampling time series}
\usage{
  embedlf(x, ...)

  \method{embedlf}{default} (x, k, m, ...)
}
\arguments{
  \item{x}{a vector}

  \item{k}{the number of high frequency lags including
  contemporaneous lag}

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
  embedlf.default
}

