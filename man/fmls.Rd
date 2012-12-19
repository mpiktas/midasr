\name{fmls}
\alias{fmls}
\title{Full MIDAS lag structure}
\usage{
  fmls(x, k, m, ...)
}
\arguments{
  \item{x}{a vector}

  \item{k}{maximum lag order}

  \item{m}{frequency ratio}

  \item{...}{further arguments}
}
\value{
  a matrix containing the lags
}
\description{
  Create a matrix of MIDAS lags, including contemporaneous
  lag up to selected order.
}
\details{
  This is a convenience function, it calls \code{link{msl}}
  to perform actual calculations.
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}
\seealso{
  mls
}

