\name{simplearma.sim}
\alias{simplearma.sim}
\title{Simulate AR(1) or MA(1) model}
\usage{
simplearma.sim(model, n, innov.sd, frequency, n.start = 300)
}
\arguments{
  \item{model}{A named vector of length one. Name is either
  "ar", or "ma" depending on which AR(1) or MA(1) process
  should be generated}

  \item{n}{the length of output series}

  \item{innov.sd}{the standard error of innovations, which
  are zero mean normal random variables}

  \item{frequency}{the frequency of the regressor, should
  be larger than one.}

  \item{n.start}{the length of the burn.in period, the
  default is 300.}
}
\value{
a time-series object of class \code{ts}
}
\description{
Simulate MIDAS regressor as a AR(1) or MA(1) time series
}
\examples{
#Generate AR(1) model with rho=0.6, with frequency 12
x <- simplearma.sim(list(ar=0.6),1500*12,1,12)
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}

