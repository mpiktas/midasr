\name{harstep}
\alias{harstep}
\title{HAR(3)-RV model MIDAS weights specification}
\usage{
harstep(p, d, m)
}
\arguments{
  \item{p}{parameters for Almon lag}

  \item{d}{number of the coefficients}

  \item{m}{the frequency, currently ignored.}
}
\value{
vector of coefficients
}
\description{
HAR(3)-RV model MIDAS weights specification
}
\details{
MIDAS weights for Heterogeneous Autoregressive model of
Realized Volatilty (HAR-RV). It is assumed that month has
20 days.
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}
\references{
Corsi, F., \emph{A Simple Approximate Long-Memory Model of
Realized Volatility}, Journal of Financial Econometrics
Vol. 7 No. 2 (2009) 174-196
}

