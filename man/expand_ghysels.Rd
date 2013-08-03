\name{expand_ghysels}
\alias{expand_ghysels}
\title{Create table of weights, lags and starting values for Ghysels weight schema}
\usage{
  expand_ghysels(weight, type = c("A", "B", "C"), kmin,
    kmax, m, start, mkmin = 0)
}
\arguments{
  \item{weight}{the names of weight functions}

  \item{type}{the type of Ghysels schema, \code{"A"},
  \code{"B"} or \code{"C"}}

  \item{kmin}{the minimum low frequency lag}

  \item{kmax}{the maximum high frequency lag}

  \item{m}{the frequency ratio}

  \item{start}{the starting values for the weights of the
  one low frequency lag}

  \item{mkmin}{the minimum low frequency lag from which to
  start the fitting}
}
\value{
  a \code{lws_table} object, a list with elements
  \code{weights}, \code{lags} and \code{starts}
}
\description{
  Create table of weights, lags and starting values for
  Ghysels weight schema, see \link{ghyselslag}
}
\details{
  Given weight function creates lags starting from
  \code{kmin} to \code{kmax} and replicates starting values
  for each low frequency lag.
}
\examples{
expand_ghysels("nealmon","A",1,2,12,c(0,0,0))
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

