\name{expand_amidas}
\alias{expand_amidas}
\title{Create table of weights, lags and starting values for Ghysels weight schema}
\usage{
  expand_amidas(weight, type = c("A", "B", "C"), from = 0,
    to, m, start)
}
\arguments{
  \item{weight}{the names of weight functions}

  \item{type}{the type of Ghysels schema, \code{"A"},
  \code{"B"} or \code{"C"}}

  \item{from}{the high frequency lags from which to start
  the fitting}

  \item{to}{to a vector of length two, containing minimum
  and maximum lags, high frequency if \code{m=1}, low
  frequency otherwise.}

  \item{m}{the frequency ratio}

  \item{start}{the starting values for the weights of the
  one low frequency lag}
}
\value{
  a \code{lws_table} object, a list with elements
  \code{weights}, \code{lags} and \code{starts}
}
\description{
  Create table of weights, lags and starting values for
  Ghysels weight schema, see \link{amweights}
}
\details{
  Given weight function creates lags starting from
  \code{kmin} to \code{kmax} and replicates starting values
  for each low frequency lag.
}
\examples{
expand_amidas("nealmon","A",0,c(1,2),12,c(0,0,0))
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

