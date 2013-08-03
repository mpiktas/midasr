\name{expand_weights_lags}
\alias{expand_weights_lags}
\title{Create table of weights, lags and starting values}
\usage{
  expand_weights_lags(weights, kmin, kmax, m = 1,
    mkmin = 0, start)
}
\arguments{
  \item{weights}{either a vector with names of the weight
  functions or a named list of weight functions}

  \item{kmin}{the minimum lag, high frequency if
  \code{m=1}, low frequency otherwise}

  \item{kmax}{the highest lag, high frequency if
  \code{m=1}, low frequency otherwise}

  \item{m}{the frequency ratio}

  \item{mkmin}{the high frequency lags from which to start
  the fitting}

  \item{start}{a named list with the starting values for
  weight functions}
}
\value{
  a \code{lws_table} object, a list with elements
  \code{weights}, \code{lags} and \code{starts}.
}
\description{
  Creates table of weights, lags and starting values
}
\details{
  For each weight function creates lags starting from
  \code{kmin} to \code{kmax}. This is a convenience
  function for easier work with the function
  \link{midas_r_ic_table}.
}
\examples{
expand_weights_lags(c("nealmon","nbeta"),4,8,1,0,start=list(nealmon=rep(0,3),nbeta=rep(0,4)))
nlmn <- expand_weights_lags("nealmon",4,8,1,0,start=list(nealmon=rep(0,3)))
nbt <- expand_weights_lags("nbeta",4,8,1,0,start=list(nbeta=rep(0,4)))

nlmn+nbt
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

