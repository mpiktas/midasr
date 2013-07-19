\name{icwlagtab}
\alias{icwlagtab}
\title{Create a weight and lag selection table for MIDAS regression model}
\usage{
  icwlagtab(x, ldata = NULL, hdata = NULL, start = NULL,
    weights, wstart, kmin, kmax = NULL,
    IC = c("AIC", "BIC"), test = c("hAh.test"),
    Ofunction = "optim", user.gradient = FALSE, ...)
}
\arguments{
  \item{x}{the formula for MIDAS regression, the lag
  selection is performed for the last MIDAS lag term in the
  formula}

  \item{ldata}{low frequency data, a \code{data.frame}
  object}

  \item{hdata}{high frequency data, a \code{data.frame}
  object}

  \item{start}{the starting values for optimisation
  excluding the starting values for the last term}

  \item{weights}{the names of weight function}

  \item{wstart}{the starting values for weight functions, a
  named list.}

  \item{kmin}{the minimum high frequency lag, defaults to
  zero.}

  \item{kmax}{the highest high frequency lag, defaults to
  square root of number of low fequency observations}

  \item{IC}{the information criterias which to compute}

  \item{test}{the names of statistical tests to perform on
  restricted model, p-values are reported in the columns of
  model selection table}

  \item{Ofunction}{see \link{midasr}}

  \item{user.gradient}{see \link{midas_r}}

  \item{...}{additional parameters to optimisation
  function, see \link{midas_r}}
}
\value{
  a \code{midas_r_icwlagtab} object which is the list with
  the following elements:

  \item{table}{the table where each row contains calculated
  information criteria for both restricted and unrestricted
  MIDAS regression model with given lag structure}
  \item{candlist}{the list containing fitted models}
  \item{IC}{the argument IC}
}
\description{
  Creates a weight and lag selection table for MIDAS
  regression model with given information criteria and
  minimum and maximum lags.
}
\details{
  This function estimates models sequentialy increasing the
  midas lag from \code{kmin} to \code{kmax} and varying the
  weights of the last term of the given formula
}
\examples{
data("USunempr")
data("USrealgdp")
y <- diff(log(USrealgdp))
x <- window(diff(USunempr),start=1949)
trend <- 1:length(y)


mwlr <- icwlagtab(y~trend+fmls(x,12,12,nealmon),weights=c("nealmon","nbeta"),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)),kmin=4,kmax=6)

mwlr
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

