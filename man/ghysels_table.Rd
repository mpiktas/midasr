\name{ghysels_table}
\alias{ghysels_table}
\title{Weight and lag selection table for MIDAS regression model with Ghysels schema}
\usage{
  ghysels_table(formula, data, weights, wstart, type,
    start = NULL, kmin = NULL, kmax = NULL,
    IC = c("AIC", "BIC"), test = c("hAh.test"),
    Ofunction = "optim", user.gradient = FALSE, ...)
}
\arguments{
  \item{formula}{the formula for MIDAS regression, the lag
  selection is performed for the last MIDAS lag term in the
  formula}

  \item{data}{a list containing data with mixed
  frequencies}

  \item{weights}{the names of weights used in Ghysels
  schema}

  \item{wstart}{the starting values for the weights of the
  firs low frequency lag}

  \item{type}{the type of Ghysels schema see
  \link{ghyselslag}, can be a vector of types}

  \item{start}{the starting values for optimisation
  excluding the starting values for the last term}

  \item{kmin}{the minimum low frequency lag, defaults to
  1.}

  \item{kmax}{the highest low frequency lag, defaults to
  square root of number of low fequency observations}

  \item{IC}{the names of information criteria which should
  be calculated}

  \item{test}{the names of statistical tests to perform on
  restricted model, p-values are reported in the columns of
  model selection table}

  \item{Ofunction}{see \link{midasr}}

  \item{user.gradient}{see \link{midas_r}}

  \item{...}{additional parameters to optimisation
  function, see \link{midas_r}}
}
\value{
  a \code{midas_r_ic_table} object which is the list with
  the following elements:

  \item{table}{the table where each row contains calculated
  information criteria for both restricted and unrestricted
  MIDAS regression model with given lag structure}
  \item{candlist}{the list containing fitted models}
  \item{IC}{the argument IC} \item{test}{the argument test}
  \item{weights}{the names of weight functions}
  \item{lags}{the lags used in models} ##'
}
\description{
  Create weight and lag selection table for the MIDAS
  regression model with Ghysels schema
}
\details{
  This function estimates models sequentialy increasing the
  midas lag from \code{kmin} to \code{kmax} and varying the
  weights of the last term of the given formula

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

tb <- ghysels_table(y~trend+fmls(x,12,12,nealmon),data=list(y=y,x=x,trend=trend),weights=c("nealmon"),wstart=list(nealmon=c(0,0,0)),start=list(trend=1),type=c("A"),kmax=4)
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

