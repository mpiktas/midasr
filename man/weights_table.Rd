\name{weights_table}
\alias{weights_table}
\title{Create a weight function selection table for MIDAS regression model}
\usage{
  weights_table(formula, data, start = NULL,
    IC = c("AIC", "BIC"), test = c("hAh.test"),
    Ofunction = "optim", user.gradient = FALSE, ...)
}
\arguments{
  \item{formula}{the formula for MIDAS regression, the lag
  selection is performed for the last MIDAS lag term in the
  formula}

  \item{data}{a list containing data with mixed
  frequencies}

  \item{start}{the starting values for optimisation}

  \item{IC}{the information criteria which to compute}

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
  \item{IC}{the argument IC}
}
\description{
  Creates a weight function selection table for MIDAS
  regression model with given information criteria and
  weight functions.
}
\details{
  This function estimates models sequentially increasing
  the midas lag from \code{kmin} to \code{kmax} of the last
  term of the given formula
}
\examples{
data("USunempr")
data("USrealgdp")
y <- diff(log(USrealgdp))
x <- window(diff(USunempr),start=1949)
trend <- 1:length(y)
mwr <- weights_table(y~trend+fmls(x,12,12,nealmon),start=list(x=list(nealmon=rep(0,3),nbeta=c(1,1,1,0))))

mwr
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

