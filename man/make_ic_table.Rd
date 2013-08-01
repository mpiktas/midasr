\name{make_ic_table}
\alias{make_ic_table}
\title{Create a weight and lag selection table for MIDAS regression model}
\usage{
  make_ic_table(formula, data = NULL, start = NULL, table,
    wstart, IC = c("AIC", "BIC"), test = c("hAh.test"),
    Ofunction = "optim", user.gradient = FALSE, ...)
}
\arguments{
  \item{formula}{the formula for MIDAS regression, the lag
  selection is performed for the last MIDAS lag term in the
  formula}

  \item{data}{a list containing data with mixed
  frequencies}

  \item{start}{the starting values for optimisation
  excluding the starting values for the last term}

  \item{table}{a list with first element a list of lag
  structures and the second element the vector of weight
  function names. The second element is recycled if its
  length does not coincide with the length of the first
  element.}

  \item{wstart}{the starting values for weight functions, a
  named list.}

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
  a \code{midas_r_ICtable} object which is the list with
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


mwlr <- make_ic_table(y~trend+fmls(x,12,12,nealmon),table=list(list(0:4,0:5,0:6),c("nealmon","nealmon","nbeta")),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)))

mwlr
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

