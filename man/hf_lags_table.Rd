\name{hf_lags_table}
\alias{hf_lags_table}
\title{Create a high frequency lag selection table for MIDAS regression model}
\usage{
hf_lags_table(formula, data, start, from, to, IC = c("AIC", "BIC"),
  test = c("hAh.test"), Ofunction = "optim", user.gradient = FALSE, ...)
}
\arguments{
  \item{formula}{the formula for MIDAS regression, the lag
  selection is performed for the last MIDAS lag term in the
  formula}

  \item{data}{a list containing data with mixed
  frequencies}

  \item{start}{the starting values for optimisation}

  \item{from}{a named list, or named vector with lag
  numbers which are the beginings of MIDAS lag structures.
  The names should correspond to the MIDAS lag terms in the
  formula for which to do the lag selection. Value NA
  indicates lag start at zero}

  \item{to}{a named list where each element is a vector
  with two elements. The first element is the lag number
  from which the lag selection starts, the second is the
  lag number at which the lag selection ends. NA indicates
  lowest (highest) lag numbers possible.}

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
a \code{midas_r_iclagtab} object which is the list with the
following elements:

\item{table}{the table where each row contains calculated
information criteria for both restricted and unrestricted
MIDAS regression model with given lag structure}
\item{candlist}{the list containing fitted models}
\item{IC}{the argument IC}
}
\description{
Creates a high frequency lag selection table for MIDAS
regression model with given information criteria and
minimum and maximum lags.
}
\details{
This function estimates models sequentially increasing the
midas lag from \code{kmin} to \code{kmax} of the last term
of the given formula
}
\examples{
data("USunempr")
data("USrealgdp")
y <- diff(log(USrealgdp))
x <- window(diff(USunempr),start=1949)
trend <- 1:length(y)

mlr <- hf_lags_table(y~trend+fmls(x,12,12,nealmon),
                     start=list(x=rep(0,3)),
                     from=c(x=0),to=list(x=c(4,4)))
mlr
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}

