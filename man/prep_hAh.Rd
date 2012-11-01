\name{prep_hAh}
\alias{prep_hAh}
\title{Calculate data for \link{hAh.test} and \link{hAhr.test}}
\usage{
  prep_hAh(x, gr, ...)
}
\arguments{
  \item{x}{\code{midas_r} object}

  \item{gr}{a function or a list of functions for each
  restriction.}

  \item{...}{arguments supplied to gradient function(s)}
}
\value{
  a list with necessary matrices
}
\description{
  Workhorse function for calculating necessary matrices for
  \link{hAh.test} and \link{hAhr.test}. Takes the same
  parameters as \link{hAh.test}
}
\author{
  Virmantas Kvedara, Vaidotas Zemlys
}
\seealso{
  hAh.test, hAhr.test
}

