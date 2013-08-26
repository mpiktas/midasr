\name{midas_r_fast}
\alias{midas_r_fast}
\title{Restricted MIDAS regression}
\usage{
  midas_r_fast(y, X, z = NULL, weight, grw = NULL, startx,
    startz = NULL, method = "BFGS", ...)
}
\arguments{
  \item{y}{model response}

  \item{X}{prepared matrix of high frequency variable lags}

  \item{z}{additional low frequency variables}

  \item{weight}{the weight function}

  \item{grw}{the gradient of weight function}

  \item{startx}{the starting values for weight function}

  \item{startz}{the starting values for additional low
  frequency variables}

  \item{method}{a method passed to \link{optim}}

  \item{...}{additional parameters to \link{optim}}
}
\value{
  a \code{midas_r} object
}
\description{
  Fast function for fitting MIDAS regression
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

