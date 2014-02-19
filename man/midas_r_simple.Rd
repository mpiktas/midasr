\name{midas_r_simple}
\alias{midas_r_simple}
\title{Restricted MIDAS regression}
\usage{
midas_r_simple(y, X, z = NULL, weight, grw = NULL, startx, startz = NULL,
  method = c("Nelder-Mead", "BFGS"), ...)
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
an object similar to \code{midas_r} object
}
\description{
Function for fitting MIDAS regression without the formula
interface
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}

