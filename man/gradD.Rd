\name{gradD}
\alias{gradD}
\title{Gradient function of coefficient function in MIDAS regression}
\usage{
  gradD(x, gr = NULL, ...)
}
\arguments{
  \item{x}{\code{\link{midas_r}} object}

  \item{gr}{a function or a list of functions for each
  restriction.}

  \item{...}{arguments supplied to gradient function(s)}
}
\value{
  a function
}
\description{
  Returns gradient function of coefficient function
}
\details{
  MIDAS regression can be written in the following form:

  \deqn{y=Xf(\eta)+u}

  This function calculates

  \deqn{\frac{\partial f}{\partial \eta}}
}
\author{
  Vaidotas Zemlys
}

