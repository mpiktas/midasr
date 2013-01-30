\name{imidas_r.imidas_r}
\alias{imidas_r.imidas_r}
\title{Restricted MIDAS regression with I(1) regressors}
\usage{
  \method{imidas_r}{imidas_r} (x, start = coef(x),
    Ofunction = x$argmap.opt$Ofunction, ...)
}
\arguments{
  \item{x}{\code{imidas_r} object}

  \item{start}{the starting values}

  \item{Ofunction}{a character string of the optimisation
  function to use. The default value is to use the function
  of previous optimisation.}

  \item{...}{further arguments to optimisation function. If
  none are supplied, the arguments of previous optimisation
  are used.}
}
\value{
  \code{imidas_r} object
}
\description{
  Reestimate the MIDAS regression with I(1) regressors with
  different starting values
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}
\seealso{
  imidas_r
}

