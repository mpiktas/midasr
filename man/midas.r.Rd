\name{midas.r}
\alias{midas.r}
\title{Restricted MIDAS regression}
\usage{
  midas.r(y, x, resfun, start, method = "BFGS",
    control.optim = list(), ...)
}
\arguments{
  \item{y}{the response variable}

  \item{x}{the predictor variable}

  \item{resfun}{the function which returns restricted
  parameters given the restriction function. The parameters
  for the restriction function must be the supplied as
  numeric vector in the first argument of the function.
  Number of lags of the regression is calculated from the
  output of this function.}

  \item{gradfun}{the gradient of the restriction function.
  Must return the matrix with dimensions \code{l x p},
  where \code{l} is the number of unrestricted coefficients
  of the regression and \code{p} is the number of
  parameters in the restriction function.}

  \item{start}{the starting values for optimisation}

  \item{method}{the method used for optimisation, see
  \code{optim} documentation. All methods are suported
  except "L-BFGS-B" and "Brent". Default method is "BFGS".}

  \item{control.optim}{a list of control parameters for
  \code{optim}.}

  \item{...}{additional parameters supplied for
  \code{resfun} and \code{gradfun}}
}
\value{
  output suitable for function HAH.test
}
\description{
  Estimate restricted MIDAS regression using non-linear
  least squares. Uses \code{optim} for optimisation.
}
\author{
  Vaidotas Zemlys
}

