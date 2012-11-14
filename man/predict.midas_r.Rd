\name{predict.midas_r}
\alias{predict.midas_r}
\title{Predict method for MIDAS regression fit}
\usage{
  \method{predict}{midas_r} (object, newldata, newhdata,
    ...)
}
\arguments{
  \item{object}{\code{\link{midas_r}} object}

  \item{newldata}{new low frequency data, must be
  \code{data.frame}}

  \item{newhdata}{new high frequency data, must be
  \code{data.frame}}

  \item{...}{additional arguments, not used}
}
\value{
  a vector of predicted values
}
\description{
  Predicted values based on MIDAS regression object
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

