\name{mmweights}
\alias{mmweights}
\title{Weights for multiplicative MIDAS regressions}
\usage{
  mmweights(p, d, m, weight = nealmon,
    type = c("A", "B", "C"))
}
\arguments{
  \item{p}{parameters for weight functions, see details.}

  \item{d}{number of lags}

  \item{m}{the frequency}

  \item{weight}{the weight function}

  \item{type}{type of structure, a string, one of A, B or
  C.}
}
\value{
  a vector of weights
}
\description{
  Produces weights for multiplicative MIDAS regression
}
\details{
  Given a weight function \eqn{w(\beta,\theta)} which has a
  property of being defined as \eqn{\beta g(\theta)} the
  following combinations are defined, corresponding to
  structure types \code{A}, \code{B} and \code{C}
  respectively:
  \deqn{(w(\beta_1,\theta_1),...,w(\beta_k,\theta_k))}
  \deqn{(w(\beta_1,\theta),...,w(\beta_k,\theta))}
  \deqn{\beta(w(1,\theta_1),...,w(1,\theta_k))}

  The starting values \eqn{p} should be supplied then as
  follows: \deqn{(\beta_1,\theta_1,...,\beta_k,\theta_k)}
  \deqn{(\beta_1,...,\beta_k,\theta)}
  \deqn{(\beta,\theta_1,...,\theta_k)}
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

