\name{ghyselslag}
\alias{ghyselslag}
\title{Special combination of weights}
\usage{
  ghyselslag(p, d, m, weight = nealmon,
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
  Combines given weight function into special combination
  used by E. Ghysels.
}
\details{
  Given a weight function \eqn{w(\beta,\theta)} which has a
  property of being defined as \eqn{\beta g(\theta)} the
  following combinations are defined.
  \deqn{(w(\beta_1,\theta_1),...,w(\beta_k,theta_k))}
  \deqn{(w(\beta_1,\theta),...,w(\beta_k,\theta))}
  \deqn{\beta(w(1,\theta_1),...,w(1,\theta_k))}
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

