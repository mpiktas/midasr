% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lagspec.R
\name{amweights}
\alias{amweights}
\title{Weights for aggregates based MIDAS regressions}
\usage{
amweights(p, d, m, weight = nealmon, type = c("A", "B", "C"))
}
\arguments{
\item{p}{parameters for weight functions, see details.}

\item{d}{number of high frequency lags}

\item{m}{the frequency}

\item{weight}{the weight function}

\item{type}{type of structure, a string, one of A, B or C.}
}
\value{
a vector of weights
}
\description{
Produces weights for aggregates based MIDAS regression
}
\details{
Suppose a weight function \eqn{w(\beta,\theta)} satisfies the following equation:
\deqn{w(\beta,\theta)=\beta g(\theta)}

The following combinations are defined, corresponding to structure types
\code{A}, \code{B} and \code{C} respectively:
\deqn{(w(\beta_1,\theta_1),...,w(\beta_k,\theta_k))}
\deqn{(w(\beta_1,\theta),...,w(\beta_k,\theta))}
\deqn{\beta(w(1,\theta),...,w(1,\theta)),}

where \eqn{k} is the number of low frequency lags, i.e. \eqn{d/m}. If the latter
value is not whole number, the error is produced.


The starting values \eqn{p} should be supplied then as follows:
\deqn{(\beta_1,\theta_1,...,\beta_k,\theta_k)}
\deqn{(\beta_1,...,\beta_k,\theta)}
\deqn{(\beta,\theta)}
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}
