\name{midas_r_np}
\alias{midas_r_np}
\title{Estimate non-parametric MIDAS regression}
\usage{
midas_r_np(x, data, lambda = NULL)
}
\arguments{
  \item{x}{formula specifying MIDAS regression}

  \item{data}{a named list containing data with mixed
  frequencies}

  \item{lambda}{smoothing parameter, defaults to
  \code{NULL}, which means that it is chosen by minimising
  AIC.}
}
\value{
a \code{midas_r_np} object
}
\description{
Estimates non-parametric MIDAS regression
}
\details{
Estimates non-parametric MIDAS regression accodring
Breitung et al.
}
\examples{
data("USunempr")
data("USrealgdp")
y <- diff(log(USrealgdp))
x <- window(diff(USunempr),start=1949)
trend <- 1:length(y)
midas_r_np(y~fmls(x,12,12)-1)
}
\author{
Vaidotas Zemlys
}
\references{
Breitung J, Roling C, Elengikal S (2013). \emph{The
statistical content and empirical testing of the MIDAS
restrictions} Working paper, URL
http://www.ect.uni-bonn.de/mitarbeiter/
joerg-breitung/npmidas.
}

