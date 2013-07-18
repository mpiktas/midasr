\name{lagsel}
\alias{lagsel}
\title{Select the model based on given information criteria}
\usage{
  lagsel(x, IC = x$IC[1],
    type = c("restricted", "unrestricted"))
}
\arguments{
  \item{x}{and output from iclagtab function}

  \item{IC}{the name of information criteria to choose}

  \item{type}{the type of MIDAS model, either restricted or
  unrestricted}
}
\value{
  the best model based on information criteria,
  \link{midas_r} object
}
\description{
  Selects the model with minimum of given information
  criteria and model type
}
\details{
  This function selects the model from the lag selection
  table for which the chosen information criteria achieves
  the smallest value
}
\examples{
data("USunempr")
data("USrealgdp")
y <- diff(log(USrealgdp))
x <- window(diff(USunempr),start=1949)
trend <- 1:length(y)

mlr <- iclagtab(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=4,kmax=6)
lagsel(mlr,"BIC","unrestricted")
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

