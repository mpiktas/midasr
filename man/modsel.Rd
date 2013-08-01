\name{modsel}
\alias{modsel}
\title{Select the model based on given information criteria}
\usage{
  modsel(x, IC = x$IC[1],
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
  This function selects the model from the model selection
  table for which the chosen information criteria achieves
  the smallest value. The function works with model tables
  produced by functions \link{lf_lags_table},
  \link{hf_lags_table} and \link{make_ic_table}.
}
\examples{
data("USunempr")
data("USrealgdp")
y <- diff(log(USrealgdp))
x <- window(diff(USunempr),start=1949)
trend <- 1:length(y)

mhfr <- hf_lags_table(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=4,kmax=6)

mlfr <- lf_lags_table(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=1,kmax=2)

mwlr <- make_ic_table(y~trend+fmls(x,12,12,nealmon),table=list(list(0:4,0:5,0:6),c("nealmon","nealmon","nbeta")),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)))

modsel(mhfr,"BIC","unrestricted")

modsel(mlfr,"BIC","unrestricted")

modsel(mwlr,"BIC","unrestricted")
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

