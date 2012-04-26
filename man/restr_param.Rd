\name{restr_param}
\alias{restr_param}
\title{Return the estimated parameters of the restriction function(s)}
\usage{
  restr_param(x, name = restr_names(x))
}
\arguments{
  \item{x}{an output from \code{\link{midas_r}}}

  \item{name}{name of the restriction function, the default
  value is the names of the restriction functions supplied
  to \code{\link{midas_r}}}
}
\value{
  a list if \code{length(name)>1}, a vector otherwise
}
\description{
  A helper function for working with output of
  \code{\link{midas_r}}. Returns the estimated parameters
  of restriction function(s)
}
\author{
  Vaidotas Zemlys
}

