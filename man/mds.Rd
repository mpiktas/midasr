\name{mds}
\alias{mds}
\title{Prepare data for MIDAS regression}
\usage{
  mds(lowfreq, highfreq)
}
\arguments{
  \item{lowfreq}{\code{data.frame} object containing low
  frequency data}

  \item{highfreq}{\code{data.frame} object containing high
  frequency data}
}
\value{
  a list with elements \code{lowfreq} and \code{highfreq}
}
\description{
  Given low and high frequency data calculate frequency
  ratio and store this information.
}
\details{
  This function converts each column of the given
  \code{data.frame} to \code{mds} object. This means adding
  attribute \code{frequency.ratio}. This function is used
  to prepare data for MIDAS regression and in general
  should not be interesting to ordinary users.
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

