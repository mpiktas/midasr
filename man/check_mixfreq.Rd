\name{check_mixfreq}
\alias{check_mixfreq}
\title{Check data for MIDAS regression}
\usage{
  check_mixfreq(lowfreq, highfreq)
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
  Given low and high frequency data check whether high
  frequency data can be converted to low frequency.
}
\details{
  If m is a frequency ratio, and n is the number of data
  points for low frequency data, then there should be n*m
  data points for high frequency data. This function checks
  whether this is the case. This function is used to
  prepare data for MIDAS regression and in general should
  not be interesting to ordinary users.
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

