% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/midas_sp.R
\name{midas_pl_plain}
\alias{midas_pl_plain}
\title{MIDAS Partialy linear non-parametric regression}
\usage{
midas_pl_plain(
  y,
  X,
  z,
  p.ar = NULL,
  weight,
  degree = 1,
  start_bws,
  start_x,
  start_ar = NULL,
  method = c("Nelder-Mead"),
  ...
)
}
\arguments{
\item{y}{model response}

\item{X}{prepared matrix of high frequency variable lags for MMM term}

\item{z}{a vector, data for the non-parametric part}

\item{p.ar}{length of AR part}

\item{weight}{the weight function}

\item{degree}{the degree of local polynomial}

\item{start_bws}{the starting values bandwith}

\item{start_x}{the starting values for weight function}

\item{start_ar}{the starting values for AR part. Should be the same length as \code{p}}

\item{method}{a method passed to \link{optim}}

\item{...}{additional parameters to \link{optim}}
}
\value{
an object similar to \code{midas_r} object
}
\description{
Function for fitting PL MIDAS regression without the formula interface
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}
