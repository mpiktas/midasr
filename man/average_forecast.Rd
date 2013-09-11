\name{average_forecast}
\alias{average_forecast}
\title{Average forecasts of MIDAS models}
\usage{
  average_forecast(modlist, data, insample, outsample,
    type = c("fixed", "recursive", "rolling"),
    fweights = c("EW", "BICW", "MSFE", "DMSFE"),
    measures = c("MSE", "MAPE", "MASE"))
}
\arguments{
  \item{modlist}{a list of \code{midas_r} objects}

  \item{data}{a list with mixed frequency data}

  \item{insample}{the low frequency indexes for in-sample
  data}

  \item{outsample}{the low frequency indexes for
  out-of-sample data}

  \item{type}{a string indicating which type of forecast to
  use.}

  \item{fweights}{names of weighting schemes}

  \item{measures}{names of accuracy measures}
}
\value{
  a list containing forecasts and tables of accuracy
  measures
}
\description{
  Average MIDAS model forecasts using specified weighting
  scheme. Produce in-sample and out-of-sample accuracy
  measures.
}
\details{
  Given the data, split it to in-sample and out-of-sample
  data. Then given the list of models, reestimate each
  model with in-sample data and produce out-of-sample
  forecast. Given the forecasts average them with the
  specified weighting scheme. Then calculate the accuracy
  measures for individual and average forecasts.

  The forecasts can be produced in 3 ways. The
  \code{"fixed"} forecast uses model estimated with
  in-sample data. The \code{"rolling"} forecast reestimates
  model each time by increasing the in-sample by one low
  frequency observation and dropping the first low
  frequency observation. These reestimated models then are
  used to produce out-of-sample forecasts. The
  \code{"recursive"} forecast differs from \code{"rolling"}
  that it does not drop observations from the beginning of
  data.
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

