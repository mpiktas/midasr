\name{combine_forecasts}
\alias{combine_forecasts}
\title{Create table for different forecast horizons}
\usage{
  combine_forecasts(formula, data, from, to, insample,
    outsample, weights, wstart, start = NULL, IC = "AIC",
    type = "restricted", test = "hAh.test",
    measures = c("MSE", "MAPE", "MASE"),
    fweights = c("EW", "BICW", "MSFE", "DMSFE"),
    fmethod = c("static", "dynamic"), ...)
}
\arguments{
  \item{formula}{initial formula for the}

  \item{data}{list of data}

  \item{from}{a named list of starts of lags from where to
  fit. Denotes the horizon}

  \item{to}{a named list for lag selections}

  \item{insample}{low-frequency in-sample indexes}

  \item{outsample}{low-frequency out-of-sample indexes}

  \item{weights}{names of weight function candidates}

  \item{wstart}{starting values for weight functions}

  \item{start}{other starting values}

  \item{IC}{name of information criteria to choose model
  from}

  \item{type}{argument to modsel}

  \item{test}{argument to modsel}

  \item{measures}{the names of goodness of fit measures}

  \item{fweights}{names of weighting schemes}

  \item{fmethod}{forecastin method, either static or
  dynamic}

  \item{...}{additional arguments for optimisation method,
  see \link{midas_r}}
}
\value{
  a list of tables, best models, in-sample data and
  out-of-sample data
}
\description{
  Creates tables for different forecast horizons and table
  for combined forecasts
}
\details{
  Divide data into in-sample and out-of-sample. Fit
  different forecasting horizons for in-sample data.
  Calculate various statistics for out-of-sample.
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

