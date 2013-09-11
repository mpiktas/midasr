\name{select_and_forecast}
\alias{select_and_forecast}
\title{Create table for different forecast horizons}
\usage{
  select_and_forecast(formula, data, from, to, insample,
    outsample, weights, wstart, start = NULL, IC = "AIC",
    seltype = c("restricted", "unrestricted"),
    test = "hAh.test",
    ftype = c("fixed", "recursive", "rolling"),
    measures = c("MSE", "MAPE", "MASE"),
    fweights = c("EW", "BICW", "MSFE", "DMSFE"), ...)
}
\arguments{
  \item{formula}{initial formula for the}

  \item{data}{list of data}

  \item{from}{a named list of starts of lags from where to
  fit. Denotes the horizon}

  \item{to}{a named list for lag selections}

  \item{insample}{the low frequency indexes for in-sample
  data}

  \item{outsample}{the low frequency indexes for
  out-of-sample data}

  \item{weights}{names of weight function candidates}

  \item{wstart}{starting values for weight functions}

  \item{start}{other starting values}

  \item{IC}{name of information criteria to choose model
  from}

  \item{seltype}{argument to modsel, \code{"restricted"}
  for model selection based on information criteria of
  restricted MIDAS model, \code{"unrestricted"} for model
  selection based on unrestricted (U-MIDAS) model.}

  \item{ftype}{which type of forecast to use.}

  \item{test}{argument to modsel}

  \item{measures}{the names of goodness of fit measures}

  \item{fweights}{names of weighting schemes}

  \item{...}{additional arguments for optimisation method,
  see \link{midas_r}}
}
\value{
  a list containing forecasts, tables of accuracy measures
  and the list with selected models
}
\description{
  Creates tables for different forecast horizons and table
  for combined forecasts
}
\details{
  Divide data into in-sample and out-of-sample. Fit
  different forecasting horizons for in-sample data.
  Calculate accuracy measures for individual and average
  forecasts.
}
\examples{
### Sets a seed for RNG ###
set.seed(1001)
## Number of low-frequency observations
n<-250
## Linear trend and higher-frequency explanatory variables (e.g. quarterly and monthly)
trend<-c(1:n)
x<-rnorm(4*n)
z<-rnorm(12*n)
## Exponential Almon polynomial constraint-consistent coefficients
fn.x <- nealmon(p=c(1,-0.5),d=8)
fn.z <- nealmon(p=c(2,0.5,-0.1),d=17)
## Simulated low-frequency series (e.g. yearly)
y<-2+0.1*trend+mls(x,0:7,4)\%*\%fn.x+mls(z,0:16,12)\%*\%fn.z+rnorm(n)
##Do not run
}
\author{
  Virmantas Kvedaras, Vaidotas Zemlys
}

