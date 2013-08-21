\name{forecast}
\alias{forecast}
\alias{forecast.midas_r}
\title{Forecast MIDAS regression}
\usage{
  forecast(object, ...)

  \method{forecast}{midas_r} (object, newdata = NULL,
    method = c("static", "dynamic"),
    insample = get_estimation_sample(object), ...)
}
\arguments{
  \item{object}{midas_r object}

  \item{newdata}{newdata}

  \item{method}{the forecasting method, either
  \code{"static"} or \code{"dynamic"}}

  \item{insample}{a list containing the historic mixed
  frequency data}

  \item{...}{additional arguments, not used}
}
\value{
  a vector of forecasts
}
\description{
  Forecasts MIDAS regression given the future values of
  regressors. For dynamic models (with lagged response
  variable) there is an option to calculate dynamic
  forecast, when forecasted values of response variable are
  substituted into the lags of response variable.
}
\details{
  Given future values of regressors this function combines
  the historical values used in the fitting the MIDAS
  regression model and calculates the forecasts.
}
\author{
  Vaidotas Zemlys
}

