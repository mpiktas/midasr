\name{forecast}
\alias{forecast}
\alias{forecast.midas_r}
\title{Forecast MIDAS regression}
\usage{
forecast(object, ...)

\method{forecast}{midas_r}(object, newdata = NULL, method = c("static",
  "dynamic"), insample = get_estimation_sample(object), ...)
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
variable) there is an option to calculate dynamic forecast,
when forecasted values of response variable are substituted
into the lags of response variable.
}
\details{
Given future values of regressors this function combines
the historical values used in the fitting the MIDAS
regression model and calculates the forecasts.
}
\examples{
data("USrealgdp")
data("USunempr")

y <- diff(log(USrealgdp))
x <- window(diff(USunempr), start = 1949)
trend <- 1:length(y)

##24 high frequency lags of x included
mr <- midas_r(y ~ trend + fmls(x, 23, 12, nealmon), start = list(x = rep(0, 3)))

##Forecast horizon
h <- 3
##Declining unemployment
xn <- rep(-0.1, 12*3)
##New trend values
trendn <- length(y) + 1:h

##Static forecasts combining historic and new high frequency data
forecast(mr, list(trend = trendn, x = xn), method = "static")

##Dynamic AR* model
mr.dyn <- midas_r(y ~ trend + mls(y, 1:2, 1, "*")
                   + fmls(x, 11, 12, nealmon),
                  start = list(x = rep(0, 3)))

forecast(mr.dyn, list(trend = trendn, x = xn), method = "dynamic")
}
\author{
Vaidotas Zemlys
}

