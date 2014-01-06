\name{predict.midas_r}
\alias{predict.midas_r}
\title{Predict method for MIDAS regression fit}
\usage{
\method{predict}{midas_r}(object, newdata, na.action = na.omit, ...)
}
\arguments{
  \item{object}{\code{\link{midas_r}} object}

  \item{newdata}{a named list containing data for mixed
  frequencies. If omitted, the fitted values are used.}

  \item{na.action}{function determining what should be done
  with missing values in \code{newdata}. The most likely
  cause of missing values is the insufficient data for the
  lagged variables. The default is to omit such missing
  values.}

  \item{...}{additional arguments, not used}
}
\value{
a vector of predicted values
}
\description{
Predicted values based on \code{midas_r} object.
}
\details{
\code{predict.midas_r} produces predicted values, obtained
by evaluating regression function in the frame
\code{newdata}. This means that the appropriate model
matrix is constructed using only the data in
\code{newdata}. This makes this function not very
convenient for forecasting purposes. If you want to supply
the new data for forecasting horizon only use the function
\link{forecast.midas_r}. Also this function produces only
static predictions, if you want dynamic forecasts use the
\link{forecast.midas_r}.
}
\examples{
data("USrealgdp")
data("USunempr")

y <- diff(log(USrealgdp))
x <- window(diff(USunempr), start = 1949)

##24 high frequency lags of x included
mr <- midas_r(y ~ fmls(x, 23, 12, nealmon), start = list(x = rep(0, 3)))

##Declining unemployment
xn <- rnorm(2 * 12, -0.1, 0.1)

##Only one predicted value, historical values discarded
predict(mr, list(x = xn))

##Historical values taken into account
forecast(mr, list(x = xn))
}
\author{
Virmantas Kvedaras, Vaidotas Zemlys
}

