##' Full MIDAS lag structure
##'
##' Create a matrix of MIDAS lags, including contemporaneous lag up to selected order.
##'
##' @param x a vector
##' @param k maximum lag order
##' @param m frequency ratio
##' @param ... further arguments
##' @return a matrix containing the lags
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @seealso mls
##' @details This is a convenience function, it calls \code{link{mls}} to perform actual calculations.
##' @export
##'
fmls <- function(x, k, m, ...) {
  mls(x, 0:k, m, ...)
}

##' MIDAS lag structure
##'
##' Create a matrix of selected MIDAS lags
##'
##' @param x a vector
##' @param k a vector of lag orders, zero denotes contemporaneous lag.
##' @param m frequency ratio
##' @param ... further arguments used in fitting MIDAS regression
##' @return a matrix containing the lags
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @details The function checks whether high frequency data is complete, i.e. \code{m} must divide \code{length(x)}.
##' @examples
##' ## Quarterly frequency data
##' x <- 1:16
##' ## Create MIDAS lag for use with yearly data
##' mls(x,0:3,4)
##'
##' ## Do not use contemporaneous lag
##' mls(x,1:3,4)
##'
##' ## Compares with embed when m=1
##' embed(x,2)
##' mls(x,0:1,1)
##' @export
mls <- function(x, k, m, ...) {
  n.x <- length(x)
  n <- n.x %/% m

  lk <- k
  k <- max(k) + 1
  mk <- min(k)
  if (mk > 0) mk <- 0

  if (n.x %% m != 0) stop("Incomplete high frequency data")
  idx <- m * (((k - 1) %/% m + 1):(n - mk))

  X <- lapply(mk:(k - 1), function(h.x) x[idx - h.x])
  X <- do.call("cbind", X)

  if (k == 1) X <- matrix(X, ncol = 1)

  colnames(X) <- paste0("X", ".", (mk + 1):k - 1, "/", "m")
  padd <- matrix(NA, nrow = n - nrow(X), ncol = ncol(X))
  res <- rbind(padd, X)
  res[, lk + 1, drop = FALSE]
}


##' MIDAS lag structure for unit root processes
##'
##' Prepares MIDAS lag structure for unit root processes
##' @param x a vector
##' @param k maximal lag order
##' @param m frequency ratio
##' @param ... further arguments used in fitting MIDAS regression
##' @return a matrix containing the first differences and the lag k+1.
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
dmls <- function(x, k, m, ...) {
  dx <- c(NA, diff(x))
  v <- fmls(dx, k, m)
  colnames(v) <- gsub("X", "DX", colnames(v))
  v
}


##' Check data for MIDAS regression
##'
##' Given mixed frequency data check whether higher frequency data can be converted to the lowest frequency.
##'
##' @param data a list containing mixed frequency data
##' @return a boolean TRUE, if mixed frequency data is conformable, FALSE if it is not.
##' @details The number of observations in higher frequency data elements should have a common divisor
##'  with the number of observations in response variable. It is always assumed that the response variable
##'  is of the lowest frequency.
##'
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
check_mixfreq <- function(data) {
  obsno <- sapply(data, function(l) ifelse(is.null(dim(l)), length(l), nrow(l)))
  m <- obsno %% min(obsno)

  sum(m > 0) == 0
}


#' MIDAS lag structure with dates
#'
#' @param x a vector, of high frequency time series. Must be zoo or ts object
#' @param k lags, a vector
#' @param y a vector of low frequency time series. Must be zoo or ts object
#' @param ... further arguments used in fitting MIDAS regression
#'
#' @return a matrix containing the lags
#' @details High frequency time series is aligned with low frequency time series using date information.
#' Then the high frequency lags are calculated.
#'
#' To align the time series the low frequency series index
#' needs to be extended by one low frequency period into the past and into the future. If supplied time series
#' object does not support extending time index, a simple heuristic is used.
#'
#' It is expected that time index for zoo objects can be converted to POSIXct format.
#'
#'
#' @author Virmantas Kvedaras, Vaidotas Zemlys-Balevičius
#' @importFrom stats lag
#' @export
#'
#' @examples
#'
#' # Example with ts objects
#' x <- ts(c(1:144), start = c(1980, 1), frequency = 12)
#' y <- ts(c(1:12), start = 1980, frequency = 1)
#'
#'
#' # msld and mls should give the same results
#'
#' m1 <- mlsd(x, 0:5, y)
#'
#' m2 <- mls(x, 0:5, 12)
#'
#' sum(abs(m1 - m2))
#'
#' # Example with zooreg
#'
#' # Convert x to zooreg object using yearmon time index
#' \dontrun{
#' xz <- zoo::as.zooreg(x)
#'
#' yz <- zoo::zoo(as.numeric(y), order.by = as.Date(paste0(1980 + 0:11, "-01-01")))
#'
#' # Heuristic works here
#' m3 <- mlsd(xz, 0:5, yz)
#'
#' sum(abs(m3 - m1))
#' }
mlsd <- function(x, k, y, ...) {
  datex <- get_datex(x)

  datey <- get_datey(y, datex)

  x <- as.numeric(x)

  ct <- cut(datex, datey, right = FALSE, labels = FALSE, include.lowest = TRUE)
  tct <- table(ct)
  uct <- unique(ct)
  nuct <- na.omit(uct)
  # We do not need the first period, but sometimes it is matched.
  # In that case it is dropped.
  id <- match(2:(length(datey) - 1), nuct)

  fhx <- function(h.x) {
    id <- h.x - k
    id[id <= 0] <- NA
    x[id]
  }
  XX <- lapply(cumsum(tct), fhx)
  X <- do.call("rbind", XX)
  colnames(X) <- paste("X", k, sep = ".")
  X[id, ]
}

get_datex <- function(x) UseMethod("get_datex")

get_datex.zoo <- function(x) {
  as.POSIXct(index(x))
}

get_datex.ts <- function(x) {
  time(x)
}

get_datey <- function(y, datex) UseMethod("get_datey")

get_datey.ts <- function(y, datex) {
  left <- time(lag(y, 1))[1]
  right <- tail(time(lag(y, -1)), n = 1)
  if (datex[1] < left) left <- datex[1]
  c(left, time(y), right) - 0.001
}

get_datey.default <- function(datey, datex) {
  ## If we get here, we assume that both datey and datex are ordered and comparable
  left <- datey[1] - (datey[2] - datey[1])
  if (datex[1] < left) left <- datex[1]
  nd <- length(datey)
  right <- datey[nd] + (datey[nd] - datey[nd - 1])
  c(left, datey, right)
}

get_datey.zoo <- function(datey, datex) {
  ## Test whether the lag extends the dates
  lagy <- lag(datey, 1)
  fd_lagy <- index(lagy)[1]
  fd_y <- index(datey)[1]
  datey_p <- as.POSIXct(index(datey))
  ## If the dates are not extended use heuristic for left and right low frequency dates
  if (fd_lagy == fd_y) {
    get_datey.default(datey_p, datex)
  } else {
    left <- as.POSIXct(fd_lagy)
    right <- as.POSIXct(tail(index(lag(datey, -1)), n = 1))
    if (datex[1] < left) left <- datex[1]
    return(c(left, datey_p, right))
  }
}
