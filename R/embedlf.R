##' Lag a mixed data sampling time series
##'
##' Compute a lagged version of mixed data sampling time series
##' 
##' @param x a vector
##' @param k the number of high frequency lags including contemporaneous lag
##' @param m frequency ratio
##' @param ... further arguments
##' @return a matrix containing the lags
##' @author Vaidotas Zemlys
##' @seealso embedlf.default 
##' @export
embedlf <- function(x,...)UseMethod("embedlf")

##' @rdname embedlf
##' @method embedlf default
##' @export
embedlf.default <- function(x, k, m, ...) {
    n.x <- length(x)
    n <- n.x %/%m

    if(n.x%%m != 0) stop("Incomplete high frequency data")
    idx <- m*(((k-1)%/%m+1):n)    
    
    X <- foreach(h.x=0:(k-1), .combine='cbind') %do% {
        x[idx-h.x]
    }

    if(k==1) X <- matrix(X,ncol=1)
    
    colnames(X) <- paste0("X", ".", 1:k-1,"/","m")
    padd <- matrix(NA,nrow=n-nrow(X),ncol=ncol(X))
    rbind(padd,X)
    
}

##' Check data for MIDAS regression
##'
##' Given low and high frequency data check whether high frequency data can be converted to low frequency.
##' 
##' @param lowfreq \code{data.frame} object containing low frequency data 
##' @param highfreq \code{data.frame} object containing high frequency data
##' @return a list with elements \code{lowfreq} and \code{highfreq}
##' @details If m is a frequency ratio, and n is the number of data points for low frequency data, then there should be n*m data points for high frequency data. This function checks whether this is the case. 
##' This function is used to prepare data for MIDAS regression and in general should not be interesting to ordinary users.
##' 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
check_mixfreq <- function(lowfreq,highfreq) {
    nl <- nrow(lowfreq)
    nh <- nrow(highfreq)

    m <- nh %/% nl

    if(nl*m!=nh) stop("Low and high frequency data sets are not conformable")   
    list(lowfreq=as.data.frame(lowfreq),highfreq=as.data.frame(highfreq))
}
