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
##' @details This is a convenience function, it calls \code{link{msl}} to perform actual calculations. 
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
##' ## Quartely frequency data
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
    n <- n.x %/%m

    lk <- k
    k <- max(k)+1    
    
    if(n.x%%m != 0) stop("Incomplete high frequency data")
    idx <- m*(((k-1)%/%m+1):n)    

    
    X <- foreach(h.x=0:(k-1), .combine='cbind') %do% {
        x[idx-h.x]
    }
    
    if(k==1) X <- matrix(X,ncol=1)
    
    colnames(X) <- paste0("X", ".", 1:k-1,"/","m")
    padd <- matrix(NA,nrow=n-nrow(X),ncol=ncol(X))
    res <- rbind(padd,X)
    res[,lk+1,drop=FALSE]
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
