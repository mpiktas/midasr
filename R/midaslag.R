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
    n <- n.x %/%m

    lk <- k
    k <- max(k)+1    
    
    if(n.x%%m != 0) stop("Incomplete high frequency data")
    idx <- m*(((k-1)%/%m+1):n)    

    X <- lapply(0:(k-1),function(h.x)x[idx-h.x])
    X <- do.call("cbind",X)
    
    if(k==1) X <- matrix(X,ncol=1)
    
    colnames(X) <- paste0("X", ".", 1:k-1,"/","m")
    padd <- matrix(NA,nrow=n-nrow(X),ncol=ncol(X))
    res <- rbind(padd,X)
    res[,lk+1,drop=FALSE]
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
dmls <- function(x,k,m,...) {
    dx <- c(NA,diff(x))
    v <- fmls(dx,k,m)
    colnames(v) <- gsub("X","DX",colnames(v))
    v
}
    

##' Check data for MIDAS regression
##'
##' Given mixed frequency data check whether higher frequency data can be converted to the lowest frequency.
##' 
##' @param data a list containing mixed frequency data
##' @return a boolean TRUE, if mixed frequency data is conformable, FALSE if it is not.
##' @details The number of observations in higher frequency data elements should have a common divisor with the number of observations in response variable. It is always assumed that the response variable is of the lowest frequency. 
##' 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
check_mixfreq <- function(data) {

    obsno <- sapply(data,function(l)ifelse(is.null(dim(l)),length(l),nrow(l)))
    m <- obsno %% min(obsno)

    sum(m>0)==0
}
