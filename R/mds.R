##' Lag a mixed data sampling time series
##'
##' Compute a lagged version of mixed data sampling time series
##' 
##' @param x a vector
##' @param ... further arguments
##' @return a matrix containing the lags
##' @author Vaidotas Zemlys
##' @seealso embedlf.default embedlf.mds
##' @export
embedlf <- function(x,...)UseMethod("embedlf")

##' Lag a mixed data sampling time series
##'
##' Compute a lagged version of mixed data sampling time series
##' 
##'
##' @param x a vector for which mixed data sampling lag has to be computed
##' @param k the number of high frequency lags plus 1.
##' @param m frequency ratio
##' @param ... further arguments
##' @return a matrix containing the lags
##' @author Vaidotas Zemlys
##' @method embedlf default
##' @seealso embedlf embedlf.mds
##' @export
embedlf.default <- function(x, k, m, ...) {
    n.x <- length(x)
    n <- n.x %/%m
             
    X <- foreach(t=ceiling(k/m):n, .combine='rbind') %do% {
        x[(m * t):(m * t - k + 1)]
    }   

    #nmx <- attr(x,"highfreqn")
    #if(is.null(nmx)) nmx <- "X"

    colnames(X) <- paste0("X", ".", 1:k-1,"/","m")

    padd <- matrix(NA,nrow=n-nrow(X),ncol=ncol(X))
    rbind(padd,X)
    
}

##' Lag a mixed data sampling time series
##'
##' Compute a lagged version of mixed data sampling time series
##' 
##'  
##' @param x an object of class \code{\link{mds}}
##' @param k the number of lags of low frequency. This can either be a number indicating the highest lag, or the vector with lags which should be included.
##' @param ... further arguments
##' @return a matrix containing lags
##' @author Vaidotas Zemlys
##' @export
##' @method embedlf mds
##' @seealso embedlf embedlf.default
embedlf.mds <- function(x,k,...) {
    m <- attr(x,"frequency.ratio")    
    n.x <- length(x)
        n <- n.x %/%m
        if(length(k)>1) {
            klags <- k
            k <- max(k)
        }
        else {
            klags <- 0:k
        }    
    if(m == 1) {
         c(rep(NA,k),x[1:(n-k)])      
    } else {
        idx <- m*c((n.x/m-n+k+1):(n.x/m))    

        X <- foreach(h.x=0:((k+1)*m-1), .combine='cbind') %do% {
            x[idx-h.x]
        }
    
        nmx <- attr(x,"varname")
        if(is.null(nmx)) nmx <- "X"
   
        colnames(X) <- paste(nmx, ".", rep(0:k, each=m), ".", rep(m:1, k+1), sep="")

###select only the lags needed
        lagn <- unlist(lapply(klags,function(no) grep(paste(nmx,".",no,"[.]",sep=""),colnames(X),value=TRUE)))
       
        padd <- matrix(NA,nrow=n-nrow(X),ncol=ncol(X))
        rbind(padd,X)
    }
}

##' Prepare data for MIDAS regression 
##'
##' Given low and high frequency data calculate frequency ratio and store this information. 
##' 
##' @param lowfreq \code{data.frame} object containing low frequency data 
##' @param highfreq \code{data.frame} object containing high frequency data
##' @return a list with elements \code{lowfreq} and \code{highfreq}
##' @details This function converts each column of the given \code{data.frame} to \code{mds} object. This means adding attribute \code{frequency.ratio}. This function is used to prepare data for MIDAS regression and in general should not be interesting to ordinary users.
##' 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
mds <- function(lowfreq,highfreq) {
    nl <- nrow(lowfreq)
    nh <- nrow(highfreq)

    m <- nh %/% nl

    if(nl*m!=nh) stop("Low and high frequency data sets are not conformable")
    lowfreq <- lapply(lowfreq,function(x) {
        class(x) <- c(class(x),"mds")
        attr(x,"frequency.ratio") <- 1
        x
    })

    highfreq <- lapply(highfreq,function(x) {
        class(x) <- c(class(x),"mds")
        attr(x,"frequency.ratio") <- m
        x
    })
    list(lowfreq=as.data.frame(lowfreq),highfreq=as.data.frame(highfreq))
}
