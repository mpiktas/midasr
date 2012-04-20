##' Mixed data sampling lag
##'
##' .. content for \details{} ..
##' @title 
##' @param x 
##' @param ... 
##' @return 
##' @author Vaidotas Zemlys
mdslag <- function(x,...)UseMethod("mdslag")

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param x 
##' @param k 
##' @param m 
##' @param ... 
##' @return 
##' @author Vaidotas Zemlys
mdslag.default <- function(x, k, m, ...) {
    n.x <- length(x)
    n <- n.x %/%m
    if(length(k)>1) {
        klags <- k
        k <- max(k)
    }
    else {
        klags <- 0:k
    }
        
    idx <- m*c((n.x/m-n+k+1):(n.x/m))    

    X <- foreach(h.x=0:((k+1)*m-1), .combine='cbind') %do% {
        x[idx-h.x]
    }   
    nmx <- attr(x,"highfreqn")
    if(is.null(nmx)) nmx <- "X"
   
    colnames(X) <- paste(nmx, ".", rep(0:k, each=m), ".", rep(m:1, k+1), sep="")

    ###select only the lags needed
    lagn <- unlist(lapply(klags,function(no) grep(paste(nmx,".",no,"[.]",sep=""),colnames(X),value=TRUE)))
       
    X[,lagn]
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param x 
##' @param k 
##' @param ... 
##' @return 
##' @author Vaidotas Zemlys
mdslag.mds <- function(x,k,...) {
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
