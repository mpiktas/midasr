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
    mk <- min(k)
    if(mk>0) mk <- 0
    
    if(n.x%%m != 0) stop("Incomplete high frequency data")
    idx <- m*(((k-1)%/%m+1):(n-mk))    

    X <- lapply(mk:(k-1),function(h.x)x[idx-h.x])
    X <- do.call("cbind",X)
    
    if(k==1) X <- matrix(X,ncol=1)
    
    colnames(X) <- paste0("X", ".", (mk+1):k-1,"/","m")
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

#' MIDAS lag structure with dates
#'
#' @param x a vector
#' @param k lags, a vector
#' @param datey low frequency dates
#' @param ... further arguments used in fitting MIDAS regression
#'
#' @return a matrix containing the first differences and the lag k+1.
#' @author Virmantas Kvedaras, Vaidotas Zemlys-Balevičius
#' @importFrom stats lag
#' @export
#'
#' @examples
#' x <- c(1:144)
#' y <- c(1:12)
#' datey <- (y-1)*12+1
#' 
#' #msld and mls should give the same results
#' 
#' m1 <- mlsd(x, 0:5, datey)
#' 
#' m2 <- mls(x, 0:5, 12)
#' 
#' sum(abs(m1 - m2))
#' 
mlsd <- function(x, k, datey, ... ) {
    
    datex <- NULL
    if (inherits(x,"ts")) datex <- time(x)
    if (inherits(x,"zoo") | inherits(x, "xts")) datex <- index(x)
    
    x <- as.numeric(x)
    if (is.null(datex)) datex <- 1:length(x)
    
    if (length(x) != length(datex)) stop("The date vector for high frequency data must be the same length as a data")
    
    ##We always assume that if we observe data for low frequency period 
    ##the high frequency data is observed until the end of that period. 
    ##We build the indexes of the high frequency observations. 
    ##For that we need to capture all the x observations, hence the minimum period
    ##must be the minimum x period. 
    ##The max period needs to be the next period after the maximum datey period.
    ##If there is any data which falls into first period, we can discard it.
    
    if (inherits(datey,"ts")) {
        datey0 <- datey
        left <- min(time(lag(datey, 1))) 
        right <- max(time(lag(datey, -1)))
        if (min(datex) < left) left <- min(datex)
        datey <- c(left, time(datey), right) - 0.001
        datey0 <- time(datey0) - 0.001
    } else {
        if (inherits(datey,"zoo") | inherits(datey, "xts")) {
            datey <- index(datey)
        }
        datey0 <- datey
        left <- datey[1] - (datey[2] - datey[1])
        if (min(datex) < left) left <- min(datex)
        nd <- length(datey)
        right <- datey[nd] + (datey[nd] - datey[nd - 1])
        datey <- c(left, datey, right)
    }
    

    ct <- cut(datex, datey, right = FALSE, labels = FALSE, include.lowest = TRUE)
    tct <- table(ct)
    uct <- unique(ct)
    nuct <- na.omit(uct)
    
    #We do not need the first period, but sometimes it is matched. 
    #In that case it is dropped.
    id <- match(2:(length(datey) - 1), nuct)
    
    #if (length(uct) != length(datey0)) {
    #    id <- id[-1]
    #}
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