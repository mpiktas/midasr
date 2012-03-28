##' Simulate AR(1) or MA(1) model
##'
##' Simulate MIDAS regressor as a AR(1) or MA(1) time series
##'  
##' @param model A named vector of length one. Name is either "ar", or "ma"
##' depending on which AR(1) or MA(1) process should be generated
##' @param n the length of output series
##' @param innov.sd the standard error of innovations, which are zero mean normal random variables
##' @param frequency the frequency of the regressor, should be larger than one.
##' @param n.start the length of the burn.in period, the default is 300.
##' @return a time-series object of class \code{ts} 
##' @author Vaidotas Zemlys 
##' @examples
##' #Generate AR(1) model with rho=0.6, with frequency 12
##' x <- simplearma.sim(list(ar=0.6),1500*12,1,12)
##' @export
simplearma.sim <- function(model,n,innov.sd,frequency,n.start=300) {
    e <- rnorm(n+n.start,sd=innov.sd)
    if(names(model)=="ar"){
        res <- filter(e,filter=model,method="recursive",sides=1)
    }
    if(names(model)=="ma"){
        res <- c(0,e[2:length(e)]+model*e[1:(length(e)-1)])
    }
    ts(res[-n.start:-1],frequency=frequency)
}

##' Simulate MIDAS regresion response variable
##'
##' Given the predictor variable and the coefficients calculate MIDAS regression response variable.
##' It is assumed that predictor variable has long enough history.
##' 
##' @param n The sample size
##' @param theta a vector with MIDAS regression coefficients 
##' @param x a \code{ts} object with MIDAS regression predictor variable
##' @param eps.sd the standard error of the regression disturbances, which are assumed to be independent normal zero mean random variables 
##' @return a \code{ts} object
##' @author Vaidotas Zemlys
##' @examples
##' ##The parameter function
##' theta.h0 <- function(p, dk) {
##'    i <- (1:dk-1)/100
##'    pol <- p[3]*i + p[4]*i^2
##'    (p[1] + p[2]*i)*exp(pol)
##' }
##'
##' ##Generate coefficients
##' theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)
##'
##' ##Plot the coefficients
##' plot(theta0)
##'
##' ##Generate the predictor variable
##' x <- simplearma.sim(list(ar=0.6),1500*12,1,12)
##'
##' ##Simulate the response variable
##' y <- midas.sim(500,theta0,x,1)
##' 
##' @export
##' @import foreach
midas.sim <- function(n,theta,x,eps.sd) {
    m <- frequency(x)
    n.x <- length(x)
    
    if(m==1) stop("The frequency of the predictor variable should be at least 2")
    if(n.x<=m*n+length(theta)-m) stop("The history of the predictor variable is not long enough, reduce the desired sample size")
   

    theta.d <- c(theta,rep(0,n.x-length(theta)))
    theta.d <- theta.d[n.x:1]

    y <- foreach(h.y = 1:n, .combine='c') %do% {
        t(x[1:(n.x-(n-h.y)*m)])%*%theta.d[((n-h.y)*m+1):n.x]
    }
    ts(y+rnorm(n,sd=eps.sd),frequency=1)    
}
