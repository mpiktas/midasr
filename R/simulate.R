##' Simulate AR(1) or MA(1) model
##'
##' Simulate MIDAS regressor as a AR(1) or MA(1) time series
##' @title 
##' @param model A named vector of length one. Name is either "ar", or "ma"
##' depending on which AR(1) or MA(1) process should be generated
##' @param n length of output series
##' @param innov.sd the standard error of innovations, which are zero mean normal random variables
##' @param frequency the frequency of the regressor, should larger than one.
##' @param n.start length of the burn.in period, the default is 300.
##' @return a time-series object of class "ts" 
##' @author Vaidotas Zemlys
simplearma.sim <- function(model,n,innov.sd,frequency,n.start=300) {
    e <- rnorm(n+n.start,sd=innov.sd)
    if(names(model)=="ar"){
        res <- filter(e,filter=model,method="recursive",sides=1)
    }
    if(names(model)=="ma"){
        res <- c(0,e[2:length(e)]+model*e[1:(length(e)-1)])
    }
    res[-n.start:-1]
}

##' Simulate MIDAS regresion response variable
##'
##' Given the predictor variable and the coefficients calculate MIDAS regression response variable.
##'
##' It is assumed that predictor variable has long enough history.
##' @title 
##' @param n 
##' @param theta 
##' @param x 
##' @param eps.sd 
##' @return 
##' @author Vaidotas Zemlys
midas.sim <- function(n,theta,x,eps.sd) {
    m <- frequency(x)
    n.x <- length(x)
    
    if(m==1) stop("The frequency of the predictor variable should be at least 2")
    if(n.x<=n*length(theta)) stop("The history of the predictor variable is not long enough, reduce the desired sample size")
   

    theta.d <- c(theta,rep(0,n.x-length(theta)))
    theta.d <- theta.d[n.x:1]

    y <- foreach(h.y = 1:n, .combine='c') %do% {
        t(x[1:(n.x-(n-h.y)*m)])%*%theta.d[((n-h.y)*m+1):n.x]
    }
    ts(y+rnorm(n,sd=eps.sd),frequency=1)    
}

##idx <- m*c((n.x/m-n+1):(n.x/m))
##X <- foreach(h.x=0:((kmax+1)*m-1), .combine='cbind')%do%{
##       x[idx-h.x]
##}
