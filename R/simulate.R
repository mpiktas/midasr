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
##' @author Virmantas Kvedaras, Vaidotas Zemlys 
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
##' 
##' @param n The sample size
##' @param theta a vector with MIDAS regression coefficients 
##' @param x a \code{ts} object with MIDAS regression predictor variable
##' @param eps.sd the standard error of the regression disturbances, which are assumed to be independent normal zero mean random variables 
##' @return a \code{ts} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
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
##' @details MIDAS regression has the following form:
##' 
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+u_t}
##'
##' or alternatively
##'
##' \deqn{y_t=\sum_{h=0}^{(k+1)m}\theta_hx_{tm-h}+u_t,}
##' where \eqn{m} is the frequency ratio and
##' \eqn{k} is the number of lags included in the regression. 
##'
##' MIDAS regression involves times series with different frequencies. In R
##' the frequency property is set when creating time series objects
##' \code{\link{ts}}. Hence the frequency ratio \eqn{m} which figures in MIDAS regression is calculated from frequency property of time series objects supplied.
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
    
    ts(y+rnorm(n,sd=eps.sd),start=end(x)[1]-n+1,frequency=1)    
}
##' Simulate autoregressive MIDAS model
##'
##' Given the predictor variable, the weights and autoregressive coefficients calculate MIDAS regression response variable.
##' 
##' @param n sample size
##' @param theta a vector with MIDAS weights for predictor variable
##' @param alpha autoregressive coefficients
##' @param x a high frequency predictor variable
##' @param eps.sd the standard error of the regression disturbances, which are assumed to be independent normal zero mean random variables
##' @param n.start number of observations to ommit for the burn.in
##' @return a \code{ts} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
##' @examples
##' theta.h0 <- function(p, dk) {
##'   i <- (1:dk-1)/100
##'   pol <- p[3]*i + p[4]*i^2
##'   (p[1] + p[2]*i)*exp(pol)
##' }
##' 
##' ##Generate coefficients
##' theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)
##'
##' ##Generate the predictor variable
##' xx <- simplearma.sim(list(ar=0.6),3000*12,1,12)
##'
##' y <- midas.auto.sim(500,theta0,c(0.5),xx,1,n.start=100)
##' x <- window(xx,start=start(y))
##' midas_r(y~mls(y,1,1)+fmls(x,4*12-1,12,theta.h0),start=list(x=c(-0.1,10,-10,-10)))
midas.auto.sim <- function(n,theta,alpha,x,eps.sd,n.start=NA) {
    m <- frequency(x)
    n.x <- length(x)

    minroots <- min(Mod(polyroot(c(1, -alpha))))
    if (minroots <= 1) 
      stop("'ar' part of model is not stationary")

    if(is.na(n.start))n.start <- 100

    nout <- n
    n <- n+n.start
                                         
    if(m==1) stop("The frequency of the predictor variable should be at least 2")
    if(n.x<=m*n+length(theta)-m) stop("The history of the predictor variable is not long enough, reduce the desired sample size")

                                  

    X <- fmls(x,length(theta)-1,m)
    xt <- as.vector(X%*%theta)
    

    xte <- ts(xt[nrow(X)-n:1+1],start=end(x)[1]-n+1,frequency=1)+rnorm(n,sd=eps.sd)

    y <- filter(xte,alpha,method="recursive")
    ts(y[-(1:n.start)],start=end(x)[1]-nout+1,frequency=1)
}
