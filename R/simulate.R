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

##' Simulate simple MIDAS regression response variable
##'
##' Given the predictor variable and the coefficients simulate MIDAS regression response variable.
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
##' ##Generate the predictor variable, leave 4 low frequency lags of data for burn-in.
##' x <- ts(arima.sim(model = list(ar = 0.6), 504*12), frequency = 12)
##'
##' ##Simulate the response variable
##' y <- midas_sim(500, x, theta0)
##'
##' x <- window(xx, start=start(y))
##' midas_r(y ~ mls(y, 1, 1) + fmls(x, 4*12-1, 12, theta.h0), start = list(x = c(-0.1, 10, -10, -10)))
##' 
##' @details MIDAS regression with one predictor variable has the following form:
##' 
##' \deqn{y_t=\sum_{j=0}^{h}\theta_jx_{tm-j}+u_t,}
##' where \eqn{m} is the frequency ratio and
##' \eqn{h} is the number of high frequency lags included in the regression. 
##'
##' MIDAS regression involves times series with different frequencies. In R
##' the frequency property is set when creating time series objects
##' \code{\link{ts}}. Hence the frequency ratio \eqn{m} which figures in MIDAS regression is calculated from frequency property of time series objects supplied.
##' @export
midas_sim <- function(n, x, theta, rand_gen = rnorm, innov = rand_gen(n, ...), ...) {
    m <- frequency(x)
    n_x <- length(x)
      
    if(n_x<=m*n+length(theta)-m) stop("The history of the predictor variable is not long enough, reduce the desired sample size")
   
    X <- fmls(x,length(theta)-1,m)
    xt <- as.vector(X%*%theta)

    ts(xt[nrow(X)-n:1+1],start=end(x)[1]-n+1,frequency=1) + innov    
}
##' Simulate simple autoregressive MIDAS model
##'
##' Given the predictor variable, the weights and autoregressive coefficients, simulate MIDAS regression response variable.
##' 
##' @param n sample size.
##' @param alpha autoregressive coefficients.
##' @param x a high frequency predictor variable.
##' @param theta a vector with MIDAS weights for predictor variable.
##' @param rand_gen a function to generate the innovations, default is the normal distribution.
##' @param innov an optional time series of innovations.
##' @param n_start number of observations to omit for the burn.in.
##' @param ... 
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
##' xx <- ts(arima.sim(model = list(ar = 0.6), 3000 *12), frequency = 12)
##' 
##' y <- midas_auto_sim(500, 0.5, xx, theta0, n_start = 200)
##' x <- window(xx, start=start(y))
##' midas_r(y ~ mls(y, 1, 1) + fmls(x, 4*12-1, 12, theta.h0), start = list(x = c(-0.1, 10, -10, -10)))
midas_auto_sim <- function(n, alpha, x, theta, rand_gen = rnorm, innov = rand_gen(n, ...), n_start = NA, ...) {
    m <- frequency(x)
    n_x <- length(x)

    minroots <- min(Mod(polyroot(c(1, -alpha))))
    if (minroots <= 1) 
      stop("'ar' part of model is not stationary")

    if(is.na(n_start))n_start <- 100

    nout <- n
    n <- n+n_start
                                         
    if(n_x<=m*n+length(theta)-m) stop("The history of the predictor variable is not long enough, reduce the desired sample size")
                                  

    X <- fmls(x,length(theta)-1,m)
    xt <- as.vector(X%*%theta)

    ##innov is called here, so its argument n already includes the n_start
    xte <- ts(xt[nrow(X)-n:1+1],start=end(x)[1]-n+1,frequency=1) + innov

    y <- filter(xte,alpha,method="recursive")
    ts(y[-(1:n_start)],start=end(x)[1]-nout+1,frequency=1)
}


##simulate.midas_r
##1. If there is no autoregression, simply add future xreg times coefficient + bootstrapped error
##2. If there is autoregression, run dynamic forecast with future xreg, past y and bootstrapped error.

##Call simulate from forecast.midas_r. 

##Maybe rewrite simulate based on formula?
##Do a bsic 

simulate.midas_r <- function(object, nsim = nrow(object$model), future=TRUE, newdata=NULL,
                             innov = NULL, bootstrap = FALSE, type = c("static", "dynamic"),
                             freqinfo = NULL, insample = NULL) {
    if (is.null(innov)) {
        if (!exists(".Random.seed", envir = .GlobalEnv)) 
            runif(1)
        if (is.null(seed)) 
            RNGstate <- .Random.seed
        else {
            R.seed <- .Random.seed
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
    }
    if(is.null(freqinfo)) {
        freqinfo <- get_frequency_info(object$terms, object$Zenv)
    }
    
    if(is.null(insample)) {
        insample <- get_estimation_sample(object)
    }
    
    if(future) {
        e <- tsboot(residuals(object), h) 
        if(type == "static") {
            res <- predict.midas_r(object, newdata = data, na.action = na.pass)
            n <- length(res)
            xm <- res[n+1-h:1]
            sim <- xm + e
        } else {
            fdata <- insample[names(freqinfo)]
            outsample <- outsample[names(outsample)!=yname]
            yna <- list(NA)
            names(yna) <- yname
            res <- rep(NA,h)        
            for(i in 1:h) {
                ##Get the data for the current low frequency lag
                hout <- mapply(function(var,m){
                    var[1:m+(i-1)*m]
                },outsample,freqinfo[names(outsample)],SIMPLIFY=FALSE)
                hout <- c(yna,hout)
                fdata <- rbind_list(fdata[names(hout)],hout)
                if(class(fdata)=="try-error")stop("Missing variables in newdata. Please supply the data for all the variables (excluding the response variable) in regression")
                rr <- predict.midas_r(object,newdata=fdata,na.action=na.pass)
                n <- length(rr)
                res[i] <- rr[n] + e[i]
                fdata[[yname]][n] <- res[i]             
            }
            sim <- res
        }
    } else {
        e <- tsboot(residuals(object))
        if(type == "static") {
            sim <- fitted(object) + e
        } else {
            
        }        
    }
    sim
}
