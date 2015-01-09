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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param object 
##' @param nsim 
##' @param future 
##' @param newdata 
##' @param innov 
##' @param method 
##' @param insample 
##' @param seed 
##' @param show_progress 
##' @return 
##' @author Vaidotas Zemlys
simulate.midas_r <- function(object, nsim = 999, future=TRUE, newdata=NULL,
                             insample = NULL,
                             method = c("static", "dynamic"),
                             innov = NULL,                            
                             seed = NULL, show_progress = TRUE) {
    method <- match.arg(method)
    yname <- all.vars(object$terms[[2]])
    
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
    
    freqinfo <- get_frequency_info(object)    
    
    if(is.null(insample)) {
        insample <- get_estimation_sample(object)
    }    

    if(future) {
        outsample <- data_to_list(newdata)
        outsample <- outsample[intersect(names(outsample),names(freqinfo))]
        firstno <- names(outsample)[1]        
        h <- length(outsample[[firstno]])%/%freqinfo[firstno]
        if(is.null(innov)) innov <- matrix(sample(residuals(object), nsim*h, replace=TRUE), nrow=nsim)
        else {
            innov <- try(matrix(innov))
            if(inherits(class(innov),"try-error")) stop("Innovations must be a matrix or an object which can be coerced to matrix")
            if(ncol(innov)!=h) stop("The number of columns of innovation matrix must coincide with number of forecasting periods, which is ", h, " based on supplied new data")
            nsim <- ncol(innov)
        }        
        if(method == "static") {
            xm <- static_forecast(object, h, insample, outsample, yname)           
            sim <- xm + innov
        } else {
            sim <- matrix(NA, ncol = h, nrow = nsim)
            if(show_progress) {
                cat("\nDynamic simulation progress:\n")    
                pb <- txtProgressBar(min=0, max=nsim, initial=0, style=3)
            }
            for(j in 1:nsim) {                
                sim[j, ] <-  dynamic_forecast(object, h, insample[names(freqinfo)], outsample, freqinfo, innov[j, ])
                if(show_progress) setTxtProgressBar(pb, j)
            }
            if(show_progress) close(pb)
        }
    } else {
        modres <- residuals(object)
        nm <- length(modres)
        if(is.null(innov)) innov <- matrix(sample(modres, nsim*nm, replace=TRUE), ncol=nsim)
        else {
            innov <- try(matrix(innov))
            if(inherits(class(innov),"try-error")) stop("Innovations must be a matrix or an object which can be coerced to matrix")
            if(nrow(innov)!=length(residuals(object))) stop("The number of rows of innovation matrix must coincide with effective sample size, which is ", nm)
            nsim <- ncol(innov)
        }       
        if(method == "static") {
            sim <- fitted(object) + innov
        } else {
            sim <- matrix(NA, ncol = nm, nrow = nsim)
            n <- length(insample[[yname]])
            dsplit <- split_data(insample,1:(n-nm),(n-nm+1):n)
            if(show_progress) {
                cat("\nDynamic simulation progress:\n")    
                pb <- txtProgressBar(min=0, max=nsim, initial=0, style=3)
            }
            for(j in 1:nsim) {                
                sim[j, ] <-  dynamic_forecast(object, h, dsplit$indata, dsplit$outdata, freqinfo, innov[j, ])
                if(show_progress)setTxtProgressBar(pb, j)
            }
            if(show_progress) close(pb)
        }        
    }
    sim
}

dynamic_forecast <- function(object, h, fdata, outsample, freqinfo, innov = rep(0,h) ) {
    yname <- names(freqinfo)[1]
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
        res[i] <- rr[n] + innov[i]
        fdata[[yname]][n] <- res[i]             
    }
    res
}

static_forecast <- function(object, h, insample, outsample, yname) {
    if(!(yname %in% names(outsample))) {
        outsample <- c(list(rep(NA,h)),outsample)
        names(outsample)[1] <- yname
    }
    data <- try(rbind_list(insample[names(outsample)],outsample))
    if(class(data)=="try-error")stop("Missing variables in newdata. Please supply the data for all the variables (excluding the response variable) in regression")
    res <- predict.midas_r(object,newdata=data,na.action=na.pass)        
    n <- length(res)
    res[n+1-h:1]
}

