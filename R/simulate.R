##' Simulate simple MIDAS regression response variable
##'
##' Given the predictor variable and the coefficients simulate MIDAS regression response variable.
##' 
##' @param n The sample size
##' @param x a \code{ts} object with MIDAS regression predictor variable
##' @param theta a vector with MIDAS regression coefficients 
##' @param rand_gen the function which generates the sample of innovations, the default is \code{\link{rnorm}} 
##' @param innov the vector with innovations, the default is NULL, i.e. innovations are generated using argument \code{rand_gen}
##' @param ... additional arguments to \code{rand_gen}.
##' @return a \code{ts} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @examples
##' ##The parameter function
##' theta_h0 <- function(p, dk) {
##'    i <- (1:dk-1)/100
##'    pol <- p[3]*i + p[4]*i^2
##'    (p[1] + p[2]*i)*exp(pol)
##' }
##'
##' ##Generate coefficients
##' theta0 <- theta_h0(c(-0.1,10,-10,-10),4*12)
##'
##' ##Plot the coefficients
##' plot(theta0)
##'
##' ##Generate the predictor variable, leave 4 low frequency lags of data for burn-in.
##' xx <- ts(arima.sim(model = list(ar = 0.6), 600 * 12), frequency = 12)
##'
##' ##Simulate the response variable
##' y <- midas_sim(500, xx, theta0)
##'
##' x <- window(xx, start=start(y))
##' midas_r(y ~ mls(y, 1, 1) + fmls(x, 4*12-1, 12, theta_h0), start = list(x = c(-0.1, 10, -10, -10)))
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
##' @param ... additional arguments to function \code{rand_gen}.
##' @return a \code{ts} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @importFrom stats filter frequency rnorm
##' @export
##' @examples
##' theta_h0 <- function(p, dk) {
##'   i <- (1:dk-1)/100
##'   pol <- p[3]*i + p[4]*i^2
##'   (p[1] + p[2]*i)*exp(pol)
##' }
##' 
##' ##Generate coefficients
##' theta0 <- theta_h0(c(-0.1,10,-10,-10),4*12)
##'
##' ##Generate the predictor variable
##' xx <- ts(arima.sim(model = list(ar = 0.6), 1000 * 12), frequency = 12)
##' 
##' y <- midas_auto_sim(500, 0.5, xx, theta0, n_start = 200)
##' x <- window(xx, start=start(y))
##' midas_r(y ~ mls(y, 1, 1) + fmls(x, 4*12-1, 12, theta_h0), start = list(x = c(-0.1, 10, -10, -10)))
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

##' Simulates one or more responses from the distribution corresponding to a fitted MIDAS regression object.
##'
##' Only the regression innovations are simulated, it is assumed that the predictor variables and coefficients are fixed. The
##' innovation distribution is simulated via bootstrap. 
##' @title Simulate MIDAS regression response
##' @param object \code{\link{midas_r}} object
##' @param nsim number of simulations
##' @param seed either NULL or an integer that will be used in a call to set.seed before simulating the time series. The default, NULL will not change the random generator state.
##' @param future logical, if \code{TRUE} forecasts are simulated, if \code{FALSE} in-sample simulation is performed.
##' @param newdata a named list containing future values of mixed frequency regressors.  The default is \code{NULL}, meaning that only in-sample data is used.
##' @param insample a list containing the historic mixed frequency data 
##' @param method the simulation method, if \code{"static"} in-sample values for dependent variable are used in autoregressive MIDAS model, if \code{"dynamic"}
##' the dependent variable values are calculated step-by-step from the initial in-sample values.
##' @param innov a matrix containing the simulated innovations. The default is \code{NULL}, meaning, that innovations are simulated from model residuals.
##' @param show_progress logical, TRUE to show progress bar, FALSE for silent evaluation
##' @param ... not used currently
##' @return a matrix of simulated responses. Each row contains a simulated response.
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @method simulate midas_r
##' @rdname simulate.midas_r
##' @examples
##' data("USrealgdp")
##' data("USunempr")
##' 
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr), start = 1949)
##' trend <- 1:length(y)
##' 
##' ##24 high frequency lags of x included
##' mr <- midas_r(y ~ trend + fmls(x, 23, 12, nealmon), start = list(x = rep(0, 3)))
##'
##' simulate(mr, nsim=10, future=FALSE)
##'
##' ##Forecast horizon
##' h <- 3
##' ##Declining unemployment
##' xn <- rep(-0.1, 12*3)
##' ##New trend values
##' trendn <- length(y) + 1:h
##'
##' simulate(mr, nsim = 10, future = TRUE, newdata = list(trend = trendn, x = xn))
##' @importFrom stats runif
##' @export
simulate.midas_r <- function(object, nsim = 999, seed = NULL, future=TRUE, newdata=NULL,
                             insample = NULL,
                             method = c("static", "dynamic"),
                             innov = NULL,                            
                             show_progress = TRUE, ...) {
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

##' @importFrom stats simulate
##' @name simulate
##' @rdname simulate.midas_r
##' @export
NULL

#' Simulate LSTR MIDAS regression model
#'
#' @param n number of observations to simulate.
#' @param m integer, frequency ratio
#' @param theta vector, restriction coefficients for high frequency variable
#' @param intercept vector of length 1, intercept for the model.
#' @param plstr vector of length 4, slope for the LSTR term and LSTR parameters 
#' @param ar.x vector, AR parameters for simulating high frequency variable
#' @param ar.y vector, AR parameters for AR part of the model
#' @param rand.gen function, a function for generating the regression innovations, default is \code{rnorm}
#' @param n.start integer, length of a 'burn-in' period. If NA, the default, a reasonable value is computed.
#' @param ... additional parameters to rand.gen
#'
#' @return a list
#' @export
#'
#' @examples
#' 
#' nnbeta <- function(p, k) nbeta(c(1,p),k)
#' 
#' dgp <- midas_lstr_sim(250, m = 12, theta = nnbeta(c(2, 4), 24), 
#'                            intercept = c(1), plstr = c(1.5, 1, log(1), 1), 
#'                            ar.x = 0.9, ar.y = 0.5, n.start = 100)

#' z <- cbind(1, mls(dgp$y, 1:2, 1))
#' colnames(z) <- c("Intercept", "y1", "y2")
#' X <- mls(dgp$x, 0:23, 12)
#'
#' lstr_mod <- midas_lstr_plain(dgp$y, X, z, nnbeta, 
#'                           start_lstr = c(1.5, 1, log(1), 1), 
#'                           start_x = c(2, 4), start_z=c(1, 0.5, 0)) 
#' 
#' coef(lstr_mod)
#' 
#' @importFrom stats filter sd
midas_lstr_sim <- function(n, m, theta, intercept, plstr, ar.x,  ar.y, 
                           rand.gen = rnorm,  n.start = NA, ...) {
    
    minroots <- min(Mod(polyroot(c(1, -ar.y))))
    
    if (minroots <= 1) stop("'ar' part of model is not stationary")
    if(is.na(n.start)) n.start <- length(ar.y) + ceiling(6/log(minroots))
    
    innov_x <- rand.gen(m*(n + n.start))
    
    x <- filter(innov_x, ar.x, method = 'recursive', init = 0)
    
    xx <- mls(x, 0:(length(theta) - 1), m)
    
    sd_x <- sd(c(xx), na.rm = TRUE)
     
    g <- intercept + lstr(xx, theta, plstr, sd_x)
    
    g[is.na(g)] <- 0
    
    innov <- rand.gen(length(g), ...)
    
    y <- filter(g + innov, ar.y, method = "recursive", init = 0)
    
    
    y <- y[-seq_len(n.start)]
    x <- x[-seq_len(n.start*m)]

    list(y = y, x = x, lstr = lstr, intercept = intercept, ar.y = ar.y)
}


#' Simulate MMM MIDAS regression model
#'
#' @param n number of observations to simulate.
#' @param m integer, frequency ratio
#' @param theta vector, restriction coefficients for high frequency variable
#' @param intercept vector of length 1, intercept for the model.
#' @param pmmm vector of length 2, slope for the MMM term and MMM parameter
#' @param ar.x vector, AR parameters for simulating high frequency variable
#' @param ar.y vector, AR parameters for AR part of the model
#' @param rand.gen function, a function for generating the regression innovations, default is \code{rnorm}
#' @param n.start integer, length of a 'burn-in' period. If NA, the default, a reasonable value is computed.
#' @param ... additional parameters to rand.gen
#'
#' @return a list
#' @export
#'
#' @examples
#' 
#' nnbeta <- function(p, k) nbeta(c(1,p),k)
#' 
#' dgp <- midas_mmm_sim(250, m = 12, theta = nnbeta(c(2, 4), 24), 
#'                            intercept = c(1), pmmm = c(1.5, 1), 
#'                            ar.x = 0.9, ar.y = 0.5, n.start = 100)

#' z <- cbind(1, mls(dgp$y, 1:2, 1))
#' colnames(z) <- c("Intercept", "y1", "y2")
#' X <- mls(dgp$x, 0:23, 12)
#'
#' mmm_mod <- midas_mmm_plain(dgp$y, X, z, nnbeta, 
#'                           start_mmm = c(1.5, 1), 
#'                           start_x = c(2, 4), start_z=c(1, 0.5, 0)) 
#' 
#' coef(mmm_mod)
#' 
#' @importFrom stats filter sd
midas_mmm_sim <- function(n, m, theta, intercept, pmmm, ar.x,  ar.y, 
                           rand.gen = rnorm,  n.start = NA, ...) {
    
    minroots <- min(Mod(polyroot(c(1, -ar.y))))
    
    if (minroots <= 1) stop("'ar' part of model is not stationary")
    if(is.na(n.start)) n.start <- length(ar.y) + ceiling(6/log(minroots))
    
    innov_x <- rand.gen(m*(n + n.start))
    
    x <- filter(innov_x, ar.x, method = 'recursive', init = 0)
    
    xx <- mls(x, 0:(length(theta) - 1), m)
    
    g <- intercept + mmm(xx, theta, pmmm)
    
    g[is.na(g)] <- 0
    
    innov <- rand.gen(length(g), ...)
    
    y <- filter(g + innov, ar.y, method = "recursive", init = 0)
    
    
    y <- y[-seq_len(n.start)]
    x <- x[-seq_len(n.start*m)]
    
    list(y = y, x = x, mmm = mmm, intercept = intercept, ar.y = ar.y)
}
