##' MIDAS regression model deviance
##'
##' Returns the deviance of a fitted MIDAS regression object
##' @param object a \code{\link{midas_r}} object
##' @param ... currently nothing
##' @return The sum of squared residuals
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname deviance.midas_r
##' @method deviance midas_r
##' @export
deviance.midas_r <- function(object,...) {
    sum(residuals(object)^2,na.rm=TRUE)
}

##' Predicted values based on \code{midas_r} object.
##'
##' \code{predict.midas_r} produces predicted values, obtained by evaluating regression function in the frame \code{newdata}. This means that the appropriate model matrix is constructed using only the data in \code{newdata}. This makes this function not very convenient for forecasting purposes. If you want to supply the new data for forecasting horizon only use the function \link{forecast.midas_r}. Also this function produces only static predictions, if you want dynamic forecasts use the \link{forecast.midas_r}.
##' 
##' @title Predict method for MIDAS regression fit
##' @param object \code{\link{midas_r}} object
##' @param newdata a named list containing data for mixed frequencies. If omitted, the in-sample values are used.
##' @param na.action function determining what should be done with missing values in \code{newdata}. The most likely cause of missing values is the insufficient data for the lagged variables. The default is to omit such missing values.
##' @param ... additional arguments, not used
##' @return a vector of predicted values
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @method predict midas_r
##' @rdname predict.midas_r
##' @examples
##' data("USrealgdp")
##' data("USunempr")
##' 
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr), start = 1949)
##' 
##' ##24 high frequency lags of x included
##' mr <- midas_r(y ~ fmls(x, 23, 12, nealmon), start = list(x = rep(0, 3)))
##' 
##' ##Declining unemployment
##' xn <- rnorm(2 * 12, -0.1, 0.1)
##' 
##' ##Only one predicted value, historical values discarded
##' predict(mr, list(x = xn))
##' 
##' ##Historical values taken into account
##' forecast(mr, list(x = xn))
##'
##' @export
predict.midas_r <- function(object, newdata, na.action = na.omit, ... ) {
    Zenv <- new.env(parent=parent.frame())
    
    if(missing(newdata))
      return(as.vector(fitted(object)))
    else {
        ee <- data_to_env(newdata)    
        ZZ <- object$Zenv    
        parent.env(ee) <- ZZ
    }
    
    assign("ee",ee,Zenv)
    cll <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("object", "newdata"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    Terms <- delete.response(terms(object))
    mf[[2L]] <- Terms    
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- na.action
    names(mf)[c(2,3,4)] <- c("formula","data","na.action")
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")

    X <- model.matrix(mt, mf) 
    as.vector(X %*% coef(object, midas = TRUE))
}
##' @export
##' @method predict midas_u
predict.midas_u <- predict.midas_r

##' @export
##' @method summary midas_r
summary.midas_r <- function(object, vcov.=vcovHAC, df=NULL, prewhite=TRUE, ...) {
    r <- as.vector(residuals(object))
    param <- coef(object)
    pnames <- names(param)
    n <- length(r)
    p <- length(param)
    rdf <- n - p
    resvar <- deviance(object)/rdf
    XD <- object$model[,-1]%*%object$gradD(coef(object))
    R <- qr.R(qr(XD))
    XDtXDinv <- chol2inv(R)
    dimnames(XDtXDinv) <- list(pnames,pnames)
        
    se <- sqrt(diag(XDtXDinv)*resvar)

    f <- as.vector(object$fitted.values)
    mss <- if (attr(object$terms, "intercept")) 
        sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)

    n <- length(r)
    p <- length(coef(object))
    rdf <- n-p
    df.int <- if (attr(object$terms, "intercept")) 1L
    else 0L

    r_squared <- mss/(mss + rss)
    adj_r_squared <- 1 - (1 - r_squared) * ((n - df.int)/rdf)
    
    if(!is.null(vcov.)) {
        set <- try(sqrt(diag(vcov.(object,prewhite=prewhite,...))))
        if(class(set)=="try-error") {
            warning("Unable to compute robust standard errors, using non-robust ones. This is an indication of problems with optimisation, please try other starting values or change optimisation method")
        } else se <- set
    }
    tval <- param/se

    #Code stolen from coeftest function
    if(is.null(df)) {
        df <- rdf
    }
    if (is.finite(df) && df > 0) {
        pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
        cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
        mthd <- "t"
    }
    else {
        pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
        cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
        mthd <- "z"
    }
    
    param <- cbind(param,se,tval,pval)
    dimnames(param) <- list(pnames, c("Estimate", "Std. Error", 
        "t value", "Pr(>|t|)"))
    ans <- list(formula=formula(object$terms),residuals=r,sigma=sqrt(resvar),
                df=c(p,rdf), cov.unscaled=XDtXDinv, call=object$call,
                coefficients=param,midas_coefficients=coef(object, midas = TRUE),
                r_squared = r_squared, adj_r_squared = adj_r_squared)
    class(ans) <- "summary.midas_r"
    ans
}

##' @export
##' @method print summary.midas_r
print.summary.midas_r <- function(x, digits=max(3, getOption("digits") - 3 ), signif.stars = getOption("show.signif.stars"), ...) {
    cat("\n Formula", deparse(formula(x)),"\n")
    df <- x$df
    rdf <- df[2L]
    cat("\n Parameters:\n")
    printCoefmat(x$coefficients,digits=digits,signif.stars=signif.stars,...)
    cat("\n Residual standard error:", format(signif(x$sigma,digits)), "on", rdf , "degrees of freedom\n")
 #   cat(" Multiple R-squared: ", formatC(x$r_squared, digits = digits))
 #   cat(",\tAdjusted R-squared: ", formatC(x$adj_r_squared, digits = digits),"\n")
    invisible(x)
}

##' @export
##' @method print midas_r
print.midas_r <- function(x, digits=max(3,getOption("digits")-3),...) {
    cat("MIDAS regression model\n")
    cat(" model:", deparse(formula(x)),"\n")
    print(coef(x),digits = digits, ...)
    cat("\n")
    cat("Function", x$argmap_opt$Ofunction, "was used for fitting\n")
    invisible(x)
}

##' @import sandwich
##' @export
##' @method estfun midas_r
estfun.midas_r <- function(x,...) {
    XD <- x$model[,-1]%*%x$gradD(coef(x))
    rval <- as.vector(residuals(x))*XD
    colnames(rval) <- names(coef(x))
    rval
}

##' @export
##' @method vcov midas_r
vcov.midas_r <- function(x,...) {
    sm <- summary(object)
    sm$cov.unscaled * sm$sigma^2
}

##' @export
##' @method bread midas_r
bread.midas_r <- function(x,...) {
    sx <- summary(x, vcov.=NULL)
    return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])))
}

##' @export
##' @method vcov midas_r
vcov.midas_r <- function(object,...) {
    ##Code stolen from stats:::vcov.nls
    sm <- summary(object)
    sm$cov.unscaled * sm$sigma^2
}


##' @export
##' @method logLik midas_r
logLik.midas_r <- function(object,...) {
    res <- residuals(object)
    N <- length(res)
    val <- -N * (log(2 * pi) + 1 - log(N)  +  log(sum(res^2)))/2
    attr(val, "df") <- 1L + length(coef(object))
    attr(val, "nobs") <- attr(val, "nall") <- N
    class(val) <- "logLik"
    val
}

##' Extracts various coefficients of MIDAS regression
##'
##' MIDAS regression has two sets of cofficients. The first set is the coefficients associated with the parameters
##' of weight functions associated with MIDAS regression terms. These are the coefficients of the NLS problem associated with MIDAS regression.
##' The second is the coefficients of the linear model, i.e  the values of weight
##' functions of terms, or so called MIDAS coefficients. By default the function returns the first set of the coefficients.
##' 
##' @title Extract coefficients of MIDAS regression
##' @param object \code{midas_r} object
##' @param midas logical, if \code{TRUE}, MIDAS coefficients are returned, if \code{FALSE} (default), coefficients of NLS problem are returned
##' @param term_names a character vector with term names. Default is \code{NULL}, which means that coefficients of all the terms are returned
##' @param ... not used currently
##' @return a vector with coefficients
##' @author Vaidotas Zemlys
##' @method coef midas_r
##' @rdname coef.midas_r
##' @examples
##'
##' #Simulate MIDAS regression
##' n<-250
##' trend<-c(1:n)
##' x<-rnorm(4*n)
##' z<-rnorm(12*n)
##' fn.x <- nealmon(p=c(1,-0.5),d=8)
##' fn.z <- nealmon(p=c(2,0.5,-0.1),d=17)
##' y<-2+0.1*trend+mls(x,0:7,4)%*%fn.x+mls(z,0:16,12)%*%fn.z+rnorm(n)
##' eqr<-midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) +
##'              mls(z, 0:16, 12, nealmon),
##'              start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
##'
##' coef(eqr)
##' coef(eqr, term_names = "x")
##' coef(eqr, midas = TRUE)
##' coef(eqr, midas = TRUE, term_names = "x")
##' 
##' @export
coef.midas_r <- function(object, midas = FALSE, term_names = NULL, ...) {
    if(is.null(term_names)) {
        if(midas) return(object$midas_coefficients)
        else return(object$coefficients)
    } else {
        if(length(setdiff(term_names,names(object$term_info)))>0) stop("Some of the term names are not present in estimated MIDAS regression")
        if(midas) {
            res <- lapply(object$term_info[term_names], function(x) object$midas_coefficients[x$midas_coef_index])
        } else {
            res <- lapply(object$term_info[term_names], function(x) object$coefficients[x$coef_index])
        }
        names(res) <- NULL
        if(length(res)==1) return(res[[1]])
        else return(unlist(res))
    }        
}

get_frequency_info <- function(x,...) UseMethod("get_frequency_info")

get_frequency_info.midas_r <- function(object) {    
    yname <- all.vars(object$terms[[2]])
    res <- sapply(object$term_info,"[[","frequency")
    ##Make sure response variable is first
    if(!(yname %in% names(res))) {
        res <- c(1,res)
        names(res)[1] <- yname
    } else {
        res <- res[c(yname,setdiff(names(res),yname))]
    }
    res[setdiff(names(res),"(Intercept)")]
}

get_frequency_info.default <- function (mt, Zenv) 
{
    vars <- as.list(attr(mt, "variables"))[-1]
    res <- lapply(vars, function(l) {
        if (length(l) > 1) {
            if (as.character(l[[1]]) %in% c("mls", "fmls", "dmls")) {
                m <- eval(l[[4]], Zenv)
                varnames <- as.character(all.vars(l[[2]]))
                list(varname = varnames, m = rep(m, length(varnames)))
            }
            else {
                varnames <- as.character(all.vars(l))
                list(varname = varnames, m = rep(1, length(varnames)))
            }
        }
        else list(varname = as.character(l), m = 1)
    })
    varn <- Reduce("c", lapply(res, "[[", "varname"))
    freq <- Reduce("c", lapply(res, "[[", "m"))
    out <- freq
    names(out) <- varn
    out[unique(names(out))]
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


##' Forecasts MIDAS regression given the future values of regressors. For dynamic models (with lagged response variable) there is an option to calculate dynamic forecast, when forecasted values of response variable are substituted into the lags of response variable.
##' 
##' Given future values of regressors this function combines the historical values used in the fitting the MIDAS regression model and calculates the forecasts.
##' 
##' @title Forecast MIDAS regression
##' @param object midas_r object
##' @param newdata a named list containing future values of mixed frequency regressors. The default is \code{NULL}, meaning that only in-sample data is used.
##' @param se logical, if \code{TRUE}, the prediction intervals are calculated
##' @param level confidence level for prediction intervals
##' @param fan if TRUE, level is set to seq(50,99,by=1). This is suitable for fan plots
##' @param npaths the number of samples for simulating prediction intervals
##' @param method the forecasting method, either \code{"static"} or \code{"dynamic"}
##' @param insample a list containing the historic mixed frequency data 
##' @param show_progress logical, if \code{TRUE}, the progress bar is shown if \code{se = TRUE}
##' @param add_ts_info logical, if \code{TRUE}, the forecast is cast as \code{ts} object. Some attempts are made to guess the correct start, by assuming that the response variable is a \code{ts} object of \code{frequency} 1. If \code{FALSE}, then the result is simply a numeric vector.
##' @param ... additional arguments to \code{simulate.midas_r}
##' @return an object of class \code{"forecast"}, a list containing following elements:
##'
##' \item{method}{the name of forecasting method: MIDAS regression, static or dynamic}
##' \item{model}{original object of class \code{midas_r}}
##' \item{mean}{point forecasts}
##' \item{lower}{lower limits for prediction intervals}
##' \item{upper}{upper limits for prediction intervals}
##' \item{fitted}{fitted values, one-step forecasts}
##' \item{residuals}{residuals from the fitted model}
##' \item{x}{the original response variable}
##'
##' The methods \code{print}, \code{summary} and \code{plot} from package \code{forecast} can be used on the object.
##' @author Vaidotas Zemlys
##' @method forecast midas_r
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
##' ##Forecast horizon
##' h <- 3
##' ##Declining unemployment
##' xn <- rep(-0.1, 12*3)
##' ##New trend values
##' trendn <- length(y) + 1:h
##' 
##' ##Static forecasts combining historic and new high frequency data
##' forecast(mr, list(trend = trendn, x = xn), method = "static")
##' 
##' ##Dynamic AR* model
##' mr.dyn <- midas_r(y ~ trend + mls(y, 1:2, 1, "*")
##'                    + fmls(x, 11, 12, nealmon),
##'                   start = list(x = rep(0, 3)))
##' 
##' forecast(mr.dyn, list(trend = trendn, x = xn), method = "dynamic")
##'
##' ##Use print, summary and plot methods from package forecast
##' 
##' fmr <- forecast(mr, list(trend = trendn, x = xn), method = "static")
##' fmr
##' summary(fmr)
##' plot(fmr)
##' 
##' @export 
forecast.midas_r <- function(object, newdata=NULL, se = FALSE, level=c(80,95),
                             fan=FALSE, npaths=999,
                             method=c("static","dynamic"), insample=get_estimation_sample(object),
                             show_progress = TRUE, add_ts_info = FALSE, ...) {
    method <- match.arg(method)
    
    pred <- point_forecast.midas_r(object, newdata = newdata, method = method, insample = insample)
    if(se) {
        sim <- simulate(object, nsim = npaths, future = TRUE, newdata = newdata, method = method, insample = insample, show_progress = show_progress, ...)
        if (fan) 
            level <- seq(51, 99, by = 3)
        else {
            if (min(level) > 0 & max(level) < 1) 
                level <- 100 * level
            else if (min(level) < 0 | max(level) > 99.99) 
                stop("Confidence limit out of range")
        }
        nint <- length(level)        
        lower <- apply(sim, 2, quantile, 0.5 - level/200, type = 8)
        upper <- apply(sim, 2, quantile, 0.5 + level/200, type = 8)
        if (nint > 1L) {
            lower <- t(lower)
            upper <- t(upper)
        }    
    } else {
        lower <- NULL
        upper <- NULL
    }
    xout <- object$model[, 1]
    if(add_ts_info) {
        xstart <- 1
        if(!is.null(rownames(object$model))) {
            xstart <- as.numeric(rownames(object$model)[1])
            if(is.na(xstart)) xstart <- 1
            }
        xout <- ts(xout, start = xstart, frequency = 1)
        pred <- ts(pred, start = end(xout)+1, frequency = 1)
        }
    return(structure(list(method = paste0("MIDAS regression forecast (",method,")"),
                          model = object,
                          level = level,
                          mean = pred,
                          lower = lower,
                          upper = upper,
                          fitted = predict(object),
                          residuals = residuals(object),                          
                          x = xout
                          ), class = "forecast"))
}

point_forecast.midas_r <- function(object, newdata=NULL, method=c("static","dynamic"), insample=get_estimation_sample(object), ...) {

    method <- match.arg(method)

    outsample <- data_to_list(newdata)
    ##Fix this, as for default it is not necessary
    insample <- data_to_list(insample)
    ##Get the name of response variable
    yname <- all.vars(object$terms[[2]])
   
    ##Get the frequency information of each variable    
    freqinfo <- get_frequency_info(object)

    if(length(setdiff(names(freqinfo),c(yname,names(outsample))))>0) stop("Missing variables in newdata. Please supply the data for all the variables (excluding the response variable) in regression")

    if(length(setdiff(names(freqinfo),c(yname,names(insample))))>0) stop("Missing variables in in-sample. Please supply the data for all the variables in regression")
    
    outsample <- outsample[intersect(names(outsample),names(freqinfo))]
    firstno <- names(outsample)[1]
 
    h <- length(outsample[[firstno]])%/%freqinfo[firstno]
        
    if(method=="static") {
        res <- static_forecast(object, h, insample, outsample, yname)
    }
    else {
        res <- dynamic_forecast(object, h, insample[names(freqinfo)], outsample, freqinfo)        
    }
    res
}

##' @export
##' @method forecast midas_u
forecast.midas_u <- forecast.midas_r

rbind_list <- function(el1,el2) {
    if(is.null(el1)) return(el2)
    if(is.null(el2)) return(el1)
    
    if(!identical(names(el1),names(el2)))stop("You can rbind only lists with identical names")

    nms <- names(el1)
    l <- list(el1,el2)
    out <- lapply(nms,function(nm)do.call("c",lapply(l,function(x)x[[nm]])))
    names(out) <- nms
    out
}

data_to_env <- function(data) {
    as.environment(data_to_list(data))
}

data_to_list <- function(data) {
    if(is.matrix(data)) data <- data.frame(data)
    if(is.data.frame(data)) {
        ee <- as.list(data)
    }
    else {
        if(is.list(data)) {
            if(is.null(names(data))) names(data) <- rep("",length(data))
            data <- mapply(function(x,nm){
                if(is.null(dim(x))) {
                    x <- list(x)
                    names(x) <- nm
                    x
                } else {
                    ##This is needed since if tseries library is not loaded as.list for mts does not work as expected                   
                    if(inherits(x,"mts")) 
                        x <- data.frame(x)
                    if(ncol(x)==1) {
                        if(!is.null(colnames(x))) {
                            if(nm=="") nm <- colnames(x)
                            else warning("Duplicate names in data. Using the one from the list")                                                                              }                        
                        x <- list(as.numeric(x[,1]))
                        names(x) <- nm
                        x
                    }                      
                    else {
                        as.list(data.frame(x))
                    }
                }
            },data,names(data),SIMPLIFY=FALSE)
            names(data) <- NULL
            ee <- do.call("c",data)
        } else {
            stop("Argument data must be a matrix, data.frame or a list")
        }
    }
    ee
}

##' Gets the data which was used to estimate MIDAS regression
##'
##' A helper function. 
##' @title Get the data which was used to etimate MIDAS regression
##' @param object \code{midas_r} object
##' @return a named list with mixed frequency data
##' @author Vaidotas Zemlys
##' @export
get_estimation_sample <- function(object) {    
    nms <- all.vars(object$terms)
    dataenv <- eval(as.name("ee"),object$Zenv)
    if(is.null(dataenv))dataenv <- object$Zenv
    
    insample <- lapply(nms,function(nm)eval(as.name(nm),dataenv))
    names(insample) <- nms
    #The weights in the mls terms come up as variables, we do not need them
    insample[!sapply(insample,is.function)]
}

##' @importFrom forecast forecast
##' @name forecast
##' @rdname forecast.midas_r
##' @export
NULL

##' Plots MIDAS coefficients of a MIDAS regression for a selected term.
##'
##' Plots MIDAS coefficients of a selected MIDAS regression term together with corresponding MIDAS coefficients and their confidence intervals
##' of unrestricted MIDAS regression
##' @title Plot MIDAS coefficients
##' @param x \code{midas_r} object
##' @param term_name the term name for which the coefficients are plotted. Default is \code{NULL}, which selects the first MIDAS term
##' @param title the title string of the graph. The default is \code{NULL} for the default title.
##' @param vcov. the covariance matrix to calculate the standard deviation of the cofficients
##' @param unrestricted the unrestricted model, the default is unrestricted model from the \code{x} object. Set NULL to plot only the weights.
##' @param ... additional arguments passed to \code{vcov.}
##' @return a data frame with restricted MIDAS coefficients, unrestricted MIDAS coefficients and lower and upper confidence interval limits. The data
##' frame is returned invisibly.
##' @author Virmantas Kvedaras, Vaidotas Zemlys
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
##' plot_midas_coef(mr)
##' @export
plot_midas_coef <- function(x, term_name=NULL, title = NULL, vcov. = sandwich, unrestricted = x$unrestricted, ...) {
    if(is.null(term_name)) {
        wt <- do.call("rbind",lapply(x$term_info,function(l)c(l$term_name,l$weight_name)))
        wt <- data.frame(wt)
        colnames(wt) <- c("term_name","weight_name")
        wt <- wt[wt$weight_name != "", ]
        if(nrow(wt)==0) stop("No terms with MIDAS weights in midas_r object")
        if(nrow(wt)>1) warning("Multiple terms with MIDAS weights, choosing the first one. Please specify the desired term name via 'term_name' argument.")
        term_name <-wt$term_name[1]
    }
    ti <- x$term_info[[term_name]]
    mcoef <- coef(x, midas = TRUE)[ti$midas_coef_index]
    k <- length(mcoef)
    
    if(is.null(unrestricted)) {
        pd <- data.frame(restricted = mcoef, unrestricted=NA, lower=NA,upper=NA)
        plot(0:(k-1), pd$restricted, col = "blue", ylab = "MIDAS coefficients", xlab="High frequency lag", type="l")
        if(is.null(title)) {
            title(main = paste0("MIDAS coefficients for term ",term_name,": ",ti$weight_name))
        } else title(main = title)     
    } else {
        ucoef <- coef(unrestricted)[ti$midas_coef_index]        
        k <- length(mcoef)
        sdval <- sqrt(diag(vcov.(unrestricted, ...)))
        sdval <- sdval[ti$midas_coef_index]
        
        pd <- data.frame(restricted=mcoef, unrestricted=ucoef, lower=ucoef - 1.96*sdval,upper=ucoef + 1.96*sdval)
    
        ylim <- range(c(pd[,1],pd[,2],pd[,3],pd[,4]))
    
        plot(0:(k-1), pd$unrestricted, col="black", ylab="MIDAS coefficients", xlab="High frequency lag", ylim = ylim)
        if(is.null(title)) {
            title(main = paste0("MIDAS coefficients for term ",term_name,": ",ti$weight_name, " vs unrestricted"))
        } else title(main = title)
        points(c(0:(k - 1)), pd$restricted, type = "l", col = "blue")    
        points(c(0:(k - 1)), pd$lower, type = "l", col = "red", lty = 2)
        points(c(0:(k - 1)), pd$upper, type = "l", col = "red", lty = 2)
    }
    invisible(pd)
}


