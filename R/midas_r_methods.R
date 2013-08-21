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
##' \code{predict.midas_r} produces predicted values, obtained by evaluating regression function in the frame \code{newdata}. This means that the apropriate model matrix is constructed using only the data in \code{newdata}. This makes this function not very convenient for forecasting purposes. If you want to supply the new data for forecasting horizon only use the function \link{forecast.midas_r}. Also this function produces only static predictions, if you want dynamic forecasts use the \link{forecast.midas_r}.
##' 
##' @title Predict method for MIDAS regression fit
##' @param object \code{\link{midas_r}} object
##' @param newdata a named list containing data for mixed frequencies. If ommitted, the fitted values are used.
##' @param na.action function determining what should be done with missing values in \code{newdata}. The most likely cause of missing values is the insufficient data for the lagged variables. The default is to ommit such missing values.
##' @param ... additional arguments, not used
##' @return a vector of predicted values
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @method predict midas_r
##' @rdname predict.midas_r
##' @export
predict.midas_r <- function(object, newdata, na.action = na.omit, ... ) {
    Zenv <- new.env(parent=parent.frame())
    
    if(missing(newdata))
      return(as.vector(fitted(object)))
    else {
        ee <- data_to_env(newdata)
        parent.env(ee) <- parent.frame()
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
    as.vector(X %*% midas_coef(object))
}

##' @export
##' @method summary midas_r
summary.midas_r <- function(object, ...) {
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
    tval <- param/se
    param <- cbind(param,se,tval,2*pt(abs(tval),rdf,lower.tail=FALSE))
    dimnames(param) <- list(pnames, c("Estimate", "Std. Error", 
        "t value", "Pr(>|t|)"))
    ans <- list(formula=formula(object$terms),residuals=r,sigma=sqrt(resvar),
                df=c(p,rdf), cov.unscaled=XDtXDinv, call=object$call,
                coefficients=param,midas.coefficients=midas_coef(object))
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
    invisible(x)
}

##' @export
##' @method print midas_r
print.midas_r <- function(x, digits=max(3,getOption("digits")-3),...) {
    cat("MIDAS regression model\n")
    cat(" model:", deparse(formula(x)),"\n")
    print(coef(x),digits = digits, ...)
    cat("\n")
    cat("Function", x$argmap.opt$Ofunction, "was used for fitting\n")
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

##' @import sandwich
##' @export
##' @method bread midas_r
bread.midas_r <- function(x,...) {
    sandwich:::bread.nls(x,...)
}

##' @export
##' @method vcov midas_r
vcov.midas_r <- function(object,...) {
    stats:::vcov.nls(object,...)
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

##' Return the coefficients of MIDAS regression
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the regression coefficients.
##' 
##' @param x an output from \code{\link{midas_r}}
##' @return vector with coefficients of MIDAS regression
##' @author Vaidotas Zemlys
##' @export
midas_coef <- function(x) {
    x$midas.coefficients
}

##' Return the estimated hyper parameters of the weight function(s)
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the estimated parameters of restriction function(s)
##' 
##' @param x an output from \code{\link{midas_r}}
##' @param name name of the restriction function, the default value is the names of the restriction functions supplied to \code{\link{midas_r}}
##' @return a list if \code{length(name)>1}, a vector otherwise
##' @author Vaidotas Zemlys
##' @export
weight_param <- function(x,name=weight_names(x)) {
    if(!any(name %in% names(x$param.map))) stop("Supply valid name(s) of the restriction function")
    res <- lapply(name,function(nm)coef(x)[x$param.map[[nm]]])
    names(res) <- name
    if(length(res)==1)res <- res[[1]]
    res
}
##' Return the restricted coefficients generated by restriction function(s)
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the restricted coefficients generated by restriction function(s)
##' 
##' @param x an output from \code{\link{midas_r}}
##' @param name name(s) of the restriction function(s), the default value is the name(s) of the restriction function(s) supplied to \code{\link{midas_r}}
##' @return a list if \code{length(name)>1}, a vector otherwise
##' @author Vaidotas Zemlys
##' @export
weight_coef <- function(x,name=weight_names(x)) {
    if(!any(name %in% names(x$param.map))) stop("Supply valid name(s) of the restriction function")
    res <- lapply(name,function(nm)x$weights[[nm]](weight_param(x,nm)))
    names(res) <- name
    if(length(res)==1)res <- res[[1]] 
    res
}
##' Return the names of restriction function(s)
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the name(s) of restriction function(s) used in call to \code{\link{midas_r}}
##' 
##' @param x an output from \code{\link{midas_r}}
##' @return a character vector
##' @author Vaidotas Zemlys
##' @export
weight_names <- function(x) {
    names(x$weights)
}
##' Return the coefficients for fmls variables
##'
##' Extracts the coefficients and returns those coefficients which name has string \code{fmls} or \code{mls} in it.
##' 
##' @param x an output from \code{\link{midas_u}}
##' @return a vector
##' @author Vaidotas Zemlys
##' @export
mls_coef <- function(x) {
    cf <- coef(x)
    cf[grep("fmls|mls|dmls",names(cf))]
}

##' Function to compute AICc information criteria for a given model
##'
##' A generic function for calculating AICc. It is implemented for \link{midas_r} package. For more examples see package AICcmodavg.
##' @title Compute AICc
##' @param mod \code{midas_r} model
##' @param ... additional parameters
##' @return a computed AICc value, a number.
##' @author Vaidotas Zemlys
##' @rdname AICc
AICc <- function(mod,...)UseMethod("AICc")

##' @export 
##' @method AICc midas_r
AICc.midas_r <- function(mod) {
    n <- length(fitted(mod))
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")
    -2 * LL + 2 * K * (n/(n - K - 1))
}

get_mls_info<- function(mt,Zenv) {
    vars <- as.list(attr(mt,"variables"))[-2:-1]
    res <- lapply(vars, function(l) {
        if(length(l)>1) {
            if(as.character(l[[1]])%in%c("mls","fmls","dmls")) {
                lags <- eval(l[[3]],Zenv)
                m <-eval(l[[4]],Zenv)              
                varname <- as.character(l[[2]])
                if(length(l)>=5)weight <- as.character(l[[5]])
                else weight <- NULL
                list(varname=varname,lags=lags,m=m,weight=weight)
            }
            else NULL
        }
        else NULL
    })    
    res[!sapply(res,is.null)]
}

##' Forecasts MIDAS regression given the future values of regressors. For dynamic models (with lagged response variable) there is an option to calculate dynamic forecast, when forecasted values of response variable are substituted into the lags of response variable.
##' 
##' Given future values of regressors this function combines the historical values used in the fitting the MIDAS regression model and calculates the forecasts.
##' 
##' @title Forecast MIDAS regression
##' @param object midas_r object
##' @param newdata newdata
##' @param method the forecasting method, either \code{"static"} or \code{"dynamic"}
##' @param insample a list containing the historic mixed frequency data 
##' @param ... additional arguments, not used
##' @return a vector of forecasts
##' @author Vaidotas Zemlys
##' @rdname forecast.midas_r
##' @export
forecast <- function(object,...) UseMethod("forecast") 

##' @rdname forecast.midas_r
##' @method forecast midas_r
##' @export
forecast.midas_r <- function(object,newdata=NULL,method=c("static","dynamic"),insample=get_estimation_sample(object),...) {

    method <- match.arg(method)
    ##Fix data_to_env to return the list if needed.
    outsample <- as.list(data_to_env(newdata))
    ##Get the name of response variable
    yname <- all.vars(object$terms[[2]])
   
    ##Get the frequency information of each variable    
    mlsinfo <- get_mls_info(object$terms,object$Zenv)
    freqinfo <- lapply(mlsinfo,with,m)
    names(freqinfo) <- sapply(mlsinfo,with,varname)

    freqinfo <- freqinfo[names(freqinfo)!=yname]
    lowfreqn <- setdiff(names(outsample),names(freqinfo))
    lowfreq <- as.list(rep(1,length(lowfreqn)))
    names(lowfreq) <- lowfreqn

    freqinfo <- c(freqinfo,lowfreq)
        
   
    if(!identical(sort(names(outsample)),sort(names(freqinfo)))) stop("Missing variables in newdata. Please supply the data for all the variables (excluding the response variable) in regression")
    freqinfo <- freqinfo[names(outsample)]
 
    h <- length(outsample[[names(freqinfo)[1]]])%/%freqinfo[[1]]
    
    
    if(method=="static") {
        if(!(yname %in% names(outsample))) {
            outsample <- c(list(rep(NA,h)),outsample)
            names(outsample)[1] <- yname
        }
        data <- try(rbind_list(insample,outsample))
        if(class(data)=="try-error")stop("Missing variables in newdata. Please supply the data for all the variables (excluding the response variable) in regression")
        res <- predict(object,newdata=data,na.action=na.pass)        
        n <- length(res)
        res[n+1-h:1]
    }
    else {
        outsample <- outsample[names(outsample)!=yname]
        yna <- list(NA)
        names(yna) <- yname
        res <- rep(NA,h)
        fdata <- insample
        for(i in 1:h) {
            ##Get the data for the current low frequency lag
            hout <- mapply(function(var,m){
                var[1:m+(i-1)*m]
            },outsample,freqinfo,SIMPLIFY=FALSE)
            hout <- c(yna,hout)
            fdata <- rbind_list(fdata,hout)
            if(class(fdata)=="try-error")stop("Missing variables in newdata. Please supply the data for all the variables (excluding the response variable) in regression")
            rr <- predict(object,newdata=fdata,na.action=na.pass)
            n <- length(rr)
            res[i] <- rr[n]
            fdata[[yname]][n] <- res[i]             
        }
        res
    }
}

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
    if(is.matrix(data)) data <- data.frame(data)
    if(is.data.frame(data)) {
        ee <- as.environment(as.list(data))
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
                    as.list(x)
                    }
            },data,names(data),SIMPLIFY=FALSE)
            names(data) <- NULL
            ee <- as.environment(do.call("c",data))
        } else {
            stop("Argument data must be a matrix, data.frame or a list")
        }
    }
    ee
}

##' Gets the data which was used to etimate MIDAS regression
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
