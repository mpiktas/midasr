##' @export
##' @method print midas_nlpr
print.midas_nlpr <- function(x, digits=max(3,getOption("digits")-3),...) {
    #Code adapted from dynln:::print.dynlm code
    model_string <-  "\nNon-linear parametric MIDAS regression model with \""
    cat(paste(model_string, class(x$lhs)[1], 
              "\" data:\n", sep = ""))
    cat(paste("Start = ", x$lhs_start, 
              ", End = ", x$lhs_end, 
              "\n", sep = ""))
    cat(" model:", deparse(formula(x)),"\n")
    print(coef(x),digits = digits, ...)
    cat("\n")
    cat("Function", x$argmap_opt$Ofunction, "was used for fitting\n")
    invisible(x)
}


##' @export
##' @importFrom stats deviance pt pnorm residuals printCoefmat
##' @method summary midas_nlpr
summary.midas_nlpr <- function(object, df=NULL, ...) {
    r <- as.vector(residuals(object))
    param <- coef(object)
    pnames <- names(param)
    n <- length(r)
    p <- length(param)
    rdf <- n - p
    
    resvar <- deviance(object)/rdf
    
    dg <- jacobian(object$rhs, param)
    V <- resvar * ginv(crossprod(dg)/nrow(dg))/n
    
    colnames(V) <- pnames
    rownames(V) <- pnames
   
    se <- sqrt(diag(V)) 
       
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
    ans <- list(formula=formula(object$terms), residuals=r, sigma=sqrt(resvar),
                df=c(p,rdf), cov.unscaled = V/resvar, call=object$call,
                coefficients=param,midas_coefficients=coef(object, midas = TRUE),
                r_squared = r_squared, adj_r_squared = adj_r_squared, lhs_start = object$lhs_start, lhs_end = object$lhs_end, class_lhs = class(object$lhs))
    class(ans) <- "summary.midas_nlpr"
    ans
}

##' Non-linear parametric MIDAS regression model deviance
##'
##' Returns the deviance of a fitted MIDAS regression object
##' @param object a \code{\link{midas_r}} object
##' @param ... currently nothing
##' @return The sum of squared residuals
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname deviance.midas_nlpr
##' @method deviance midas_nlpr
##' @export
deviance.midas_nlpr <- function(object,...) {
    sum(residuals(object)^2,na.rm=TRUE)
}

##' @export
##' @method print summary.midas_nlpr
print.summary.midas_nlpr <- function(x, digits=max(3, getOption("digits") - 3 ), signif.stars = getOption("show.signif.stars"), ...) {
    cat(paste("\nNon linear parametric MIDAS regression model with \"", x$class_lhs[1], 
              "\" data:\n", sep = ""))
    cat(paste("Start = ", x$lhs_start, 
              ", End = ", x$lhs_end, 
              "\n", sep = ""))
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

##' Predicted values based on \code{midas_nlpr} object.
##'
##' \code{predict.midas_nlpr} produces predicted values, obtained by evaluating regression function in the frame \code{newdata}. This means that the appropriate model matrix is constructed using only the data in \code{newdata}. This makes this function not very convenient for forecasting purposes. If you want to supply the new data for forecasting horizon only use the function \link{forecast.midas_r}. Also this function produces only static predictions, if you want dynamic forecasts use the \link{forecast.midas_r}.
##' 
##' @title Predict method for non-linear parametric MIDAS regression fit
##' @param object \code{\link{midas_nlpr}} object
##' @param newdata a named list containing data for mixed frequencies. If omitted, the in-sample values are used.
##' @param na.action function determining what should be done with missing values in \code{newdata}. The most likely cause of missing values is the insufficient data for the lagged variables. The default is to omit such missing values.
##' @param ... additional arguments, not used
##' @return a vector of predicted values
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @method predict midas_nlpr
##' @rdname predict.midas_nlpr
##' @export
predict.midas_nlpr <- function(object, newdata, na.action = na.omit, ... ) {
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
    
    rhs <- object$rhs
    
    environment(rhs)$X <- X
    xind1 <- environment(rhs)$xind1
    environment(rhs)$X1 <- X[, xind1]
    
    as.vector(rhs(coef(object)))
}
