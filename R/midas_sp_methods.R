##' @export
##' @method print midas_sp
print.midas_sp <- function(x, digits=max(3,getOption("digits")-3),...) {
    #Code adapted from dynln:::print.dynlm code
    model_string <-  "\nSemi-parametric MIDAS regression model with \""
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
##' @method summary midas_sp
summary.midas_sp <- function(object, df=NULL, ...) {
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
    mss <- sum((f - mean(f))^2)
    rss <- sum(r^2)
    
    n <- length(r)
    p <- length(coef(object))
    rdf <- n-p
    df.int <- 0L
    
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
    class(ans) <- "summary.midas_sp"
    ans
}

##' Semi-parametric MIDAS regression model deviance
##'
##' Returns the deviance of a fitted MIDAS regression object
##' @param object a \code{\link{midas_r}} object
##' @param ... currently nothing
##' @return The sum of squared residuals
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname deviance.midas_sp
##' @method deviance midas_sp
##' @export
deviance.midas_sp <- function(object,...) {
    sum(residuals(object)^2,na.rm=TRUE)
}

##' @export
##' @method print summary.midas_sp
print.summary.midas_sp <- function(x, digits=max(3, getOption("digits") - 3 ), signif.stars = getOption("show.signif.stars"), ...) {
    cat(paste("\nSemi-parametric MIDAS regression model with \"", x$class_lhs[1], 
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

##' @export
##' @method extract midas_sp
extract.midas_sp <- function(model, include.rsquared = TRUE, include.adjrs = TRUE, 
                               include.nobs = TRUE, include.rmse = TRUE, ...) {
    extract.midas_r(model, include.rsquared, include.adjrs, include.nobs, include.rmse, ...)
}
