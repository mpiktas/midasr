##' MIDAS regression model deviance
##'
##' Returns the deviance of a fitted MIDAS regression object
##' @param object a \code{\link{midas_r}} object
##' @param ... currently nothing
##' @return The sum of squared residuals
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname deviance
##' @method deviance midas_r
deviance.midas_r <- function(object,...) {
    sum(residuals(object)^2,na.rm=TRUE)
}

##' Predict method for MIDAS regression fit
##'
##' Predicted values based on MIDAS regression object
##' @param object \code{\link{midas_r}} object
##' @param newldata new low frequency data, must be \code{data.frame}
##' @param newhdata new high frequency data, must be \code{data.frame}
##' @param ... additional arguments, not used
##' @return a vector of predicted values
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @method predict midas_r
##' @export
predict.midas_r <- function(object, newldata, newhdata, ... ) {
    if(missing(newldata) | missing(newhdata))
      return(as.vector(fitted(object)))

    if(!is.data.frame(newldata)| !is.data.frame(newhdata)) stop("New data must supplied as data.frame")
    
    Zenv <- new.env(parent=parent.frame())

    ee <- as.environment(c(as.list(newldata),as.list(newhdata)))
    parent.env(ee) <- parent.frame()
    
    assign("ee",ee,Zenv)
    cll <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("object", "newldata"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf[[2L]] <- formula(object)
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- as.name("na.omit")
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


