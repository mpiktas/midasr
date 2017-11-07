##' @export
##' @method print midas_r
print.midas_qr <- function(x, digits=max(3,getOption("digits")-3), ...) {
    print.midas_r(x)
}

##' @export
##' @method coef midas_r
coef.midas_qr <- function(object, midas = FALSE, term_names = NULL, ...) {
    if(is.null(term_names)) coef.midas_r(object, midas, ...)
    else {
        if(length(object$tau) > 1) {
            mc <- object$midas_coefficients
            if(length(setdiff(term_names,names(object$term_info)))>0) stop("Some of the term names are not present in estimated MIDAS regression")
            if(midas) {
                res <- lapply(object$term_info[term_names], function(x) object$midas_coefficients[x$midas_coef_index,])
            } else {
                res <- lapply(object$term_info[term_names], function(x) object$coefficients[x$coef_index,])
            }
            names(res) <- NULL
            if(length(res) ==1) return(res[[1]])
            else return(unlist(res))
        } else coef.midas_r(object, midas, term_names, ...)
    }        
 
}

## Code is borrowed from summary.nlrq in quantreg
##' @export
summary.midas_qr <- function(object, ... ) {
    y <- as.vector(object$residuals)
    X <- object$model[, -1] %*% object$gradD(coef(object))
    tau <- object$tau
    
    f <- summary(rq(y ~ X - 1, tau), se = "boot", covariance = TRUE, ...)
    f$coefficients[, 1] <- coef(object)
    f$coefficients[, 3] <- f$coefficients[, 1]/f$coefficients[, 2]
    f$coefficients[, 4] <- if (f$rdf > 0) 
        2 * (1 - pt(abs(f$coef[, 3]), f$rdf))
    dimnames(f$coefficients)[[1]] <- names(coef(object))
    f$call <- object$call
    f$tau <- tau
    
    ans <- list(formula=formula(object$terms), residuals=y,
                call=object$call,
                f = f,
                tau = object$tau,
                coefficients=coef(f), midas_coefficients=coef(object, midas = TRUE),
                lhs_start = object$lhs_start, lhs_end = object$lhs_end, class_lhs = class(object$lhs))
    
    class(ans) <- "summary.midas_qr"
    ans
}

##' @export
print.summary.midas_qr <- function(x, digits =max(3, getOption("digits") - 3 ), 
                                   signif.stars = getOption("show.signif.stars"), ...)  {
    
    cat(paste("\nMIDAS quantile regression model with \"", x$class_lhs[1], 
              "\" data:\n", sep = ""))
    cat(paste("Start = ", x$lhs_start, 
              ", End = ", x$lhs_end, 
              "\n", sep = ""))
    cat("\n Formula", deparse(formula(x)),"\n")
    
    cat("\n tau:", x$tau, "\n")
    cat("\n Parameters:\n")
    printCoefmat(coef(x), digits = digits, signif.stars = signif.stars,...)
    invisible(x)   
}