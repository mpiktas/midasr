##' @export
##' @method print midas_r
print.midas_qr <- function(x, digits=max(3,getOption("digits") - 3), ...) {
    print.midas_r(x)
}

##' @export
##' @method coef midas_r
coef.midas_qr <- function(object, midas = FALSE, term_names = NULL, ...) {
    if (is.null(term_names)) coef.midas_r(object, midas, ...)
    else {
        if (length(object$tau) > 1) {
            mc <- object$midas_coefficients
            if (length(setdiff(term_names,names(object$term_info))) > 0) 
                stop("Some of the term names are not present in estimated MIDAS regression")
            if (midas) {
                res <- lapply(object$term_info[term_names], function(x) object$midas_coefficients[x$midas_coef_index,])
            } else {
                res <- lapply(object$term_info[term_names], function(x) object$coefficients[x$coef_index,])
            }
            names(res) <- NULL
            if (length(res) == 1) return(res[[1]])
            else return(unlist(res))
        } else coef.midas_r(object, midas, term_names, ...)
    }        
 
}

## Code is borrowed from summary.nlrq in quantreg
##' @export
summary.midas_qr <- function(object, ... ) {
    y <- object$residuals
    tau <- object$tau
    cfs <- coef(object)
    
    do_f <- function(yy, cf, tau) {
        XX <- object$model[, -1] %*% object$gradD(cf)
        f <- summary(rq(yy ~ XX - 1, tau), se = "boot", covariance = TRUE, ...)
        f$coefficients[, 1] <- cf
        f$coefficients[, 3] <- f$coefficients[, 1]/f$coefficients[, 2]
        f$coefficients[, 4] <- if (f$rdf > 0) 2 * (1 - pt(abs(f$coef[, 3]), f$rdf))
        dimnames(f$coefficients)[[1]] <- names(cf)
        f$call <- object$call
        f$tau <- tau
        f
    }
    
    if (length(tau) > 1) {
        all_f <- vector(mode = "list", length(tau))
        for (i in 1:length(tau)) {
            cff <- as.vector(cfs[,i])
            names(cff) <- rownames(cfs)
            all_f[[i]] <- do_f(as.vector(y[, i]), cff, tau[i])
        }
        cf <- lapply(all_f, coef)
    } else {
        all_f <- do_f(as.vector(y), cfs, tau)
        cf <- coef(all_f)
    }
    
    ans <- list(formula = formula(object$terms), residuals = y,
                call = object$call,
                f = all_f,
                tau = object$tau,
                coefficients = cf, midas_coefficients = coef(object, midas = TRUE),
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
    cat("\n Parameters:\n")
    cf <- coef(x)
    if (length(x$tau) > 1) {
        for (i in 1:length(x$tau)) {
            cat("\n tau:", x$tau[i], "\n")    
            printCoefmat(cf[[i]], digits = digits, signif.stars = signif.stars,...)
        }
    } else {
        cat("\n tau:", x$tau, "\n")    
        printCoefmat(coef(x), digits = digits, signif.stars = signif.stars,...)
    }
    invisible(x)   
}