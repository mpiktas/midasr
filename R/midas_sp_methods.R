##' @export
##' @method fitted midas_sp
fitted.midas_sp <- function(object, ...) {
    fit <- midas_sp_fit(object, Xeval = NULL, Zeval = NULL)
    as.numeric(fit$fitted.values)
}

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
    cat(" model:", deparse(x$call$formula),"\n")
    print(coef(x),digits = digits, ...)
    cat("\n")
    cat("Function", x$argmap_opt$Ofunction, "was used for fitting\n")
    invisible(x)
}


##' @export
##' @importFrom stats deviance pt pnorm residuals printCoefmat cor
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
    
    r_squared <- cor(f, object$model[, 1])
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
    ans <- list(formula = formula(object$call$formula), residuals=r, sigma=sqrt(resvar),
                df=c(p,rdf), cov.unscaled = V/resvar, call=object$call,
                coefficients=param,
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
    cat("\n Formula:", deparse(x$formula),"\n")
    df <- x$df
    rdf <- df[2L]
    cat("\n Parameters:\n")
    printCoefmat(x$coefficients,digits=digits,signif.stars=signif.stars,...)
    cat("\n Residual standard error:", format(signif(x$sigma,digits)), "on", rdf , "degrees of freedom\n")
    #   cat(" Multiple R-squared: ", formatC(x$r_squared, digits = digits))
    #   cat(",\tAdjusted R-squared: ", formatC(x$adj_r_squared, digits = digits),"\n")
    invisible(x)
}


##' @include midas_r_methods.R
setMethod("extract", signature = className("midas_sp", "midasr"), definition = extract.midas_r)

##' @method update midas_sp
##' @importFrom stats getCall setNames
##' @export
update.midas_sp <- function(object, formula.,..., evaluate = TRUE) {
    if (is.null(call <- getCall(object))) 
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) 
        call$formula <- update(Formula(object), formula.)
    
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call))) 
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }        
    }
    
    ##1. If no start is supplied update the start from the call
    ##2. If start is supplied intersect it with already fitted values.
    
    cf <- coef(object)
    ustart <- lapply(object$term_info,function(x)cf[x[["coef_index"]]])
    
    redo <- FALSE
    if(!("start" %in% names(extras))) {        
        if(!("start" %in% names(call) && is.null(call$start))) {
            call$start <- ustart
            object$start_opt <- cf
        }
    } else {
        cstart <- eval(call$start,object$Zenv)
        ustart[names(cstart)] <- cstart
        call$start <- ustart
        object$start_opt <- unlist(ustart)
    } 
    
    if (evaluate) {
        if(!missing(formula.) || "data" %in% names(extras)  || "weight_gradients" %in% names(extras) || redo) {
            eval(call, parent.frame())
        } else {
            ##If we got here, we assume that we do not need to reevaluate terms.
            if(!is.null(extras$Ofunction)) {
                Ofunction <- eval(extras$Ofunction)
                extras$Ofunction <- NULL
            } else Ofunction <- object$argmap_opt$Ofunction            
            dotargnm <- names(extras)
            if (length(dotargnm) > 0) {
                offending <- dotargnm[!dotargnm %in% names(formals(Ofunction))]
                if (length(offending) > 0) {
                    stop(paste("The function ", Ofunction, " does not have the following arguments: ", 
                               paste(offending, collapse = ", "), sep = ""))
                }
            }
            else {
                extras <- NULL
            }
            if (Ofunction != object$argmap_opt$Ofunction) {
                argmap <- c(list(Ofunction = Ofunction), extras)
            }
            else {
                argmap <- object$argmap_opt
                argmap$Ofunction <- NULL
                argnm <- union(names(argmap), names(extras))
                marg <- vector("list", length(argnm))
                names(marg) <- argnm
                marg[names(extras)] <- extras
                oldarg <- setdiff(names(argmap), names(extras))
                marg[oldarg] <- argmap[oldarg]
                argmap <- c(list(Ofunction = Ofunction), marg)
            }
            object$call <- call
            object$argmap_opt <- argmap
            midas_nlpr.fit(object)
        }
    }
    else call
}


##' Extracts various coefficients of MIDAS regression
##'
##' MIDAS regression has two sets of cofficients. The first set is the coefficients associated with the parameters
##' of weight functions associated with MIDAS regression terms. These are the coefficients of the NLS problem associated with MIDAS regression.
##' The second is the coefficients of the linear model, i.e  the values of weight
##' functions of terms, or so called MIDAS coefficients. By default the function returns the first set of the coefficients.
##' 
##' @title Extract coefficients of MIDAS regression
##' @param object \code{midas_nlpr} object
##' @param type one of plain, midas, or nlpr. Returns appropriate coefficients.
##' @param term_names a character vector with term names. Default is \code{NULL}, which means that coefficients of all the terms are returned
##' @param ... not used currently
##' @return a vector with coefficients
##' @author Vaidotas Zemlys
##' @method coef midas_sp
##' @rdname coef.midas_sp
##' @export
coef.midas_sp <- function(object, type = c("plain", "midas", "bw"), term_names = NULL, ...) {
    type <- match.arg(type)
    if (is.null(term_names) & type == "midas") stop("Please provide a term name to get midas type coefficients")
    if (is.null(term_names)) {
        if(type == "bw") {
            object$coefficients[1:length(object$bws)]
        }
        else return(object$coefficients)
    } else {
        if (type == "bw") {
            warning("Bandwith setting is the same for all terms") 
            coef(object, type = "bw", term_names = NULL)
        }
        if (length(setdiff(term_names,names(object$term_info))) > 0) stop("Some of the term names are not present in estimated MIDAS regression")
        if (type == "plain") {
            res <- lapply(object$term_info[term_names], function(x) {
                if (!is.null(x$coef_index)) object$coefficients[x$coef_index]
            })
        }  
        if (type == "midas") {
            res <- lapply(object$term_info[term_names], function(x) {
                if (!is.null(x$coef_index))  x$weight(object$coefficients[x$coef_index])
            })
        }
        names(res) <- NULL
        if (length(res) == 1) return(res[[1]])
        else return(unlist(res))
        
    }
}

##' @include midas_r_methods.R
##' @include midas_nlpr_methods.R
##' @method plot_midas_coef midas_sp
##' @export
plot_midas_coef.midas_sp <- plot_midas_coef.midas_nlpr

##' Plot non-parametric part of the single index MIDAS regression
##'
##' Plot non-parametric part of the single index MIDAS regression
##' of unrestricted MIDAS regression
##' @param x \code{midas_r} object
##' @param term_name the term name for which the coefficients are plotted. Default is \code{NULL}, which selects the first MIDAS term
##' @param title the title string of the graph. The default is \code{NULL} for the default title.
##' @param compare the parameters for weight function to compare with the model, default is NULL
##' @param ... not used
##' @return a data frame with restricted MIDAS coefficients, unrestricted MIDAS coefficients and lower and upper confidence interval limits. The data
##' frame is returned invisibly.
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @importFrom graphics plot points
##' @importFrom numDeriv jacobian
##' @importFrom stats na.omit approx density
##' @export
plot_sp <- function(x, term_name, title = NULL,  compare = NULL, ... ) {
    
    if (length(x$bws) > 1) stop("Plotting is currently supported for univariate non-parametric functions only")
    ti <- x$term_info[[term_name]]
    
    if (ti$term_type != "Z") stop("The term name supplied does not correspond to non-parametric term") 
    
    gfun <- function(p) {
        if(is.null(x$model_matrix_map$X)) { 
            xi <- 0 
        } else {
            xi <- x$model[, x$model_matrix_map$X] %*% x$coef_X(p)
        }
        Z <- x$model[, x$model_matrix_map$Z, drop = FALSE]
        if (is.null(ti$coef_index)) zi <- Z[, ti$midas_coef_index, drop = FALSE]
        else zi <- Z[, ti$midas_coef_index] %*% ti$weight(p[ti$coef_index])
        if (ncol(zi) == 1) zi <- as.numeric(zi)
        u <- x$model[, 1] - xi
        list(xi = xi, zi = zi, u = u, g = cv_np(u, zi, p[1:length(x$bws)], x$degree))
    }
    
    
    gg <- gfun(coef(x))
    
    kappa <- 0.5/pi # for standard Gaussian kernel
    dns <- density(gg$zi, na.rm = T) # choose a better bandwith here !!!
    f_z <- approx(dns$x, dns$y, xout = gg$zi)$y
    se_np <- sqrt(kappa*deviance(x)/f_z/(exp(coef(x)[1])*x$nobs))
    
    
    ozi <- order(gg$zi)
    
    se_npo <- se_np[ozi] 
    xio <- gg$zi[ozi]
    trmo <-  gg$g[ozi]
    pd <- data.frame(xi = xio,  term = trmo, compare = NA, lower = trmo - 1.96*se_npo, upper = trmo + 1.96*se_npo)
    
    if (!is.null(compare)) {
       pd$compare <- compare(pd$xi)
    }
    
    ylim <- range(na.omit(c(pd$term,pd$compare, pd$lower, pd$upper)))
    
    plot(pd$xi, pd$term, col = "blue", ylab = "Estimated non-parametric function", xlab = "MIDAS aggregate", ylim = ylim, type = "l")
    
    points(pd$xi, pd$compare, type = "l", col="black")
    points(pd$xi, pd$lower, type = "l", col = "grey", lty = 2)
    points(pd$xi, pd$upper, type = "l", col = "grey", lty = 2)
    
    if (is.null(title)) {
        title(main = paste0("Non-parametric estimate for term ",term_name," with weight: ",ti$weight_name))
    } else title(main = title) 
    
    invisible(pd)
    
}
