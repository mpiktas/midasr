##' @importFrom stats fitted
##' @name fitted
NULL

##' @export
##' @method fitted midas_nlpr
fitted.midas_nlpr <- function(object, ...) {
    as.numeric(object$rhs(coef(object)))
}

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
    
    r_squared <- R2.np(object$model[,1], f)
    
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
                r_squared = r_squared, lhs_start = object$lhs_start, lhs_end = object$lhs_end, class_lhs = class(object$lhs))
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
##' @method coef midas_nlpr
##' @rdname coef.midas_nlpr
##' @export
coef.midas_nlpr <- function(object, type = c("plain", "midas", "nlpr"), term_names = NULL, ...) {
    type <- match.arg(type)
    if (is.null(term_names) & type != "plain") stop("Please provide a term name to get midas or nlpr type coefficients")
    if (is.null(term_names)) {
        return(object$coefficients)
    } else {
        if (length(setdiff(term_names,names(object$term_info))) > 0) stop("Some of the term names are not present in estimated MIDAS regression")
        if (type == "plain") {
            res <- lapply(object$term_info[term_names], function(x) {
                if (is.null(x[["nlpr"]])) object$coefficients[x$coef_index]
                else object$coefficients[x$coef_index][x$param_map$r]
            })
        }  
        if (type == "midas") {
            res <- lapply(object$term_info[term_names], function(x) {
                 if (is.null(x[["nlpr"]]))  x$weight(object$coefficients[x$coef_index])
                 else x$weight(object$coefficients[x$coef_index][x$param_map$r])
            })
        }
        if (type == "nlpr") {
            cf <- coef(object, type = "plain", term_name = NULL)
            res <- lapply(object$term_info[term_names], function(x) {
                if (is.null(x[["nlpr"]]))  stop("The term ", x$term_name, " is not a non-linear parametric term")
                else cf[x$coef_index][x$param_map$nlpr]
            })
        }
        
        
        names(res) <- NULL
        if (length(res) == 1) return(res[[1]])
        else return(unlist(res))
            
    }
}

##' @include midas_r_methods.R
setMethod("extract", signature = className("midas_nlpr", "midasr"), definition = extract.midas_r)


##' Plots MIDAS coefficients of a MIDAS regression for a selected term.
##'
##' Plots MIDAS coefficients of a selected MIDAS regression term together with corresponding MIDAS coefficients and their confidence intervals
##' of unrestricted MIDAS regression
##' @title Plot MIDAS coefficients
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
##' @importFrom stats na.omit
##' @method plot_midas_coef midas_nlpr
##' @rdname plot_midas_coef.midas_nlpr
##' @export
plot_midas_coef.midas_nlpr <- function(x, term_name = NULL, title = NULL,  compare = NULL, ...) {
    if(is.null(term_name)) {
        wt <- do.call("rbind",lapply(x$term_info,function(l)c(l$term_name,l$weight_name)))
        wt <- data.frame(wt, stringsAsFactors = FALSE)
        colnames(wt) <- c("term_name","weight_name")
        wt <- wt[wt$weight_name != "", ]
        if (nrow(wt) == 0) stop("No terms with MIDAS weights in midas_r object")
        if (nrow(wt) > 1) warning("Multiple terms with MIDAS weights, choosing the first one. Please specify the desired term name via 'term_name' argument.")
        term_name <- wt$term_name[1]
    }
    ti <- x$term_info[[term_name]]
   
    mcoef <- coef(x, type = "midas", term_name = term_name)
    pcoef <- coef(x, type = "plain", term_name = term_name)
    sx <- summary(x)
    V <- sx$cov.unscaled * sx$sigma^2
    ti <- x$term_info[[term_name]]
    if (is.null(ti$nlpr)) {
        Vind <- ti$coef_index
    } else {
        Vind <- ti$coef_index[ti$param_map$r]
    }
    grad_r <- jacobian(ti$weight, pcoef)
    V_s <- V[Vind, Vind]
    se_r <- sqrt(diag(grad_r %*% V_s %*% t(grad_r)))
    pd <- data.frame(restricted = mcoef, compare = NA, lower = mcoef - 1.96*se_r, upper = mcoef + 1.96*se_r, lag_struct = ti$lag_structure)
   
    if (!is.null(compare)) {
        pd$compare <- ti$weight(compare)
    }
    
    ylim <- range(na.omit(c(pd[,1],pd[,2], pd[,3],pd[,4])))
    plot(pd$lag_struct, pd$restricted, col = "blue", ylab = "MIDAS coefficients", xlab = "High frequency lag", ylim = ylim)
    
    points(pd$lag_struct, pd$compare, col = "black") 
    points(pd$lag_struct, pd$lower, type = "l", col = "grey", lty = 2)
    points(pd$lag_struct, pd$upper, type = "l", col = "grey", lty = 2)
    
    if (is.null(title)) {
        title(main = paste0("MIDAS coefficients for term ",term_name,": ",ti$weight_name))
    } else title(main = title) 
    
    invisible(pd)
}

##' Plots logistic function for LSTR MIDAS regression
##'
##' Plots logistic function for LSTR MIDSAS regression
##' of unrestricted MIDAS regression
##' @title Plot MIDAS coefficients
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
##' @importFrom stats na.omit
##' @export
plot_lstr <- function(x, term_name, title = NULL,  compare = NULL, ... ) {
    ti <- x$term_info[[term_name]]
    
    if (is.null(ti) && is.null(ti$nlpr)) stop("Please provide the name of the term which is nlpr term")
    
    sx <- summary(x)
    cf <- coef(x)
    r_map <- ti$coef_index[ti$param_map$r]
    n_map <- ti$coef_index[ti$param_map$nlpr[-2:-1]]
    Vind <- c(r_map,n_map)  
  
    tfun <- function(p) {
        ri <- 1:length(r_map)
        pr <- p[ri]
        pn <- p[-ri]
        as.vector(lstr_G(x$model[, ti$midas_coef_index + 1], ti$weight(pr), pn, ti$sd_x))
    }
    
    sx <- summary(x)
    V <- sx$cov.unscaled * sx$sigma^2
    
    
    grad_r <- jacobian(tfun, cf[Vind])
    
    V_s <- V[Vind, Vind]
    se_r <- sqrt(diag(grad_r %*% V_s %*% t(grad_r)))
    
    xi <- as.vector(x$model[, ti$midas_coef_index + 1] %*% ti$weight(cf[ti$coef_index][ti$param_map$r]))
    
    ixi <- order(xi)
    
    nt <- tfun(cf[Vind])
    
    xi_o <- xi[ixi]
    nt_o <- nt[ixi]
    se_ro <- se_r[ixi]
    
    pd <- data.frame(xi = xi_o, term = nt_o, compare = NA, lower = nt_o - 1.96*se_ro, upper = nt_o + 1.96*se_ro, xi2 = NA)
   
    if (!is.null(compare)) {
        op <- rep(NA, length(Vind))
        op[1:length(r_map)] <- compare$r
        op[-(1:length(r_map))] <- compare$lstr
        xi2 <- as.vector(x$model[, ti$midas_coef_index + 1] %*% ti$weight(op[1:length(r_map)]))
        ixi2 <- order(xi2)
        pd$compare <- tfun(op)[ixi2]
        
        pd$xi2 <-  xi2[ixi2]
    }
    
    ylim <- range(na.omit(c(pd$term,pd$compare, pd$lower, pd$upper)))
    
    plot(pd$xi, pd$term, col = "blue", ylab = "Logistic function", xlab = "MIDAS aggregate", ylim = ylim, type = "l")
    
    points(pd$xi2, pd$compare, col = "black", type = "l") 
    points(pd$xi, pd$lower, type = "l", col = "grey", lty = 2)
    points(pd$xi, pd$upper, type = "l", col = "grey", lty = 2)
    
    if (is.null(title)) {
        title(main = paste0("Logistic function of term ",term_name," with weight: ",ti$weight_name))
    } else title(main = title) 
    
    invisible(pd)
    
}

R2.np <- function(y, fit) {
    yc <- y - mean(y)
    fitc <- fit - mean(y)
    (sum(yc*fitc)^2)/(sum(yc^2)*sum(fitc^2))
}
