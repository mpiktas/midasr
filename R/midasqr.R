##' Restricted MIDAS quantile regression
##'
##' Estimate restricted MIDAS quantile regression using nonlinear quantile regression
##'
##' @param formula formula for restricted MIDAS regression or \code{midas_qr} object. Formula must include \code{\link{mls}} function
##' @param data a named list containing data with mixed frequencies
##' @param tau quantile
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are the arguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param weight_gradients a named list containing gradient functions of weights. The weight gradient function must return the matrix with dimensions
##' \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the number of coefficients in unrestricted and restricted regressions correspondingly.
##' The names of the list should coincide with the names of weights used in formula.
##' The default value is NULL, which means that the numeric approximation of weight function gradient is calculated. If the argument is not NULL, but the
##' name of the weight used in formula is not present, it is assumed that there exists an R function which has  
##' the name of the weight function appended with \code{_gradient}. 
##' @param guess_start, logical, if TRUE tries certain strategy to improve starting values 
##' @param ... additional arguments supplied to optimisation function
##' @return a \code{midas_r} object which is the list with the following elements:
##' 
##' \item{coefficients}{the estimates of parameters of restrictions}
##' \item{midas_coefficients}{the estimates of MIDAS coefficients of MIDAS regression}
##' \item{model}{model data}
##' \item{unrestricted}{unrestricted regression estimated using \code{\link{midas_u}}}
##' \item{term_info}{the named list. Each element is a list with the information about the term, such as its frequency, function for weights, gradient function of weights, etc.}
##' \item{fn0}{optimisation function for non-linear least squares problem solved in restricted MIDAS regression}
##' \item{rhs}{the function which evaluates the right-hand side of the MIDAS regression}
##' \item{gen_midas_coef}{the function which generates the MIDAS coefficients of MIDAS regression}
##' \item{opt}{the output of optimisation procedure}
##' \item{argmap_opt}{the list containing the name of optimisation function together with arguments for optimisation function}
##' \item{start_opt}{the starting values used in optimisation}
##' \item{start_list}{the starting values as a list}
##' \item{call}{the call to the function}
##' \item{terms}{terms object}
##' \item{gradient}{gradient of NLS objective function}
##' \item{hessian}{hessian of NLS objective function}
##' \item{gradD}{gradient function of MIDAS weight functions} 
##' \item{Zenv}{the environment in which data is placed}
##' \item{use_gradient}{TRUE if user supplied gradient is used, FALSE otherwise}
##' \item{nobs}{the number of effective observations}
##' \item{convergence}{the convergence message}
##' \item{fitted.values}{the fitted values of MIDAS regression}
##' \item{residuals}{the residuals of MIDAS regression}
##' 
##' @examples 
##' ##Take the same example as in midas_r
##' 
##' theta_h0 <- function(p, dk, ...) {
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
##' ##Generate the predictor variable
##' xx <- ts(arima.sim(model = list(ar = 0.6), 600 * 12), frequency = 12)
##'
##' ##Simulate the response variable
##' y <- midas_sim(500, xx, theta0)
##'
##' x <- window(xx, start=start(y))
##' 
##' ##Fit quantile regression. All the coefficients except intercept should be constant.
##' ##Intercept coefficient should correspond to quantile function of regression errors.
##' mr <- midas_qr(y~fmls(x,4*12-1,12,theta_h0), tau = c(0.1, 0.5, 0.9),
##'               list(y=y,x=x),
##'               start=list(x=c(-0.1,10,-10,-10)))
##'               
##' mr
##' @author Vaidotas Zemlys-Balevicius
##' @rdname midas_qr
##' @import quantreg
##' @export
midas_qr <- function(formula, data, tau = 0.5, start, Ofunction="nlrq", weight_gradients=NULL, guess_start = TRUE, ...) {
    Zenv <- new.env(parent=environment(formula))
    
    if(missing(data)) {
        ee <- NULL
    }
    else {
        ee <- data_to_env(data)
        parent.env(ee) <- parent.frame()
    }
    
    if(missing(start)) {
        stop("Please supply starting values.")
    } 
    
    assign("ee",ee,Zenv)
    formula <- as.formula(formula)
    cl <- match.call()    
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- formula
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- as.name("na.omit")
    names(mf)[c(2,3,4)] <- c("formula","data","na.action")
    
    ##Add check for dynamic terms. 
    ##They are not supported as they do not make sense in quantile regression
    itr <- checkARstar(terms(eval(mf[[2]], Zenv)))
    if(!is.null(itr$lagsTable)) mf[[2]] <- itr$x
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
    args <- list(...)
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)
   
    if(is.null(ee)) { 
        yy <- eval(formula[[2]], Zenv)
    }else {
        yy <- eval(formula[[2]], ee)
    }
    
    y_index <- 1:length(yy) 
    if(!is.null(attr(mf, "na.action"))) {
        y_index <- y_index[-attr(mf, "na.action")]
    }
    if(length(y_index)>1) {
        if(sum(abs(diff(y_index) - 1))>0) warning("There are NAs in the middle of the time series")                
    }
    
    ysave <- yy[y_index]
    
    if(inherits(yy, "ts")) {
        class(ysave) <- class(yy)
        attr(ysave, "tsp") <- c(time(yy)[range(y_index)], frequency(yy))
    }
    
    if(inherits(yy,c("zoo","ts"))) {
        y_start <- index2char(index(ysave)[1], frequency(ysave))
        y_end <- index2char(index(ysave)[length(ysave)], frequency(ysave))
    } else {
        y_start <- y_index[1]
        y_end <- y_index[length(y_index)]
    }
    
    eps <- .Machine$double.eps^(2/3)
    if (any(tau < 0) || any(tau > 1)) 
        stop("invalid tau:  taus should be >= 0 and <= 1")
    
    if (any(tau == 0)) 
            tau[tau == 0] <- eps
    if (any(tau == 1)) 
            tau[tau == 1] <- 1 - eps
    
    prepmd <- prepmidas_r(y,X,mt,Zenv,cl,args,start,Ofunction,weight_gradients,itr$lagsTable, guess_start = guess_start, tau = tau)
    
    prepmd <- c(prepmd, list(lhs = ysave, lhs_start = y_start, lhs_end = y_end ))
    
    class(prepmd) <- c("midas_qr")
    
    midas_qr.fit(prepmd)    
    
}

midas_qr.fit <- function(x) {
    args <- x$argmap_opt
    function.opt <- args$Ofunction
    args$Ofunction <- NULL
    if(function.opt=="nlrq") {
        rhs <- x$rhs
        if(x$use_gradient) {
            orhs <- rhs
            rhs <- function(p) {
                res <- orhs(p)
                attr(res,"gradient") <- x$model[,-1]%*%x$gradD(p)
                res
            }
        }
        z <- x$model[,1]
        args$formula <- formula(z~rhs(p))
        args$start <- list(p=x$start_opt)
        res <- vector(mode = "list", length = length(x$tau))
        for(i in 1:length(x$tau)) {
            args$tau <- x$tau[i]
            opt <- try(do.call("nlrq",args),silent = TRUE)
                if(inherits(opt,"try-error")) {
                    stop("The optimisation algorithm of MIDAS regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
                }
            par <- coef(opt)
            names(par) <- names(coef(x))
            res[[i]] <- list(opt = opt, par = par, convergence = opt$convInfo$stopCode)
        }
    }
    if (function.opt == "dry_run") {
        opt <- NULL
        res_template <- list(opt = opt, par = x$start_opt, convergence = "Dry run, no optimisation done")
        res <- vector(mode = "list", length = length(x$tau))
        for(i in 1:length(x$tau)) {
            res[[i]] <- res_template
        }
        
    }
    if (length(res) == 1) {
        x$opt <- res[[1]]$opt
        x$coefficients <- res[[1]]$par
        names(res[[1]]$par) <- NULL
        x$midas_coefficients <- x$gen_midas_coef(res[[1]]$par)
        x$fitted.values <- as.vector(x$model[,-1]%*%x$midas_coefficients)
        x$residuals <- as.vector(x$model[,1]-x$fitted.values)
    } else {
        x$opt <- lapply(res, "[[", "opt")
        x$convergence <- sapply(res, "[[", "convergence")
        x$coefficients <- sapply(res, "[[", "par")
        colnames(x$coefficients) <- x$tau
        par <- sapply(res, "[[", "par")
        rownames(par) <- NULL
        x$midas_coefficients <- apply(par,2,x$gen_midas_coef)
        colnames(x$midas_coefficients) <- x$tau
        x$fitted.values <- x$model[,-1]%*%x$midas_coefficients
        x$residuals <- x$model[,1]-x$fitted.values
    }
    x
    
}
