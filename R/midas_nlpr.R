##' Non-linear parametric MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares.
##'
##' @param formula formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param data a named list containing data with mixed frequencies
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen 
##' R function. Other elements of the list are the arguments passed to this function.  The default optimisation function is \code{\link{optim}} with arguments
##'  \code{method="Nelder-Mead"} and \code{control=list(maxit=5000)}. Other supported functions are \code{\link{nls}}, \code{\link{optimx}}.
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
##' \item{nobs}{the number of effective observations}
##' \item{convergence}{the convergence message}
##' \item{fitted.values}{the fitted values of MIDAS regression}
##' \item{residuals}{the residuals of MIDAS regression}
##' 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname midas_nlpr
##' @details Given MIDAS regression:
##'
##' \deqn{y_t = \sum_{j=1}^p\alpha_jy_{t-j} +\sum_{i=0}^{k}\sum_{j=0}^{l_i}\beta_{j}^{(i)}x_{tm_i-j}^{(i)} + u_t,}
##' 
##' estimate the parameters of the restriction
##'
##' \deqn{\beta_j^{(i)}=g^{(i)}(j,\lambda).}
##'
##' Such model is a generalisation of so called ADL-MIDAS regression. It is not required that all the coefficients should be restricted, i.e the function \eqn{g^{(i)}}
##' might be an identity function. Model with no restrictions is called U-MIDAS model. The regressors \eqn{x_\tau^{(i)}} must be of higher
##' (or of the same) frequency as the dependent variable \eqn{y_t}. 
##'
##' @importFrom stats as.formula formula model.matrix model.response terms lsfit time
##' @importFrom zoo index index2char
##' @export
midas_nlpr <- function(formula, data, start, Ofunction="optim", ...) {

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

    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
    args <- list(...)
    y <- as.numeric(model.response(mf, "numeric"))
    X <- model.matrix(mt, mf)
    
    #Save ts/zoo information
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
    
    prepmd <- prep_midas_nlpr(y,X,mt,Zenv,cl,args,start,Ofunction)
    
    prepmd <- c(prepmd, list(lhs = ysave, lhs_start = y_start, lhs_end = y_end))
    
    class(prepmd) <- "midas_nlpr"
    
    midas_nlpr.fit(prepmd)    
}

##' @method update midas_nlpr
##' @importFrom stats getCall update.formula setNames
##' @export
update.midas_nlpr <- function(object, formula.,..., evaluate = TRUE) {
    if (is.null(call <- getCall(object))) 
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) 
        call$formula <- update.formula(formula(object), formula.)
          
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


##' Fit restricted MIDAS regression
##'
##' Workhorse function for fitting restricted MIDAS regression
##'  
##' @param x \code{midas_r} object
##' @return \code{\link{midas_r}} object
##' @author Vaidotas Zemlys
midas_nlpr.fit <- function(x) {
    args <- x$argmap_opt
    function.opt <- args$Ofunction
    args$Ofunction <- NULL
   
    if(!(function.opt %in% c("optim","spg","optimx","nls","dry_run"))) 
        stop("The optimisation function is not in the supported functions list. Please see the midasr:::midas_nlpr.fit code for the supported function list")
    if(function.opt == "optim" | function.opt =="spg") {  
        args$par <- x$start_opt
        args$fn <- x$fn0
        opt <- try(do.call(function.opt,args),silent=TRUE)
        if(inherits(opt,"try-error")) {
            stop("The optimisation algorithm of MIDAS regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        par <- opt$par
        names(par) <- names(coef(x))
        x$convergence <- opt$convergence
    }
    if(function.opt=="optimx") {  
        args$par <- x$start_opt
        args$fn <- x$fn0
        opt <- try(do.call(function.opt,args),silent=TRUE)
        if(inherits(opt,"try-error")) {
            stop("The optimisation algorithm of MIDAS regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        bmet <- which.min(opt$value)
        par <- as.numeric(opt[bmet,1:length(args$par)])        
        names(par) <- names(coef(x))
        x$convergence <- opt$convcode[bmet]
    }    
    if(function.opt=="nls") {
        rhs <- x$rhs
        y <- x$model[,1]
        args$formula <- formula(y~rhs(p))
        args$start <- list(p=x$start_opt)
        opt <- try(do.call("nls",args),silent=TRUE)
        if(inherits(opt,"try-error")) {
            stop("The optimisation algorithm of MIDAS regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        par <- coef(opt)
        names(par) <- names(coef(x))
        x$convergence <- opt$convInfo$stopCode
    }
    if(function.opt == "dry_run") {
        opt <- NULL
        par <- x$start_opt
    }
    
    x$opt <- opt
    x$coefficients <- par
    if (inherits(x, "midas_sp")) {
        bws <- par[1:length(x$bws)]
        x$bws <- bws
    }
    names(par) <- NULL
    x$fitted.values <- x$rhs(par)
    x$residuals <- as.vector(x$model[,1]-x$fitted.values)
    x
}

## Prepare necessary objects for fitting of the non-linear parametric MIDAS regression
##
## y the response
## X the model matrix
## mt the terms of the formula
## Zenv the environment to evaluate the formula
## cl call of the function
## args additional argument
## start starting values
## Ofunction the optimisation function
## weight_gradients a list of gradient functions for weights
## lagsTable the lagstable from checkARstar
## unrestricted the unrestricted model
## guess_start if TRUE, get the initial values for non-MIDAS terms via OLS, if FALSE, initialize them with zero.
prep_midas_nlpr <- function(y, X, mt, Zenv, cl, args, start, Ofunction,  guess_start = TRUE) {

    
    start <- start[!sapply(start,is.null)]
    
    if(!is.null(args$guess_start)) {
        guess_start <- args$guess_start
        args$guess_start <- NULL
    }    
    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    rfd <- lapply(terms.lhs, dterm_nlpr, Zenv = Zenv)
    
    if (attr(mt,"intercept")==1)  {
        intc <- dterm_nlpr(expression(1), Zenv)
        intc$term_name <- "(Intercept)"
        rfd <- c(list(intc), rfd)
    }
    
    term_names <- sapply(rfd,"[[","term_name")
    if (length(setdiff(names(start), intersect(names(start), term_names))) > 0) {
        stop("The names for the starting values should coincide with terms in formula")
    }
    names(rfd) <- term_names
    
    rfd[names(start)] <- mapply(function(tmi, st) {
        tmi[["full_start"]] <- unlist(st)
        if (is.list(st)) {
            #This is for pretty names, remove all the previous names
            for (i in 1:length(st))names(st[[i]]) <- NULL
            if (!("r" %in% names(st))) stop("The starting values for the restriction should be in an element named r")
            if (setdiff(names(st),c("lstr","mmm")) != "r") stop("The starting values for nlpr term should be in an element named either lstr or mmm")
            tmi[["start"]] <- st[["r"]]
            nlpr_name <- setdiff(names(st),"r")
            tmi[["nlpr"]] <- eval(as.name(nlpr_name))
            bi <- build_indices(cumsum(sapply(st,length)), names(st))
            names(bi)[names(bi) == nlpr_name] <- "nlpr"
            tmi[["param_map"]] <- bi
            tmi
        } else {
            names(st) <- NULL
        }
        tmi[["full_start"]] <- unlist(st)
        tmi
    }, rfd[names(start)], start, SIMPLIFY = FALSE)    
    
    nlpr_terms <- names(which(sapply(rfd, function(l) !is.null(l[["nlpr"]]))))
    
    rf <- lapply(rfd,"[[","weight")
    
    fake_start_default <- lapply(rfd,"[[","start")

    fake_pinds <- build_indices_list(fake_start_default)

    fake_coef2 <- function(p) {              
            pp <- lapply(fake_pinds,function(x)p[x])     
            res <- mapply(function(fun,param)fun(param),rf,pp,SIMPLIFY=FALSE)
            return(res)
    }
    
    fake_midas_coef <- fake_coef2(unlist(fake_start_default)) 

    if(sum(is.na(unlist(fake_midas_coef)))>0) stop("Check your starting values, NA in midas coefficients") 
    
    npx <- cumsum(sapply(fake_midas_coef,length))
    xinds <- build_indices(npx,names(fake_start_default))

    start_default <- lapply(rfd, "[[", "full_start")
    pinds <- build_indices_list(start_default)
    
    pinds1 <- pinds[setdiff(names(pinds), nlpr_terms)]
    rf1 <- rf[setdiff(names(pinds), nlpr_terms)]
    
    rfd <- mapply(function(tmi, xind, pind){
        tmi[["xind"]] <- xind
        tmi[["pind"]] <- pind
        tmi[["sd_x"]] <- sd(X[,xind], na.rm = TRUE)
        tmi
    }, rfd, xinds, pinds, SIMPLIFY = FALSE)
    
    usual_terms <- 
    
    coef_list <- function(p, pi, rfl) {
        pp <- lapply(pi,function(x)p[x])     
        res <- mapply(function(fun,param)fun(param), rfl, pp, SIMPLIFY = FALSE)
        return(res)
    }
    coef1 <- function(p) unlist(coef_list(p, pi = pinds1, rfl = rf1))
   
    do_nlpr_term <- function(p, tfun, xind, wfun, pind, param_map, sd_x) {
        pr <- p[pind][param_map[["r"]]]
        pn <- p[pind][param_map[["nlpr"]]] 
        tfun(X[, xind], wfun(pr), pn, sd_x)
    }
    
    xind1 <- unlist(xinds[setdiff(names(xinds), nlpr_terms)])
    X1 <- X[, xind1, drop = FALSE]
    
    
    rhs <- function(p) { 
        T2 <- lapply(rfd[nlpr_terms],function(l) {
            do_nlpr_term(p, l[["nlpr"]], l[["xind"]], l[["weight"]], l[["pind"]], l[["param_map"]], l[["sd_x"]])
        })
        cf1 <- coef1(p)
        X1 %*% cf1 + Reduce("+", T2)
    }
  
    fn0 <- function(p,...) {
        r <- y - rhs(p)
        sum(r^2)
    }
    
    hess <- function(x)numDeriv::hessian(fn0,x)
      
    
    control <- c(list(Ofunction = Ofunction),args)
   
    ##The default method is "Nelder-Mead" and number of maximum iterations is 5000
    if (!("method" %in% names(control)) & Ofunction == "optim") {        
        control$method <- "Nelder-Mead"
        if (is.null(control$maxit)) control$maxit <- 5000
    } 
    #Do a rename to conform to midas_r
    term_info <- lapply(rfd, function(l) {
        nm <- names(l)
        nm[nm  == "xind"] <- "midas_coef_index"
        nm[nm == "pind"] <- "coef_index"
        names(l) <- nm
        l
    })
    
    
    list(coefficients = unlist(start_default),
         model = cbind(y,X),         
         fn0 = fn0,
         rhs = rhs,
         opt = NULL,
         argmap_opt = control,
         start_opt = unlist(start_default),
         start_list = start,
         call = cl,
         terms = mt,
         hessian = hess,
         Zenv = Zenv,
         term_info = term_info,
         nobs = nrow(X))   
}

##' LSTR (Logistic Smooth TRansition)  MIDAS regression
##'
##' Function for fitting LSTR MIDAS regression without the formula interface
##' @param y model response
##' @param X prepared matrix of high frequency variable lags for LSTR term
##' @param z additional low frequency variables
##' @param weight the weight function
##' @param start_lstr the starting values for lstr term
##' @param start_x the starting values for weight function
##' @param start_z the starting values for additional low frequency variables
##' @param method a method passed to \link{optimx}
##' @param ... additional parameters to \link{optimx}
##' @return an object similar to \code{midas_r} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @import numDeriv
##' @import optimx
##' @importFrom stats na.omit sd
##'
##' @export
##' 
midas_lstr_plain <- function(y, X, z = NULL, weight, start_lstr, start_x, start_z = NULL, method = c("Nelder-Mead", "BFGS"), ...) {
    d <- ncol(X)
   
    if(!is.null(z) && !is.matrix(z)) z <- matrix(z, ncol=1)
    
    model <- na.omit(cbind(y,X,z))
    
    y <- model[,1]
    X <- model[, 2:(ncol(X) + 1)]
    if (is.null(z)) { 
        z <- 0
    } else {
        z <- model[, (ncol(X) + 2):ncol(model)]
    }
    n <- nrow(model)
    
    sx <- length(start_x)
    sd_x <- sd(c(X))
    
    rhs <- function(p) {
        plstr <- p[1:4]
        pr <- p[5:(4 + sx)]
        if (is.null(z)) pz <- 0 else {
            pz <- p[(5 + sx):length(p)]
        }
        lstr(X, weight(pr, ncol(X)), plstr, sd_x) + z %*% pz  
    }
    
    fn0 <- function(p) {
        sum((y - rhs(p))^2)
    }
    
    start <- c(start_lstr, start_x, start_z)
    opt <- optimx(start, fn0, method = method,...)
    bmet <- which.min(opt$value)
    par <- as.numeric(opt[bmet, 1:length(start)])   
    call <- match.call()
    fitted.values <- as.vector(y - rhs(par))
    names(par) <- c(paste0("lstr",1:length(start_lstr)), paste0("x", 1:length(start_x)), paste0("z", 1:length(start_z)))
    list(coefficients = par,
         midas_coefficients = weight(par[5:(4 + sx)],ncol(X)),
         lstr_coefficients = par[1:4],
         model = model,
         weights = weight,
         fn0 = fn0,
         rhs = rhs,
         opt = opt,
         call = call,
         hessian = function(x)numDeriv::hessian(fn0,x),
         fitted.values = fitted.values,
         residuals = as.vector(y - fitted.values),
         start = start)
             
}

##' MMM (Mean-Min-Max)  MIDAS regression
##'
##' Function for fitting MMM MIDAS regression without the formula interface
##' @param y model response
##' @param X prepared matrix of high frequency variable lags for MMM term
##' @param z additional low frequency variables
##' @param weight the weight function
##' @param start_mmm the starting values for MMM term
##' @param start_x the starting values for weight function
##' @param start_z the starting values for additional low frequency variables
##' @param method a method passed to \link{optimx}
##' @param ... additional parameters to \link{optimx}
##' @return an object similar to \code{midas_r} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @import numDeriv
##' @import optimx
##' @importFrom stats na.omit sd
##'
##' @export
##' 
midas_mmm_plain <- function(y, X, z = NULL, weight, start_mmm, start_x, start_z = NULL, method = c("Nelder-Mead", "BFGS"), ...) {
    d <- ncol(X)
    
    if(!is.null(z) && !is.matrix(z)) z <- matrix(z, ncol=1)
    
    model <- na.omit(cbind(y,X,z))
    
    y <- model[,1]
    X <- model[, 2:(ncol(X) + 1)]
    if (is.null(z)) { 
        z <- 0
    } else {
        z <- model[, (ncol(X) + 2):ncol(model)]
    }
    n <- nrow(model)
    
    sx <- length(start_x)
    sd_x <- sd(c(X))
    
    rhs <- function(p) {
        pmmm <- p[1:2]
        pr <- p[3:(2 + sx)]
        if (is.null(z)) pz <- 0 else {
            pz <- p[(3 + sx):length(p)]
        }
        mmm(X, weight(pr, ncol(X)), pmmm) + z %*% pz  
    }
    
    fn0 <- function(p) {
        sum((y - rhs(p))^2)
    }
    
    start <- c(start_mmm, start_x, start_z)
    opt <- optimx(start, fn0, method = method,...)
    bmet <- which.min(opt$value)
    par <- as.numeric(opt[bmet, 1:length(start)])   
    call <- match.call()
    fitted.values <- as.vector(y - rhs(par))
    names(par) <- c(paste0("mm",1:length(start_mmm)), paste0("x", 1:length(start_x)), paste0("z", 1:length(start_z)))
    list(coefficients = par,
         midas_coefficients = weight(par[3:(2 + sx)],ncol(X)),
         mmm_coefficients = par[1:2],
         model = model,
         weights = weight,
         fn0 = fn0,
         rhs = rhs,
         opt = opt,
         call = call,
         hessian = function(x)numDeriv::hessian(fn0,x),
         fitted.values = fitted.values,
         residuals = as.vector(y - fitted.values),
         start = start)
    
}

#' Compute LSTR term for high frequency variable
#'
#' @param X matrix, high frequency variable embedded in low frequency, output of mls
#' @param theta vector, restriction coefficients for high frequency variable
#' @param beta vector of length 4, parameters for LSTR term, slope and 3 LSTR parameters
#' @param sd_x vector of length 1, defaults to standard deviation of X.
#'
#' @return a vector
#' @export
#'
lstr <- function(X, theta, beta, sd_x = sd(c(X), na.rm = TRUE)) {
    xx <- X %*% theta
    b <- -exp(beta[3])*(xx - beta[4])/sd_x
    G <- 1/(1 + exp(b))
    beta[1]*xx*(1 + beta[2]*G) 
}

lstr_G <- function(X, theta, beta, sd_x = sd(c(X), na.rm = TRUE)) {
    xx <- X %*% theta
    b <- -exp(beta[1])*(xx - beta[2])/sd_x
    1/(1 + exp(b))
}
#' Compute MMM term for high frequency variable
#'
#' @param X matrix, high frequency variable embedded in low frequency, output of mls
#' @param theta vector, restriction coefficients for high frequency variable
#' @param beta vector of length 2, parameters for MMM term, slope and MMM parameter.
#' @param ..., currently not used
#'
#' @return a vector
#' @export
#'
mmm <- function(X, theta, beta, ...) {
    mtr <- exp(beta[2]*X)
    mtr_denom <- apply(mtr, 1, sum)
    mmm_term <- ncol(X)*X*mtr/mtr_denom
    
    beta[1]*mmm_term %*% theta 
}

dterm_nlpr <- function(fr, Zenv) {
    term_name <- as.character(fr)[1]
    weight_name <- ""
    rf <- function(p)p
    start <- 0
    freq <- 1
    lagstruct <- 0
    if(term_name %in% c("mls", "fmls", "dmls","mlsd")) {
        type <- term_name
        term_name <- as.character(fr[[2]])
        
        wpos <- 5
        if(type == "mlsd") {
            freq <- NA
        } else {
            freq <- eval(fr[[4]], Zenv)
        }
        
        lags <- eval(fr[[3]], Zenv)
        nol <- switch(type,
                      fmls = lags+1,
                      dmls = lags+1,
                      mls = length(lags),
                      mlsd = length(lags)
        )
        lagstruct <- switch(type,
                            fmls = 0:lags,
                            dmls = 0:lags,
                            mls = lags,
                            mlsd = lags
        )
        start <- rep(0, nol)
        if(length(fr) > wpos - 1) {
            mf <- fr[-wpos]
            mf[[1]] <- fr[[wpos]]
            weight_name <- as.character(fr[[wpos]])
            
            noarg <- length(formals(eval(fr[[wpos]], Zenv)))
            if(noarg < 2) stop("The weight function must have at least two arguments")            
            mf <- mf[1:min(length(mf), noarg + 1)]
            if(length(mf)>3) {
                start_eval <- 4
                if(type == "mlsd") start_eval <- 5
                if(length(mf)>=start_eval) {
                    for(j in start_eval:length(mf)) {
                        mf[[j]] <- eval(mf[[j]], Zenv)
                    }
                }
            }
            mf[[3]] <- nol
            rf <- function(p) {
                mf[[2]] <- p
                eval(mf,Zenv)
            }
        }
    }
    list(weight = rf,
         term_name = term_name,
         start = start,
         full_start = start, 
         weight_name = weight_name,
         frequency = freq,
         lag_structure = lagstruct
    )
}

