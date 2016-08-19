##' Estimate unrestricted MIDAS regression
##'
##' Estimate unrestricted MIDAS regression using OLS. This function is a wrapper for \code{lm}.
##' 
##' @param formula MIDAS regression model formula
##' @param data  a named list containing data with mixed frequencies
##' @param ... further arguments, which could be passed to \code{\link{lm}} function.
##' @return \code{\link{lm}} object.
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Kvedaras V., Zemlys, V. \emph{Testing the functional constraints on parameters in regressions with variables of different frequency} Economics Letters 116 (2012) 250-254 
##' @examples
##' ##The parameter function
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
##' ##Do not run
##' #plot(theta0)
##'
##' ##' ##Generate the predictor variable
##' xx <- ts(arima.sim(model = list(ar = 0.6), 600 * 12), frequency = 12)
##'
##' ##Simulate the response variable
##' y <- midas_sim(500, xx, theta0)
##'
##' x <- window(xx, start=start(y))
##' 
##' ##Create low frequency data.frame
##' ldt <- data.frame(y=y,trend=1:length(y))
##'
##' ##Create high frequency data.frame
##'
##' hdt <- data.frame(x=window(x, start=start(y)))
##' 
##' ##Fit unrestricted model
##' mu <- midas_u(y~fmls(x,2,12)-1, list(ldt, hdt))
##'
##' ##Include intercept and trend in regression
##'
##' mu_it <- midas_u(y~fmls(x,2,12)+trend, list(ldt, hdt))
##'
##' ##Pass data as partialy named list
##'
##' mu_it <- midas_u(y~fmls(x,2,12)+trend, list(ldt, x=hdt$x))
##' 
##' @details MIDAS regression has the following form:
##' 
##' \deqn{y_t = \sum_{j=1}^p\alpha_jy_{t-j} +\sum_{i=0}^{k}\sum_{j=0}^{l_i}\beta_{j}^{(i)}x_{tm_i-j}^{(i)} + u_t,}
##'
##' where \eqn{x_\tau^{(i)}}, \eqn{i=0,...k} are regressors of higher (or similar) frequency than \eqn{y_t}. 
##' Given certain assumptions the coefficients can be estimated using usual OLS and they have the familiar properties associated with simple linear regression.
##'
##' @export
midas_u <- function(formula, data ,...) {
    Zenv <- new.env(parent=environment(formula))

    if(missing(data)) {
        ee <- NULL
    }
    else {
        ee <- data_to_env(data)
        parent.env(ee) <- parent.frame()
    }
    
    assign("ee",ee,Zenv)    
    mf <- match.call(expand.dots = TRUE)
    mf <- mf[-4]
    mf[[1L]] <- as.name("lm")
    mf[[3L]] <- as.name("ee")   
   
    if(is.null(ee)) { 
        yy <- eval(formula[[2]], Zenv)
    }else {
        yy <- eval(formula[[2]], ee)
    }
    
    if(inherits(yy, "ts")) {
        y_index <- 1:length(yy) 
        if(!is.null(attr(mf, "na.action"))) {
            y_index <- y_index[-attr(mf, "na.action")]
        }
        if(length(y_index)>1) {
            if(sum(abs(diff(y_index) - 1))>0) warning("There are NAs in the middle of the time series")                
        }
        ysave <- yy[y_index]
        class(ysave) <- class(yy)
        attr(ysave, "tsp") <- c(time(yy)[range(y_index)], frequency(yy))
    } else {
        ysave <- yy
    }
    
    out <- eval(mf,Zenv)
    out$Zenv <- Zenv
    out$midas_coefficients <- out$coefficients
    out$lhs <- ysave
    class(out) <- c("midas_u",class(out))
    
    out
}

##' Restricted MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares.
##'
##' @param formula formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param data a named list containing data with mixed frequencies
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are the arguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param weight_gradients a named list containing gradient functions of weights. The weight gradient function must return the matrix with dimensions
##' \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the number of coefficients in unrestricted and restricted regressions correspondingly.
##' The names of the list should coincide with the names of weights used in formula.
##' The default value is NULL, which means that the numeric approximation of weight function gradient is calculated. If the argument is not NULL, but the
##' name of the weight used in formula is not present, it is assumed that there exists an R function which has  
##' the name of the weight function appended with \code{_gradient}. 
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
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Clements, M. and Galvao, A., \emph{Macroeconomic Forecasting With Mixed-Frequency Data: Forecasting Output Growth in the United States}, Journal of Business and Economic Statistics, Vol.26 (No.4), (2008) 546-554
##' @rdname midas_r
##' @examples
##' ##The parameter function
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
##' ##Fit restricted model
##' mr <- midas_r(y~fmls(x,4*12-1,12,theta_h0)-1,
##'               list(y=y,x=x),
##'               start=list(x=c(-0.1,10,-10,-10)))
##'
##' ##Include intercept and trend in regression
##' mr_it <- midas_r(y~fmls(x,4*12-1,12,theta_h0)+trend,
##'                  list(data.frame(y=y,trend=1:500),x=x),
##'                  start=list(x=c(-0.1,10,-10,-10)))
##' 
##' data("USrealgdp")
##' data("USunempr")
##' 
##' y.ar <- diff(log(USrealgdp))
##' xx <- window(diff(USunempr), start = 1949)
##' trend <- 1:length(y.ar)
##' 
##' ##Fit AR(1) model
##' mr_ar <- midas_r(y.ar ~ trend + mls(y.ar, 1, 1) +
##'                  fmls(xx, 11, 12, nealmon),
##'                  start = list(xx = rep(0, 3)))
##' 
##' ##First order MIDAS-AR* restricted model 
##' mr_arstar <-  midas_r(y.ar ~ trend + mls(y.ar, 1, 1, "*")
##'                      + fmls(xx, 11, 12, nealmon),
##'                      start = list(xx = rep(0, 3)))
##'
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
##' MIDAS-AR* (a model with a common factor, see (Clements and Galvao, 2008)) can be estimated by specifying additional argument, see an example.
##'
##' The restriction function must return the restricted coefficients of
##' the MIDAS regression.
##'
##' @importFrom stats as.formula formula model.matrix model.response terms lsfit time
##' @importFrom zoo index index2char
##' @export
midas_r <- function(formula, data, start, Ofunction="optim", weight_gradients=NULL,...) {

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

    itr <- checkARstar(terms(eval(mf[[2]], Zenv)))
    if(!is.null(itr$lagsTable)) 
      mf[[2]] <- itr$x
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
    args <- list(...)
    y <- model.response(mf, "numeric")
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
    
    prepmd <- prepmidas_r(y,X,mt,Zenv,cl,args,start,Ofunction,weight_gradients,itr$lagsTable)
    
    prepmd <- c(prepmd, list(lhs = ysave, lhs_start = y_start, lhs_end = y_end))
    
    class(prepmd) <- "midas_r"
    
    midas_r.fit(prepmd)    
}

##' @method update midas_r
##' @importFrom stats getCall update.formula setNames
##' @export
update.midas_r <- function(object, formula.,..., evaluate = TRUE) {
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
        if(is.null(extras$start)) {
            ##If start is null, we want to fit unrestricted midas model, this means that we need to call midas_r
            call["start"] <- list(NULL)                        
            redo <- TRUE
        } else {
            cstart <- eval(call$start,object$Zenv)
            ustart[names(cstart)] <- cstart
            call$start <- ustart
            object$start_opt <- unlist(ustart)
        }
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
            midas_r.fit(object)
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
midas_r.fit <- function(x) {
    args <- x$argmap_opt
    function.opt <- args$Ofunction
    args$Ofunction <- NULL
   
    if(!(function.opt %in% c("optim","spg","optimx","lm","nls","dry_run"))) 
        stop("The optimisation function is not in the supported functions list. Please see the midasr:::midas_r.fit code for the supported function list")
    if(function.opt == "optim" | function.opt =="spg") {  
        args$par <- x$start_opt
        args$fn <- x$fn0
        if(x$use_gradient) {
            args$gr <- x$gradient
        }
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
        if(x$use_gradient) {
            args$gr <- x$gradient
        }
        opt <- try(do.call(function.opt,args),silent=TRUE)
        if(inherits(opt,"try-error")) {
            stop("The optimisation algorithm of MIDAS regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        bmet <- which.min(opt$value)
        par <- as.numeric(opt[bmet,1:length(args$par)])        
        names(par) <- names(coef(x))
        x$convergence <- opt$convcode[bmet]
    }    
    if(function.opt=="lm") {
        if(is.null(x$unrestricted))stop("Not possible to estimate MIDAS model, more parameters than observations")
        par <- coef(x$unrestricted)
        names(par) <- names(coef(x))
        opt <- NULL
        x$convergence <- 0
    }
    if(function.opt=="nls") {
        rhs <- x$rhs
        if(x$use_gradient) {
            orhs <- rhs
            rhs <- function(p) {
                res <- orhs(p)
                attr(res,"gradient") <- x$model[,-1]%*%x$gradD(p)
                res
            }
        }
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
    names(par) <- NULL
    x$midas_coefficients <- x$gen_midas_coef(par)
    x$fitted.values <- as.vector(x$model[,-1]%*%x$midas_coefficients)
    x$residuals <- as.vector(x$model[,1]-x$fitted.values)
    x
}

## Prepare necessary objects for fitting of the MIDAS regression
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
## Vaidotas Zemlys
prepmidas_r <- function(y, X, mt, Zenv, cl, args, start, Ofunction, weight_gradients, lagsTable, unrestricted = NULL, guess_start = TRUE, tau = NULL) {

    start <- start[!sapply(start,is.null)]
    if(is.null(weight_gradients)) use_gradient <- FALSE
    else use_gradient=TRUE

    if(!is.null(args$guess_start)) {
        guess_start <- args$guess_start
        args$guess_start <- NULL
    }    
    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    
    dterm <- function(fr, ltb = NULL) {
        term_name <- as.character(fr)[1]
        weight_name <- ""
        rf <- function(p)p
        grf <- function(p)diag(1)
        start <- 0
        freq <- 1
        lagstruct <- 0
        normalized <- FALSE
        if(term_name %in% c("mls", "fmls", "dmls")) {
            type <- term_name
            term_name <- as.character(fr[[2]])
            
            freq <- eval(fr[[4]], Zenv)
            lags <- eval(fr[[3]], Zenv)
            nol <- switch(type,
                          fmls = lags+1,
                          dmls = lags+1,
                          mls = length(lags)
            )
            lagstruct <- switch(type,
                                fmls = 0:lags,
                                dmls = 0:lags,
                                mls = lags
            )
            start <- rep(0, nol)
            grf <- function(p)diag(nol)
            if(length(fr) > 4 && fr[[5]] != "*") {
                mf <- fr[-5]
                mf[[1]] <- fr[[5]]
                weight_name <- as.character(fr[[5]])
                
                ##Since we allow stars and other stuff in mls, maybe allow to 
                ##specify the multiplicative property in a call to mls?
                
                if(weight_name %in% c("nealmon","nbeta","nbetaMT","gompertzp","nakagamip","lcauchyp")) {
                    normalized <- TRUE
                } else {
                    norm_attr <- attr(eval(as.name(weight_name),Zenv),"normalized")
                    if(is.null(norm_attr)) {
                        normalized <- FALSE
                    } else  {
                        if(!is.logical(norm_attr)) stop("The attribute normalized of MIDAS weight function must be logical!")
                        normalized <- norm_attr
                    }
                }
                noarg <- length(formals(eval(fr[[5]], Zenv)))
                if(noarg<2)stop("The weight function must have at least two arguments")            
                mf <- mf[1:min(length(mf), noarg + 1)]
                if(length(mf)>3) {
                    for(j in 4:length(mf)) {
                        mf[[j]] <- eval(mf[[j]], Zenv)
                    }
                }
                mf[[3]] <- ifelse(is.null(ltb), nol, sum(ltb[, 1]))
                rf <- function(p) {
                    mf[[2]] <- p
                    eval(mf,Zenv)
                }
                if(use_gradient) {
                    gmf <- mf
                    if(weight_name %in% names(weight_gradients)) {
                        weight_gradient_name <- paste0(as.character(fr[[2]]),"_tmp_gradient_fun")
                        gmf[[1]] <- as.name(weight_gradient_name)
                        assign(weight_gradient_name,weight_gradients[[weight_name]],Zenv)
                    } else {                                        
                        gmf[[1]] <- as.name(paste0(weight_name,"_gradient"))
                    }
                    grf <- function(p) {
                        gmf[[2]] <- p
                        eval(gmf,Zenv)
                    }
                } else grf <- NULL
            }
        }
        list(weight = rf,
             term_name = term_name,
             gradient = grf,
             start = start,                    
             weight_name = weight_name,
             frequency = freq,
             lag_structure = lagstruct,
             normalized = normalized
        )
    }
    if(is.null(lagsTable)){
        ltb <- rep(list(NULL), length(terms.lhs))
    } else {
        ltb <- lagsTable
        if(attr(mt,"intercept")==1) {
            ltb <- ltb[-1]   
        }
    }
    rfd <- mapply(dterm, terms.lhs, ltb, SIMPLIFY = FALSE)
    
    if (attr(mt,"intercept")==1)  {
        intc <- dterm(expression(1))
        intc$term_name <- "(Intercept)"
        rfd <- c(list(intc), rfd)
    }
    
    rf <- lapply(rfd,"[[","weight")
    names(rf) <- sapply(rfd,"[[","term_name")
    
    ##Note this is a bit of misnomer. Variable weight_names is actualy a vector of term names which have MIDAS weights.
    ##It *is not* the same as actual name of weight function. This is a left-over from the old code. 
    weight_names <- sapply(rfd,"[[","weight_name")
    weight_inds <- which(weight_names!="")
    weight_names <- names(rf)[weight_names!=""]
    
    
    start_default <- lapply(rfd,"[[","start")
    names(start_default) <- names(rf)

    ##If there are no weight functions, we have unrestricted MIDAS model.
    if(length(weight_names) == 0)Ofunction <- "lm"
    else {
        if(is.null(start)) {            
            cl$formula <- update_weights(cl$formula,setNames(lapply(1:length(weight_names), function(x)""), weight_names))
            warning("Since the start = NULL, it is assumed that U-MIDAS model is fitted")
            return(eval(cl,Zenv))            
        } else {
            if(any(!weight_names%in% names(start)))stop("Starting values for weight parameters must be supplied")
        }
    }
    
    start_default[names(start)] <- start
            
    np <- cumsum(sapply(start_default,length))

    build_indices <- function(ci,nm) {
        inds <- cbind(c(1,ci[-length(ci)]+1),ci)
        inds <- apply(inds,1,function(x)list(x[1]:x[2]))
        inds <- lapply(inds,function(x)x[[1]])
        names(inds) <- nm
        inds
    }
    
    pinds <- build_indices(np,names(start_default))

    for(i in 1:length(start_default))names(start_default[[i]]) <- NULL

    
    if(!is.null(lagsTable)) {
        inones <- function(ones, intro) {
            ones[ones == 1] <- intro
            ones
        }
        yname <- all.vars(mt[[2]])
        nms <- names(pinds)
        all_coef2 <- function(p) {
            pp <- lapply(pinds, function(x) p[x])
            cr <- c(1, -p[pinds[[yname]]])
            res <- mapply(function(fun, cf, tb) {
                restr <- fun(cf)
                if(is.null(tb)) {
                    restr
                } else {
                    mltp <- rowSums(apply(tb, 2, inones, restr) %*% diag(cr))
                    mltp[rowSums(tb)!= 0]
                }
            }, rf, pp, lagsTable, SIMPLIFY = FALSE)
            return(res)
        }
    } else {    
        all_coef2 <- function(p) {              
            pp <- lapply(pinds,function(x)p[x])     
            res <- mapply(function(fun,param)fun(param),rf,pp,SIMPLIFY=FALSE)
            return(res)
        }
    }
    
    
    initial_midas_coef <- all_coef2(unlist(start_default)) 

    if(sum(is.na(unlist(initial_midas_coef)))>0) stop("Check your starting values, NA in midas coefficients") 
    
    npx <- cumsum(sapply(initial_midas_coef,length))
    xinds <- build_indices(npx,names(start_default))        
     
    if(length(weight_names)>0 && guess_start) {
        wi <- rep(FALSE,length(rf))
        wi[weight_inds] <- TRUE
        Xstart <- mapply(function(cf,inds,iswhgt) {        
            if(iswhgt) {
                X[, inds, drop = FALSE] %*% cf
            }
            else X[, inds, drop = FALSE]
        }, initial_midas_coef, xinds,wi,SIMPLIFY=FALSE)

        npxx <- cumsum(sapply(Xstart,function(x) {
            ifelse(is.null(dim(x)),1,ncol(x))
            }))
        xxinds <- build_indices(npxx,names(start_default))
        XX <- do.call("cbind",Xstart)
        ###If the starting values for the weight restriction are all zeros, then the weighted explanatory variable is zero.
        ###In this case lsfit gives a warning about colinear matrix, which we can ignore.
        prec <- suppressWarnings(lsfit(XX,y,intercept=FALSE))
        lmstart <- lapply(xxinds,function(x)coef(prec)[x])
        names(lmstart) <- names(xxinds)
        for(i in 1:length(lmstart))names(lmstart[[i]]) <- NULL

        nms <- !(names(start_default) %in% names(start))
        start_default[nms] <- lmstart[nms]
        
        
        term_norm_info <- sapply(rfd, "[[", "normalized")
        names(term_norm_info) <- names(rf) 
        if(sum(term_norm_info)>0) {
            for(nm in names(term_norm_info)[term_norm_info]) {
                if(abs(start_default[[nm]][1]-1)<1e-5) {
                    start_default[[nm]][1] <- lmstart[nm]
                } else {
                    warning("For normalized weights please supply 1 as a first coefficient. Leaving the starting value unchanged.")
                }
            }
        }
    }
    
    starto <- unlist(start_default)
    ##This is workaround for AR* model
    all_coef <- function(p) unlist(all_coef2(p)) 
    mdsrhs <- function(p) {       
        coefs <- all_coef(p)
        X%*%coefs
    }
  
    #aa <- try(mdsrhs(starto))
    
    fn0 <- function(p,...) {
        r <- y - mdsrhs(p)
        sum(r^2)
    }
    if(!is.null(tau)) {
        fn0 <- function(p,...) {
            r <- y - mdsrhs(p)
            sum(tau * pmax(r, 0) + (tau - 1) * pmin(r,0))
        }
    }
    
    if(!use_gradient) {
        gradD <- function(p)jacobian(all_coef,p)
        gr <- function(p)grad(fn0,p)
    }
    else {
        grf <- sapply(rfd,"[[","gradient")
        ##Calculate the initial value to get the idea about the dimensions
        pp0 <- lapply(pinds,function(xx)starto[xx])            
        grmat0 <- mapply(function(fun,param)fun(param),grf,pp0,SIMPLIFY=FALSE)
        colnos <- sapply(grmat0,ncol)
        rownos <- sapply(grmat0,nrow)
        np <- length(colnos)
        ccol <- cumsum(colnos)
        rrow <- cumsum(rownos)
        pindm <- cbind(c(1,rrow[-np]+1),rrow,
                       c(1,ccol[-np]+1),ccol)
        pindm <- apply(pindm,1,function(x)list(row=x[1]:x[2],col=x[3]:x[4]))
        if(is.null(lagsTable)) {
          gradD <- function(p) {
            pp <- lapply(pinds,function(x)p[x])
            grmat <- mapply(function(fun,param)fun(param),grf,pp,SIMPLIFY=FALSE)
            if(length(grmat)==1) {
              res <- grmat[[1]]
            }
            else {
              res <- matrix(0,nrow=sum(rownos),ncol=sum(colnos))
              for(j in 1:length(grmat)) {
                ind <- pindm[[j]]
                res[ind$row,ind$col] <- grmat[[j]]                    
              }
            }
            res
          }
        } else {
            expandD <- function(grm, ltb, cr) {
                if(is.null(ltb)) {
                    return(grm)
                } else {        
                    el <- lapply(data.frame(ltb), inones2, grm)
                    mltp <- Reduce("+", mapply(`*`, el, cr, SIMPLIFY = FALSE))
                    return(mltp[rowSums(ltb)!= 0, ])
                } 
            }
            
            inones2 <- function(ones, intro) {
                m <- matrix(0, nrow = length(ones), ncol = ncol(intro))
                if(sum(ones) != nrow(intro)) stop("Wrong gradient for AR* term")
                m[ones == 1, ] <- intro
                m
            }
            
            expandD2 <- function(fun, param, ltb, nparam = 1){
                cf <- fun(param)
                if(is.null(ltb)) return(matrix(0, nrow = length(cf), ncol = nparam))
                else {
                    mltp <- -apply(ltb, 2, inones, cf) 
                    return(mltp[rowSums(ltb) != 0, -1, drop = FALSE])
                }    
            }
            dind <- which(names(pinds)==yname)
            cr <- c(1, -starto[pinds[[dind]]])
            pp <- lapply(pinds, function(x) starto[x])
            grmat1 <- mapply(function(fun, param) fun(param), grf, pp, SIMPLIFY = FALSE)
            egrmat1 <- mapply(expandD, grmat1, lagsTable, SIMPLIFY = FALSE, MoreArgs = list(cr))
            colnos <- sapply(egrmat1,ncol)
            rownos <- sapply(egrmat1,nrow)
            np <- length(colnos)
            ccol <- cumsum(colnos)
            rrow <- cumsum(rownos)
            pindm <- cbind(c(1,rrow[-np]+1),rrow,
                           c(1,ccol[-np]+1),ccol)
            pindm <- apply(pindm,1,function(x)list(row=x[1]:x[2],col=x[3]:x[4]))
            
            gradD <- function(p) {   
                cr <- c(1, -p[pinds[[dind]]])
                pp <- lapply(pinds, function(x) p[x])
                grmat <- mapply(function(fun, param) fun(param), grf, pp, SIMPLIFY = FALSE)
                egrmat <- mapply(expandD, grmat, lagsTable, SIMPLIFY = FALSE, MoreArgs = list(cr))
                res <- matrix(0,nrow=sum(rownos),ncol=sum(colnos))
                gr_star <- do.call("rbind",mapply(expandD2, rf, pp , lagsTable, SIMPLIFY = FALSE, MoreArgs = list(length(pinds[[dind]]))))
                res[, pinds[[dind]]] <- gr_star
                for(j in 1:length(egrmat)) {
                    ind <- pindm[[j]]
                    res[ind$row,ind$col] <- egrmat[[j]]
                }
                res 
            }
        }
        gr <- function(p) {
             XD <- X%*%gradD(p)
             resid <- y - X %*% all_coef(p)
             as.vector(-2*apply(as.vector(resid)*XD,2,sum))             
        }
        ##Seems to work
    }
    hess <- function(x)numDeriv::hessian(fn0,x)
      
    if(is.null(unrestricted)) {        
        if(ncol(X)<nrow(X)) {
            if(attr(mt,"intercept")==1) {
                unrestricted <- lm(y~.,data=data.frame(cbind(y,X[,-1]),check.names=FALSE))
            } else {
                unrestricted <- lm(y~.-1,data=data.frame(cbind(y,X),check.names=FALSE))
            }
            
        }
    }

    control <- c(list(Ofunction=Ofunction),args)
    ##Override default method of optim. Use BFGS instead of Nelder-Mead
    if(!("method"%in% names(control)) & Ofunction=="optim") {        
        control$method <- "BFGS"
    }    
    term_info <- rfd
    names(term_info) <- sapply(term_info,"[[","term_name")
    term_info <- mapply(function(term,pind,xind){
        term$start <- NULL
        term$coef_index <- pind
        term$midas_coef_index <- xind
        term
    },term_info,pinds[names(term_info)],xinds[names(term_info)],SIMPLIFY=FALSE)
    
    if(!is.null(tau))  {
        ##At the moment do not calculate the gradient and hessian for 
        ##quantile regression, as it does not make sense
        gr <- NULL
        hess <- NULL
    }
    list(coefficients=starto,
         midas_coefficients=all_coef(starto),
         model=cbind(y,X),         
         unrestricted=unrestricted,
         term_info=term_info,
         fn0=fn0,
         rhs=mdsrhs,
         gen_midas_coef = all_coef,
         opt=NULL,
         argmap_opt=control,
         start_opt=starto,
         start_list=start,
         call=cl,
         terms=mt,
         gradient=gr,
         hessian=hess,
         gradD=gradD,
         Zenv=Zenv,
         use_gradient=use_gradient,
         nobs=nrow(X))   
}

##' Restricted MIDAS regression
##'
##' Function for fitting MIDAS regression without the formula interface
##' @param y model response
##' @param X prepared matrix of high frequency variable lags
##' @param z additional low frequency variables
##' @param weight the weight function
##' @param grw the gradient of weight function
##' @param startx the starting values for weight function
##' @param startz the starting values for additional low frequency variables
##' @param method a method passed to \link{optimx}
##' @param ... additional parameters to \link{optimx}
##' @return an object similar to \code{midas_r} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @import numDeriv
##' @import optimx
##' @importFrom stats na.omit
##' @examples
##' 
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##'
##' X<-fmls(x,11,12)
##'
##' midas_r_simple(y,X,trend,weight=nealmon,startx=c(0,0,0))
##' @export
##' 
midas_r_simple <- function(y,X,z=NULL,weight,grw=NULL,startx,startz=NULL,method=c("Nelder-Mead","BFGS"),...) {
    d <- ncol(X)
    nw <- length(startx)
    if(!is.matrix(z))z <- matrix(z,ncol=1)
    model <- na.omit(cbind(y,X,z))
    y <- model[,1]
    XX <- model[,-1]
    
    if(is.null(z)) {        
        all_coef <- function(p) {
            weight(p,d)
        }
        gradD <- function(p)grw(p,d)
        start <- startx
    }
    else {
        all_coef <- function(p) {
            c(weight(p[1:nw],d),p[-nw:-1])
        }        
        nz <- ncol(z)
        if(is.null(startz)) {            
            ZZ <- model[,1+1:d]%*%weight(startx,d)
            Z <- model[,(d+2):ncol(model)]
            prec <- suppressWarnings(lsfit(cbind(Z,ZZ),y,intercept=FALSE))
            startz <- coef(prec)[1:nz]
        }
        if(!is.null(grw)) {
            gradD <- function(p) {
                ww <- grw(p[1:nw],d)               
                zr <- matrix(0,nrow=d,ncol=nz)
                zb <- matrix(0,nrow=nz,ncol=nw)
                rbind(cbind(ww,zr),cbind(zb,diag(nz)))
            }
        }
        else gradD <- NULL
        start <- c(startx,startz)
    }
   
    n <- nrow(model)
    fn0 <- function(p) {
        sum((y-XX%*%all_coef(p))^2)
    }
  
    if(is.null(grw)) {
        gradD <- function(p)jacobian(all_coef,p)
        gr <- function(p)grad(fn0,p)
        gr0 <- NULL
    }
    else {
        gr <- function(p) {
             XD <- XX %*% gradD(p)
             resid <- y - XX %*% all_coef(p)
             as.vector(-2*apply(as.vector(resid)*XD,2,sum)) 
        }
        gr0 <- gr
    }
    opt <- optimx(start,fn0,gr0,method=method,...)
    bmet <- which.min(opt$value)
    par <- as.numeric(opt[bmet, 1:length(start)])   
    call <- match.call()
    fitted.values <- as.vector(XX%*%all_coef(par))
    list(coefficients=par,
         midas_coefficients=all_coef(par),
         model=model,
         weights=weight,
         fn0=fn0,    
         opt=opt,
         call=call,
         gradient=gr,
         hessian=function(x)numDeriv::hessian(fn0,x),
         gradD=gradD,
         fitted.values=fitted.values,
         residuals=as.vector(y-fitted.values))
             
}
##' Updates weights in a expression with MIDAS term
##'
##' For a MIDAS term \code{fmls(x, 6, 1, nealmon)} change weight \code{nealmon} to another weight.
##' @title Updates weights in MIDAS regression formula
##' @param expr expression with MIDAS term
##' @param tb a named list with redefined weights
##' @return an expression with changed weights
##' @author Vaidotas Zemlys
##' @export
##' @examples
##'
##' update_weights(y~trend+mls(x,0:7,4,nealmon)+mls(z,0:16,12,nealmon),list(x = "nbeta", z = ""))
##' 
update_weights <- function(expr,tb) {
    if(length(expr)==3) {
        expr[[2]] <- update_weights(expr[[2]],tb)
        expr[[3]] <- update_weights(expr[[3]],tb)
    }
    if(length(expr)==5) {
        fun <- as.character(expr[[1]])
        if(fun[[1]] %in% c("fmls","mls","dmls")) {
            term_name <- as.character(expr[[2]])
            if(term_name %in% names(tb)) {
                if(is.null(tb[[term_name]])|| tb[[term_name]] == "") {
                    expr <- expr[1:4]
                } else expr[[5]] <- as.name(tb[[term_name]])                    
            }
        }
        else return(expr)
    }
    return(expr)
}

## Check whether the MIDAS model is MIDAS-AR* model
##
## authored by Julius Vainora
checkARstar <- function(trms) {
  vars <- as.list(attr(trms, "variables"))[-2:-1]
  env <- environment(trms)
  idx <- which(sapply(vars, function(y) if(length(y) >= 2) y[[2]]) == trms[[2]])
  
  lagsTable <- NULL
  if(length(idx) > 0 && length(vars[[idx]]) >= 5 && vars[[idx]][[5]] == "*") {    
    fs <- lapply(sapply(vars, function(y) if(length(y) >= 4) y[[4]]), eval, env)
    if(length(unique(unlist(fs))) > 1) {
      ## mls for y is assumed
      lags <- eval(vars[[idx]][[3]], env)
      push <- lapply(fs, "*", lags)
      
      lagsTable <- lapply(1:length(vars), function(w) {
        z <- vars[[w]]
        if(length(z) >= 4 && eval(z[[4]], env) != 1) {
          l <- eval(z[[3]], env)
          if(length(l) == 1 & as.character(z[1]) %in% c("fmls", "dmls"))
            l <- 0:l
          tp <- matrix(0, ncol = length(lags) + 1, nrow = max(l) + max(push[[w]]) + 1)
          tp[l + 1, 1] <- 1
          for(r in 2:ncol(tp))
            tp[l + push[[w]][r - 1] + 1, r] <- 1
          tp
        }
      })
      
      shortSeq <- function(s) {
        wt <- which(!diff(s) == 1)
        idx <- c(1, 1 + c(wt, wt - 1), length(s))
        ams <- s[intersect(1:length(s), idx)]
        fc <- cumsum(c(TRUE, !round(diff(ams) / 2 + head(ams, -1)) %in% s))
        out <- lapply(split(ams, fc), function(x) if(length(x) == 2) 
          do.call("call", c(":", as.list(x))) else x)
        names(out) <- NULL; out
      }
      
      vars <- lapply(1:length(vars), function(w) {
        z <- vars[[w]]
        
        if(length(z) >= 4 && eval(z[[4]], env) != 1) {
          fun <- as.character(z[1])
          l <- eval(z[[3]], env)
          if(fun %in% c("fmls", "dmls")) {
            if(length(l) == 1)
              l <- 0:l
            else
              stop("fmls and dmls are not used with a vector of lag orders")                    
          }          
          nl <- sort(unique(l + rep(c(0, push[[w]]), each = length(l))))
          if(fun == "mls")
            z[[3]] <- do.call("call", c("c", shortSeq(nl)))
          else if(all(diff(nl) == 1))
            z[3] <- max(nl)
          else if(fun == "fmls"){
            z[1] <- call("mls")
            z[[3]] <- do.call("call", c("c", shortSeq(nl)))
          } else
            # Problem in case of dmls and not full lag vector
            stop("Use fmls or mls instead of dmls")          
        }; z
      })
      icp <- attr(trms, "intercept") == 1
      trms <- formula(paste(trms[[2]], "~", paste(vars, collapse = " + ")), env)
      if(!icp)
        trms <- update.formula(trms, . ~ . -1)
      else
        lagsTable <- c(list(NULL), lagsTable)
    }
  }
  list(x = trms, lagsTable = lagsTable)
}
