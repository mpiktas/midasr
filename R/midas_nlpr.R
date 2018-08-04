##' Non-linear parametric MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares.
##'
##' @param formula formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param data a named list containing data with mixed frequencies
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are the arguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
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
##' MIDAS-AR* (a model with a common factor, see (Clements and Galvao, 2008)) can be estimated by specifying additional argument, see an example.
##'
##' The restriction function must return the restricted coefficients of
##' the MIDAS regression.
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
    
    prepmd <- prep_midas_nlpr(y,X,mt,Zenv,cl,args,start,Ofunction,weight_gradients,itr$lagsTable)
    
    prepmd <- c(prepmd, list(lhs = ysave, lhs_start = y_start, lhs_end = y_end))
    
    class(prepmd) <- "midas_nlpr"
    
    midas_r.fit(prepmd)    
}

##' @method update midas_r
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
midas_nlpr.fit <- function(x) {
    args <- x$argmap_opt
    function.opt <- args$Ofunction
    args$Ofunction <- NULL
   
    if(!(function.opt %in% c("optim","spg","optimx","nls","dry_run"))) 
        stop("The optimisation function is not in the supported functions list. Please see the midasr:::midas_nlpr.fit code for the supported function list")
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
prep_midas_nlpr <- function(y, X, mt, Zenv, cl, args, start, Ofunction,  guess_start = TRUE) {

    start <- start[!sapply(start,is.null)]
    
    if(!is.null(args$guess_start)) {
        guess_start <- args$guess_start
        args$guess_start <- NULL
    }    
    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    
    dterm <- function(fr, ltb = NULL) {
        term_name <- as.character(fr)[1]
        weight_name <- ""
        rf <- function(p)p
        start <- 0
        freq <- 1
        lagstruct <- 0
        if(term_name %in% c("mls", "fmls", "dmls","mlsd")) {
            type <- term_name
            term_name <- as.character(fr[[2]])
            
            if(type == "mlsd") {
                wpos <- 6
                freq <- NA
            } else {
                wpos <- 5
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
                if(wpos == 6) {
                   # We need to exclude date information here
                   mf <- fr[-((wpos - 1):wpos)]
                   mf[[4]] <- NA
                } else {
                    mf <- fr[-wpos]
                }
                mf[[1]] <- fr[[wpos]]
                weight_name <- as.character(fr[[wpos]])
                
                noarg <- length(formals(eval(fr[[wpos]], Zenv)))
                if(noarg < 2) stop("The weight function must have at least two arguments")            
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
            }
        }
        list(weight = rf,
             term_name = term_name,
             start = start,                    
             weight_name = weight_name,
             frequency = freq,
             lag_structure = lagstruct
        )
    }
    
        
    ltb <- rep(list(NULL), length(terms.lhs))
    
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

    if(any(!weight_names%in% names(start)))stop("Starting values for weight parameters must be supplied")
        
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

    
    all_coef2 <- function(p) {              
            pp <- lapply(pinds,function(x)p[x])     
            res <- mapply(function(fun,param)fun(param),rf,pp,SIMPLIFY=FALSE)
            return(res)
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
        
        for(ww in which(wi)) {
            normalized <- FALSE
            if(rfd[[ww]]$weight_name %in% c("nealmon","nbeta","nbetaMT","gompertzp","nakagamip","lcauchyp")) normalized <- TRUE
            else {
                normalized <- is_weight_normalized(rf[[ww]], start_default[[ww]])
            }
            if(normalized) {
                start_default[[ww]][1] <- lmstart[[ww]]
            }
        }
    }
    
    starto <- unlist(start_default)
    
    all_coef <- function(p) unlist(all_coef2(p)) 
    
    mdsrhs <- function(p) {       
        coefs <- all_coef(p)
        X%*%coefs
    }
  
    fn0 <- function(p,...) {
        r <- y - mdsrhs(p)
        sum(r^2)
    }
    
    hess <- function(x)numDeriv::hessian(fn0,x)
      
    
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
    
    
    list(coefficients=starto,
         midas_coefficients=all_coef(starto),
         model=cbind(y,X),         
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

##' Non-linear parametric  MIDAS regression
##'
##' Function for fitting non-linear parametric MIDAS regression without the formula interface
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
##'
##' @export
##' 
midas_nlpr_simple <- function(y,X,z=NULL,weight,startx,startz=NULL,method=c("Nelder-Mead","BFGS"),...) {
    d <- ncol(X)
    nw <- length(startx)
   
    if(!is.null(z) && !is.matrix(z)) z <- matrix(z,ncol=1)
    model <- na.omit(cbind(y,X,z))
    y <- model[,1]
    XX <- model[,-1]
    
    if(is.null(z)) {        
        all_coef <- function(p) {
            weight(p,d)
        }
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
        start <- c(startx,startz)
    }
   
    n <- nrow(model)
    fn0 <- function(p) {
        sum((y-XX%*%all_coef(p))^2)
    }
  
    
    gradD <- function(p)jacobian(all_coef,p)
    gr <- function(p)grad(fn0,p)
    gr0 <- NULL
    
   
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
