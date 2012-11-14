##' Estimate unrestricted MIDAS regression
##'
##' Estimate unrestricted MIDAS regression using OLS. This function is a wrapper for \code{lm}.
##' 
##' @param formula MIDAS regression model formula
##' @param ldata low frequency data
##' @param hdata high frequency data
##' @param ... further arguments, which could be passed to \code{\link{lm}} function.
##' @return \code{\link{lm}} object.
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Kvedaras V., Zemlys, V. \emph{Testing the functional constraints on parameters in regressions with variables of different frequency} Economics Letters 116 (2012) 250-254 
##' @examples
##' ##The parameter function
##' theta.h0 <- function(p, dk) {
##'    i <- (1:dk-1)/100
##'    pol <- p[3]*i + p[4]*i^2
##'    (p[1] + p[2]*i)*exp(pol)
##' }
##'
##' ##Generate coefficients
##' theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)
##'
##' ##Plot the coefficients
##' plot(theta0)
##'
##' ##Generate the predictor variable
##' x <- simplearma.sim(list(ar=0.6),1500*12,1,12)
##'
##' ##Simulate the response variable
##' y <- midas.sim(500,theta0,x,1)
##'
##' ##Create low frequency data.frame
##' ldt <- data.frame(y=y,trend=1:length(y))
##'
##' ##Create high frequency data.frame
##'
##' hdt <- data.frame(x=window(x, start=start(y)))
##' 
##' ##Fit unrestricted model
##' mu <- midas_u(y~embedlf(x,3,12)-1, ldt, hdt)
##'
##' ##Include intercept and trend in regression
##'
##' mu.it <- midas_u(y~embedlf(x,3,12)+trend, ldt, hdt)
##' 
##' @details MIDAS regression has the following form:
##' 
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+\mathbf{z_t}\mathbf{\beta}+u_t}
##'
##' or alternatively
##'
##' \deqn{y_t=\sum_{h=0}^{(k+1)m}\theta_hx_{tm-h}+\mathbf{z_t}\mathbf{\beta}+u_t,}
##' where \eqn{m} is the frequency ratio and
##' \eqn{k} is the number of lags included in the regression. 
##'
##' Given certain assumptions the coefficients can be estimated using usual OLS and they have the familiar properties associated with simple linear regression.
##'
##' MIDAS regression involves times series with different frequencies.
##' 
##' @export
midas_u <- function(formula, ldata=NULL, hdata=NULL,...) {
    Zenv <- new.env(parent=environment(formula))

    if(missing(ldata) | missing(hdata)) {
        ee <- NULL
    }
    else {
        data <- check_mixfreq(ldata,hdata)

        ee <- as.environment(c(as.list(data$lowfreq),as.list(data$highfreq)))
        parent.env(ee) <- parent.frame()
    }
    
    assign("ee",ee,Zenv)    
    mf <- match.call(expand.dots = TRUE)
    mf <- mf[-4]
    mf[[1L]] <- as.name("lm")
    mf[[3L]] <- as.name("ee")   
    names(mf)[3] <- "data"    
    eval(mf,Zenv)
}
##' Restricted MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares.
##'
##' @param x either formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{embedlf}} function
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param control the list with information which R function to use for optimisation. The list must have element named \code{Rfunction} which contains character string of chosen R function. Other elements of the list are thearguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param ... additional arguments supplied to restriction function(s)
##' @return a \code{midas_r} object which is the list with the following elements:
##' 
##' \item{coefficients}{the estimates of parameters of restrictions}
##' \item{midas.coefficientas}{the estimates of restricted coefficients of MIDAS regression}
##' \item{model}{model data}
##' \item{restrictions}{the restriction function(s) used in estimation.}
##' \item{unrestricted}{unrestricted regression estimated using \code{\link{midas_u}}}
##' \item{param.map}{parameter map for optimisation function}
##' \item{fn0}{optimisation function for non-linear least squares problem solved in restricted MIDAS regression}
##' \item{rhs}{the function which evaluates the right-hand side of the MIDAS regression}
##' \item{allcoef}{the function which evaluates the restricted coefficientsof MIDAS regression}
##' \item{opt}{the output of optimisation procedure}
##' \item{argmap.opt}{the list containing the name of optimisation function together with arguments for optimisation function}
##' \item{start.opt}{the starting values used in optimisation}
##' 
##' 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname midas_r
##' @seealso midas_r.midas_r
##' @examples
##' ##The parameter function
##' theta.h0 <- function(p, dk) {
##'    i <- (1:dk-1)/100
##'    pol <- p[3]*i + p[4]*i^2
##'    (p[1] + p[2]*i)*exp(pol)
##' }
##'
##' ##Generate coefficients
##' theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)
##'
##' ##Plot the coefficients
##' plot(theta0)
##'
##' ##Generate the predictor variable
##' x <- simplearma.sim(list(ar=0.6),1500*12,1,12)
##'
##' ##Simulate the response variable
##' y <- midas.sim(500,theta0,x,1)
##'
##' ##Remove unnecessary history of x
##' x <- window(x,start=start(y))
##' 
##' ##Fit restricted model
##' mr <- midas_r(y~embedlf(x,4*12,12,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)
##'
##' ##Include intercept and trend in regression
##'
##' mr.it <- midas_r(y~embedlf(x,4*12,12,theta.h0)+trend,data.frame(y=y,trend=1:500),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)
##' 
##' @details Given MIDAS regression:
##'
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+\mathbf{z_t}\beta+u_t}
##'
##' estimate the parameters of the restriction
##'
##' \deqn{\theta_h=g(h,\lambda),}
##' where \eqn{h=0,...,(k+1)m}, together with coefficients \eqn{\beta} corresponding to additional low frequency regressors.
##'
##' MIDAS regression involves times series with different frequencies. 
##'
##' The restriction function must return the restricted coefficients of
##' the MIDAS regression.
##' 
##' @export
midas_r <- function(x,...)UseMethod("midas_r")

is.midas_r <- function(x) inherits(x,"midas_r")

#' @rdname midas_r
#' @method midas_r default
#' @export
midas_r.default <- function(x, ldata=NULL, hdata=NULL, start, control=list(Rfunction="optim",method="BFGS"), ...) {

    Zenv <- new.env(parent=environment(x))
      
    if(missing(ldata)|missing(hdata)) {
        ee <- NULL
    }
    else {
        data <- check_mixfreq(ldata,hdata)

        ee <- as.environment(c(as.list(data$lowfreq),as.list(data$highfreq)))
        parent.env(ee) <- parent.frame()
    }
    assign("ee",ee,Zenv)
    cll <- match.call()
    mf <- match.call(expand.dots = FALSE)
    ##Fix this!!
    m <- match(c("x", "ldata"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- as.name("na.omit")
    names(mf)[c(2,3,4)] <- c("formula","data","na.action")
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
  
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)    
    
    cl <- list(...)

    ##High frequency variables can enter to formula
    ##only within embedlf function
    terms.lhs <- attr(mt,"term.labels") 
    if(attr(mt,"intercept")==1) terms.lhs <- c("(Intercept)",terms.lhs)
  
    rf <- lapply(terms.lhs,function(fr) {     
        if(length(grep("embedlf",fr))>0) {
            fr <- parse(text=fr)[[1]]
            ff <- eval(fr[[5]],parent.frame())
            rf.arg <- formals(ff)
            class(rf.arg) <- "list"
            rf.argnm <-  intersect(names(rf.arg)[-1],names(cl))
            rf.arg[rf.argnm] <- cl[rf.argnm]            
            rf.name <- as.character(fr[[5]])
            rf <- function(p) {
                rf.arg[[1]] <- p
                do.call(rf.name,rf.arg)
            }
            return(rf)
        }
        else {
            return(function(p)p)
        }
    })

    names(rf) <- sapply(terms.lhs,function(fr) {
        if(length(grep("embedlf",fr))>0) {
            fr <- parse(text=fr)[[1]]
            as.character(fr[[5]])
        }
        else {
            fr
        }
    })
    
    restr.name <- setdiff(names(rf), terms.lhs)
    
    start_default<- as.list(rep(0,length(rf)))
    names(start_default) <- names(rf)

    if(any(!restr.name%in% names(start)))stop("Starting values are not supplied for restriction function(s)")
    
    start_default[names(start)] <- start

    restr.name <- setdiff(names(start_default), terms.lhs)
    
    restr.no <- sum(sapply(start_default[restr.name], length))
    
    np <- cumsum(sapply(start_default,length))
    
    pinds <- cbind(c(1,np[-length(np)]+1),np)
    pinds <- apply(pinds,1,function(x)list(x[1]:x[2]))
    pinds <- lapply(pinds,function(x)x[[1]])
    names(pinds) <- names(start_default)

    for(i in 1:length(start_default))names(start_default[[i]]) <- NULL
    
    starto <- unlist(start_default)
   
    all_coef <- function(p) {              
        pp <- lapply(pinds,function(x)p[x])     
        res <- mapply(function(fun,param)fun(param),rf,pp,SIMPLIFY=FALSE)
        unlist(res)
    }   
    
    mdsrhs <- function(p) {       
        coefs <- all_coef(p)
        X%*%coefs
    }
  
    aa <- try(mdsrhs(starto))
    
    fn0 <- function(p,...) {
        r <- y - mdsrhs(p)
        sum(r^2)
    }
    

    unrestricted <- NULL
    if(ncol(X)<nrow(X)) unrestricted <- lm(y~.-1,data=data.frame(cbind(y,X),check.names=FALSE))
    
    prepmd <- list(coefficients=starto,
                   midas.coefficients=all_coef(starto),
                   model=cbind(y,X),
                   restrictions=rf[restr.name],
                   unrestricted=unrestricted,
                   param.map=pinds,
                   fn0=fn0,
                   rhs=mdsrhs,
                   allcoef=all_coef,
                   opt=NULL,
                   argmap.opt=control,
                   start.opt=starto,
                   call=cll,
                   terms=mt)
    class(prepmd) <- "midas_r"
    midas_r.fit(prepmd)    
}

##' Restricted MIDAS regression
##'
##' Reestimate the MIDAS regression with different starting values
##' 
##' @param x \code{midas_r} object 
##' @param start the starting values
##' @param Rfunction a character string of the optimisation function to use. The default value is to use the function of previous optimisation.
##' @param ... further arguments to optimisation function. If none are supplied, the arguments of previous optimisation are used.
##' @return \code{midas_r} object
##' @method midas_r midas_r
##' @seealso midas_r
##' @author Vaidotas Zemlys
##' @export
midas_r.midas_r <- function(x,start=coef(x),Rfunction=x$argmap.opt$Rfunction,...) {
   
    oarg <- list(...)
    cl <- match.call()
    dotargnm <- names(oarg)
    
    ##Perform check whether arguments are ok and eval them
    if(length(dotargnm)>0) {
        offending <- dotargnm[!dotargnm %in% names(formals(Rfunction))]
        if(length(offending)>0)  {
            stop(paste("The function ",Rfunction," does not have the following arguments: ", paste(offending,collapse=", "),sep=""))
        }
    }
    else {
        oarg <- NULL
    }

    if(Rfunction!=x$argmap.opt$Rfunction) {
        argmap <- c(list(Rfunction=Rfunction),oarg)
    }              
    else {
         argmap <- x$argmap.opt
         argmap$func <- NULL
         argnm <- union(names(argmap),names(oarg))
         marg <- vector("list",length(argnm))
         names(marg) <- argnm
         ##New supplied arguments override the old ones
         marg[names(oarg)] <- oarg
         ##Already set arguments are left intact
         oldarg <- setdiff(names(argmap),names(oarg))
         marg[oldarg] <- argmap[oldarg]
         argmap <- c(list(Rfunction=Rfunction),marg)
    }
    
    x$start.opt <- start
    x$argmap.opt <- argmap
    x$call <- cl
    midas_r.fit(x)
}

##' Fit restricted MIDAS regression
##'
##' Workhorse function for fitting restricted MIDAS regression
##'  
##' @param x \code{midas_r} object
##' @return \code{\link{midas_r}} object
##' @author Vaidotas Zemlys
midas_r.fit <- function(x) {
    args <- x$argmap.opt
    function.opt <- args$Rfunction
    args$Rfunction <- NULL
    if(function.opt=="optim" | function.opt=="spg") {  
        args$p <- x$start.opt
        args$fn <- x$fn0
        opt <- do.call(function.opt,args)
        par <- opt$par
        names(par) <- names(coef(x))
    }
    if(function.opt=="nls") {
        rhs <- x$rhs
        y <- x$model[,1]
        args$formula <- formula(y~rhs(p))
        args$start <- list(p=x$start.opt)
        opt <- do.call("nls",args)
        par <- coef(opt)
        names(par) <- names(coef(x))
    }
    x$opt <- opt
    x$coefficients <- par
    x$midas.coefficients <- x$allcoef(par)
    x$fitted.values <- as.vector(x$model[,-1]%*%x$midas.coefficients)
    x$residuals <- as.vector(x$model[,1]-x$fitted.values)
    x
}

##' MIDAS regression model deviance
##'
##' Returns the deviance of a fitted MIDAS regression object
##' @param object a \code{\link{midas_r}} object
##' @param ... currently nothing
##' @return The sum of squared residuals
##' @author Virmantas Kvedaras.Vaidotas Zemlys
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

summary.midas_r <- function(x) {
    r <- as.vector(residuals(x))
    param <- coef(x)
    pnames <- names(param)
    n <- length(r)
    p <- length(param)
    rdf <- n - p
    resvar <- deviance(x)/rdf
    XD <- x$model[,-1]%*%gradD(x)(coef(x))
    R <- qr.R(qr(XD))
    XDtXDinv <- chol2inv(R)
    dimnames(XDtXDinv) <- list(pnames,pnames)
    se <- sqrt(diag(XDtXDinv)*resvar)
    tval <- param/se
    param <- cbind(param,se,tval,2*pt(abs(tval),rdf,lower.tail=FALSE))
    dimnames(param) <- list(pnames, c("Estimate", "Std. Error", 
        "t value", "Pr(>|t|)"))
    ans <- list(formula=formula(x$terms),residuals=r,sigma=sqrt(resvar),
                df=c(p,rdf), cov.unscaled=XDtXDinv, call=x$call,
                coefficients=param,midas.coefficients=midas_coef(x))
    class(ans) <- "summary.midas_r"
    ans
}

print.summary.midas_r <- function(x, digits=max(3, getOption("digits") - 3 ), signif.stars = getOption("show.signif.stars"), ...) {
    cat("\n Formula", deparse(formula(x)),"\n")
    df <- x$df
    rdf <- df[2L]
    cat("\n Parameters:\n")
    printCoefmat(x$coefficients,digits=digits,signif.stars=signif.stars,...)
    cat("\n Residual standard error:", format(signif(x$sigma,digits)), "on", rdf , "degrees of freedom\n")
    invisible(x)
}

print.midas_r <- function(x, digits=max(3,getOption("digits")-3),...) {
    cat("MIDAS regression model\n")
    cat(" model:", deparse(formula(x)),"\n")
    print(coef(x),digits = digits, ...)
    cat("\n")
    cat("Function ", x$argmap.opt$func, " was used for fitting\n")
    invisible(x)
}


##' Calculate numerical approximation of gradient of RSS of
##' MIDAS regression 
##'
##' Calculate nummerical approximation of gradient of residual sum of
##' squares function of estimated MIDAS regression
##' 
##' @param x a \code{\link{midas_r}} object
##' @param ... further arguments to function \code{\link{grad}} from package \code{numDeriv}
##' @return numeric vector
##' @author Vaidotas Zemlys
##' @seealso midas_r.formula
##' @import numDeriv
##' @rdname gradient
##' @export
gradient <- function(x,...)UseMethod("gradient")

#' @rdname gradient
#' @method gradient midas_r
#' @export
gradient.midas_r <- function(x,...) {
    grad(x$fn0,coef(x),...)
}

##' Calculate numerical approximation of gradient of RSS of
##' MIDAS regression 
##'
##' Calculate nummerical approximation of gradient of residual sum of
##' squares function of estimated MIDAS regression
##' 
##' @param x a \code{\link{midas_r}} object
##' @param ... further arguments to function \code{\link{hessian}} from package \code{numDeriv}
##' @return numeric vector
##' @author Vaidotas Zemlys
##' @rdname hessian 
##' @seealso midas_r
##' @import numDeriv
##' @export
hessian <- function(x,...)UseMethod("hessian")

#' @rdname hessian
#' @method hessian midas_r
#' @import numDeriv
#' @export
hessian.midas_r <- function(x,...) {
   numDeriv::hessian(x$fn0,coef(x),...)
}

##' Check whether non-linear least squares restricted MIDAS regression problem has converged
##'
##' Computes the gradient and hessian of the optimisation function of restricted MIDAS regression and checks whether the conditions of local optimum are met. Numerical estimates are used.
##' @param x \code{\link{midas_r}} object
##' @param tol a tolerance, values below the tolerance are considered zero
##' @return a list with gradient, hessian of optimisation function and convergence message
##' @rdname deriv_tests
##' @seealso midas_r
##' @export
##' @author Vaidotas Zemlys
deriv_tests<- function(x,tol=1e-6) UseMethod("deriv_tests")

#' @rdname deriv_tests
#' @method deriv_tests midas_r
#' @export
deriv_tests.midas_r <- function(x,tol=1e-6) {
    gr <- gradient(x)
    hess <- hessian(x)
    egh <- eigen(hess)$values

    first <- sum(abs(gr))<tol
    second <- !any(egh<0) & !any(abs(egh)<tol)
    list(first=first,second=second,gradient=gr,eigenval=egh)
}

##' Return the coefficients of MIDAS regression
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the regression coefficients.
##' 
##' @param x an output from \code{\link{midas_r}}
##' @return vector with coefficients of MIDAS regression
##' @author Vaidotas Zemlys
##' @export
midas_coef <- function(x) {
    x$midas.coefficients
}
##' Return the estimated parameters of the restriction function(s)
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the estimated parameters of restriction function(s)
##' 
##' @param x an output from \code{\link{midas_r}}
##' @param name name of the restriction function, the default value is the names of the restriction functions supplied to \code{\link{midas_r}}
##' @return a list if \code{length(name)>1}, a vector otherwise
##' @author Vaidotas Zemlys
##' @export
restr_param <- function(x,name=restr_names(x)) {
    if(!any(name %in% names(x$param.map))) stop("Supply valid name(s) of the restriction function")
    res <- lapply(name,function(nm)coef(x)[x$param.map[[nm]]])
    names(res) <- name
    if(length(res)==1)res <- res[[1]]
    res
}
##' Return the restricted coefficients generated by restriction function(s)
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the restricted coefficients generated by restriction function(s)
##' 
##' @param x an output from \code{\link{midas_r}}
##' @param name name(s) of the restriction function(s), the default value is the name(s) of the restriction function(s) supplied to \code{\link{midas_r}}
##' @return a list if \code{length(name)>1}, a vector otherwise
##' @author Vaidotas Zemlys
##' @export
restr_coef <- function(x,name=restr_names(x)) {
    if(!any(name %in% names(x$param.map))) stop("Supply valid name(s) of the restriction function")
    res <- lapply(name,function(nm)x$restrictions[[nm]](restr_param(x,nm)))
    names(res) <- name
    if(length(res)==1)res <- res[[1]] 
    res
}
##' Return the names of restriction function(s)
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the name(s) of restriction function(s) used in call to \code{\link{midas_r}}
##' 
##' @param x an output from \code{\link{midas_r}}
##' @return a character vector
##' @author Vaidotas Zemlys
##' @export
restr_names <- function(x) {
    names(x$restrictions)
}
##' Return the coefficients for embedlf variables
##'
##' Extracts the coefficients and returns those coefficients which name has string \code{embedlf} in it.
##' 
##' @param x an output from \code{\link{midas_u}}
##' @return a vector
##' @author Vaidotas Zemlys
##' @export
embedlf_coef <- function(x) {
    cf <- coef(x)
    cf[grep("embedlf",names(cf))]
}

##' Test restrictions on coefficients of MIDAS regression
##'
##' Perform a test whether the restriction on MIDAS regression coefficients holds.
##' 
##' @param x MIDAS regression model with restricted coefficients, estimated with \code{\link{midas_r}}
##' @param gr the gradient of the restriction function. Must return the matrix with dimensions \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the numbers of coefficients in unrestricted and restricted regressions correspondingly. Default value is \code{NULL}, which means that the numeric approximation of gradient is calculated.
##' @param ... the arguments supplied to gradient function, if \code{gr} is not \code{NULL}
##' @return a \code{htest} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Kvedaras V., Zemlys, V. \emph{Testing the functional constraints on parameters in regressions with variables of different frequency} Economics Letters 116 (2012) 250-254
##' @seealso hAhr.test
##' @examples
##' ##The parameter function
##' theta.h0 <- function(p, dk) {
##'    i <- (1:dk-1)
##'    (p[1] + p[2]*i)*exp(p[3]*i + p[4]*i^2)
##' }
##'
##' ##Generate coefficients
##' theta0 <- theta.h0(c(-0.1,0.1,-0.1,-0.001),4*12)
##'
##' ##Plot the coefficients
##' plot(theta0)
##'
##' ##Generate the predictor variable
##' set.seed(13)
##' x <- simplearma.sim(list(ar=0.6),1500*12,1,12)
##'
##' ##Simulate the response variable
##' y <- midas.sim(500,theta0,x,1)
##'
##' ##Remove unnecessary history of x
##' x <- window(x,start=start(y))
##' 
##' ##Fit restricted model
##' mr <- midas_r(y~embedlf(x,4*12,12,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,0.1,-0.1,-0.001)),dk=4*12)
##' 
##' ##The gradient function
##' grad.h0<-function(p, dk) {
##'    i <- (1:dk-1)
##'    a <- exp(p[3]*i + p[4]*i^2)
##'    cbind(a, a*i, a*i*(p[1]+p[2]*i), a*i^2*(p[1]+p[2]*i))
##' }
##'
##' ##Perform test (the expected result should be the acceptance of null)
##'
##' hAh.test(mr,gr=grad.h0,dk=4*12)
##'
##' ##Use numerical gradient instead of supplied one 
##' hAh.test(mr)
##'
##' 
##' @details  Given MIDAS regression:
##'
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+u_t}
##'
##' test the null hypothesis that the following restriction holds:
##'
##' \deqn{\theta_h=g(h,\lambda),}
##' where \eqn{h=0,...,(k+1)m}. 
##' @export
##' @import MASS
##' @import numDeriv
hAh.test <- function(x,gr=NULL,...) {

    prep <- prep_hAh(x,gr,...)
    
    unrestricted <- x$unrestricted

    se2 <- sum(residuals(unrestricted)^2)/(nrow(x$model)-prep$dk)
    A0 <- (diag(prep$dk)-prep$P%*%tcrossprod(prep$Delta.0,prep$P))/se2        
    STATISTIC <- t(prep$h.0)%*%A0%*%prep$h.0
    
    names(STATISTIC) <- "hAh"
    METHOD <- "hAh restriction test"
    PARAMETER <- prep$dk-length(coef(x))
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")
}
##' Test restrictions on coefficients of MIDAS regression using robust version of the test
##'
##' Perform a test whether the restriction on MIDAS regression coefficients holds.
##' @param x MIDAS regression model with restricted coefficients, estimated with \code{\link{midas_r}}
##' @param PHI the "meat" covariance matrix, defaults to \code{vcovHAC(x$unrestricted, sandwich=FALSE)}
##' @param gr the gradient of the restriction function. Must return the matrix with dimensions \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the numbers of coefficients in unrestricted and restricted regressions correspondingly. Default value is \code{NULL}, which means that the numeric approximation of gradient is calculated.
##' @param ... the arguments supplied to gradient function, if \code{gr} is not \code{NULL}
##' @return a \code{htest} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Kvedaras V., Zemlys, V. \emph{The statistical content and empirical testing of the MIDAS restrictions}
##' @seealso hAh.test
##' @examples
##'##The parameter function
##' theta.h0 <- function(p, dk) {
##'    i <- (1:dk-1)
##'    (p[1] + p[2]*i)*exp(p[3]*i + p[4]*i^2)
##' }
##'
##' ##Generate coefficients
##' theta0 <- theta.h0(c(-0.1,0.1,-0.1,-0.001),4*12)
##'
##' ##Plot the coefficients
##' plot(theta0)
##'
##' ##Generate the predictor variable
##' set.seed(13)
##' x <- simplearma.sim(list(ar=0.6),1500*12,1,12)
##'
##' ##Simulate the response variable
##' y <- midas.sim(500,theta0,x,1)
##'
##' ##Remove unnecessary history of x
##' x <- window(x,start=start(y))
##' 
##' ##Fit restricted model
##' mr <- midas_r(y~embedlf(x,4*12,12,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,0.1,-0.1,-0.001)),dk=4*12)
##' 
##' ##The gradient function
##' grad.h0<-function(p, dk) {
##'    i <- (1:dk-1)
##'    a <- exp(p[3]*i + p[4]*i^2)
##'    cbind(a, a*i, a*i*(p[1]+p[2]*i), a*i^2*(p[1]+p[2]*i))
##' }
##'
##' ##Perform test (the expected result should be the acceptance of null)
##'
##' hAhr.test(mr,gr=grad.h0,dk=4*12)
##'
##' ##Use numerical gradient instead of supplied one 
##' hAhr.test(mr)
##' 
##' @details  Given MIDAS regression:
##'
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+u_t}
##'
##' test the null hypothesis that the following restriction holds:
##'
##' \deqn{\theta_h=g(h,\lambda),}
##' where \eqn{h=0,...,(k+1)m}. 
##' @export
##' @import MASS
##' @import numDeriv
##' @import sandwich
hAhr.test <- function(x,PHI=vcovHAC(x$unrestricted,sandwich=FALSE),gr=NULL,...) {
    prep <- prep_hAh(x,gr,...)
    
    unrestricted <- x$unrestricted

    nyx <- nrow(x$model)
    nkx <- ncol(x$model)-1
    II <- diag(nkx)-prep$XtX %*% prep$Delta.0
    A0 <- ginv(nyx * ginv(t(prep$P)) %*% II %*% PHI %*% t(II) %*% ginv(prep$P))

    STATISTIC <- t(prep$h.0)%*%A0%*%prep$h.0
    
    names(STATISTIC) <- "hAhr"
    METHOD <- "hAh restriction test (robust version)"
    PARAMETER <- prep$dk-length(coef(x))
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")
}
##' Calculate data for \link{hAh.test} and \link{hAhr.test}
##'
##' Workhorse function for calculating necessary matrices for \link{hAh.test} and \link{hAhr.test}. Takes the same parameters as \link{hAh.test}
##' @param x \code{midas_r} object
##' @param gr a function or a list of functions for each restriction.
##' @param ... arguments supplied to gradient function(s)
##' @return a list with necessary matrices
##' @author Virmantas Kvedara, Vaidotas Zemlys
##' @seealso hAh.test, hAhr.test
prep_hAh <- function(x,gr=NULL,...) {

    unrestricted <- x$unrestricted
    if(is.null(unrestricted))stop("Unrestricted model cannot be estimated due to the lack of degrees of freedom, testing the restriction is not possible")
    
    D0 <- gradD(x,gr=gr,...)(coef(x))

    X <- x$model[,-1]
    
    XtX <- crossprod(X)

    dk <- ncol(XtX)    

    if(nrow(D0) != dk)stop("The gradient dimensions are incorrect. Number of rows does not equal number of unrestricted coefficients")
    
    P <- chol(XtX)

    cfur <- coef(unrestricted)
   
    h.0 <- P%*%(cfur-x$midas.coefficients)

    Delta.0 <- D0%*%tcrossprod(ginv(crossprod(D0,XtX)%*%D0),D0)
    
    list(P=P,XtX=XtX,dk=dk,Delta.0=Delta.0,h.0=h.0)
}
##' Gradient function of coefficient function in MIDAS regression
##'
##' Returns gradient function of coefficient function
##' 
##' @param x \code{\link{midas_r}} object
##' @param gr a function or a list of functions for each restriction.
##' @param ... arguments supplied to gradient function(s)
##' @return a function
##' @author Vaidotas Zemlys
##' @details MIDAS regression can be written in the following form:
##'
##' \deqn{y=Xf(\eta)+u}
##'
##' This function calculates
##'
##' \deqn{\frac{\partial f}{\partial \eta}}
##' 
gradD <- function(x,gr=NULL,...) {
    rf <- lapply(x$param.map,function(x)return(function(p,...)p)) 
    names(rf) <- names(x$param.map)

    restr.name <- names(x$restrictions)   
    
    if(is.null(gr)) {
        rf[names(x$restrictions)] <- x$restrictions
        all_coef <- function(p) {
             pp <- lapply(x$param.map,function(xx)p[xx])     
             res <- mapply(function(fun,param)fun(param),rf,pp,SIMPLIFY=FALSE)
             unlist(res)
         }

        gr <- function(p,...)jacobian(all_coef,p)
    }
    else {
        if(is.function(gr) & length(x$restrictions)==1) {
            grf <- list(gr)
            names(grf) <- names(x$restrictions)
        }
        else {
            if(is.null(names(gr)) | !is.list(gr)) stop("Argument gr must be a named list")
            if(any(!names(gr) %in% names(x$restrictions))) stop("The names of the list must coincide with the names of the restriction functions")
            grf <- lapply(x$param.map,function(x)return(function(p)return(matrix(1))))
            grf[names(x$restrictions)] <- gr[names(x$restrictions)]
        }
        pp0 <- lapply(x$param.map,function(xx)coef(x)[xx])            
        grmat0 <- mapply(function(fun,param,MoreArgs=...)fun(param,...),grf,pp0,SIMPLIFY=FALSE)
        colnos <- sapply(grmat0,ncol)
        rownos <- sapply(grmat0,nrow)
        np <- length(colnos)
        ccol <- cumsum(colnos)
        rrow <- cumsum(rownos)
        pinds <- cbind(c(1,rrow[-np]+1),rrow,
                       c(1,ccol[-np]+1),ccol)
        pinds <- apply(pinds,1,function(x)list(row=x[1]:x[2],col=x[3]:x[4]))
                
        gr <- function(p,...) {
            pp <- lapply(x$param.map,function(x)p[x])
            grmat <- mapply(function(fun,param,MoreArgs=...)fun(param,...),grf,pp,SIMPLIFY=FALSE)
            if(length(grmat)==1) {
                res <- grmat[[1]]
            }
            else {
                res <- matrix(0,nrow=sum(rownos),ncol=sum(colnos))
                mapply(function(m,ind){res[ind$row,ind$col] <- m},
                       grmat,pinds)
                   
            }
            res
        }
    }
    function(p)gr(p,...)
}
