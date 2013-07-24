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
##' ##Do not run
##' #plot(theta0)
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
##' mu <- midas_u(y~fmls(x,2,12)-1, ldt, hdt)
##'
##' ##Include intercept and trend in regression
##'
##' mu.it <- midas_u(y~fmls(x,2,12)+trend, ldt, hdt)
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
##' @param x either formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are the arguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param user.gradient the default value is FALSE, which means that the numeric approximation of weight function gradient is calculated. If TRUE  it is assumed that the R function for weight function gradient has the name of the weight function appended with \code{.gradient}. This function must return the matrix with dimensions \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the numbers of coefficients in unrestricted and restricted regressions correspondingly. 
##' @param ... additional arguments supplied to optimisation function
##' @return a \code{midas_r} object which is the list with the following elements:
##' 
##' \item{coefficients}{the estimates of parameters of restrictions}
##' \item{midas.coefficientas}{the estimates of restricted coefficients of MIDAS regression}
##' \item{model}{model data}
##' \item{weights}{the MIDAS weights used in estimation.}
##' \item{unrestricted}{unrestricted regression estimated using \code{\link{midas_u}}}
##' \item{param.map}{parameter map for optimisation function}
##' \item{fn0}{optimisation function for non-linear least squares problem solved in restricted MIDAS regression}
##' \item{rhs}{the function which evaluates the right-hand side of the MIDAS regression}
##' \item{allcoef}{the function which evaluates the restricted coefficients of MIDAS regression}
##' \item{opt}{the output of optimisation procedure}
##' \item{argmap.opt}{the list containing the name of optimisation function together with arguments for optimisation function}
##' \item{start.opt}{the starting values used in optimisation}
##' \item{call}{the call to the function}
##' \item{terms}{terms object}
##' \item{gradient}{gradient of NLS objective function}
##' \item{hessian}{hessian of NLS objective function}
##' \item{Zenv}{the environment in which data is placed}
##' \item{user.gradient}{the value of supplied argument user.gradient}
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
##' mr <- midas_r(y~fmls(x,4*12-1,12,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(x=c(-0.1,10,-10,-10)))
##'
##' ##Include intercept and trend in regression
##'
##' mr.it <- midas_r(y~fmls(x,4*12-1,12,theta.h0)+trend,data.frame(y=y,trend=1:500),data.frame(x=x),start=list(x=c(-0.1,10,-10,-10)))
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
midas_r.default <- function(x, ldata=NULL, hdata=NULL, start, Ofunction="optim", user.gradient=FALSE,...) {

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
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    ##Fix this!!
    m <- match(c("x", "ldata"), names(mf), 0L)
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

    prepmd <- prepmidas_r(y,X,mt,Zenv,cl,args,start,Ofunction,user.gradient,itr$lagsTable)
    
    class(prepmd) <- "midas_r"
    midas_r.fit(prepmd)    
}

##' Restricted MIDAS regression
##'
##' Reestimate the MIDAS regression with different starting values
##' 
##' @param x \code{midas_r} object 
##' @param start the starting values
##' @param Ofunction a character string of the optimisation function to use. The default value is to use the function of previous optimisation.
##' @param ... further arguments to optimisation function. If none are supplied, the arguments of previous optimisation are used.
##' @return \code{midas_r} object
##' @method midas_r midas_r
##' @seealso midas_r
##' @author Vaidotas Zemlys
##' @export
midas_r.midas_r <- function(x,start=coef(x),Ofunction=x$argmap.opt$Ofunction,...) {
   
    oarg <- list(...)
    cl <- match.call()
    dotargnm <- names(oarg)
    
    ##Perform check whether arguments are ok and eval them
    if(length(dotargnm)>0) {
        offending <- dotargnm[!dotargnm %in% names(formals(Ofunction))]
        if(length(offending)>0)  {
            stop(paste("The function ",Ofunction," does not have the following arguments: ", paste(offending,collapse=", "),sep=""))
        }
    }
    else {
        oarg <- NULL
    }

    if(Ofunction!=x$argmap.opt$Ofunction) {
        argmap <- c(list(Ofunction=Ofunction),oarg)
    }              
    else {
         argmap <- x$argmap.opt
         argmap$Ofunction <- NULL
         argnm <- union(names(argmap),names(oarg))
         marg <- vector("list",length(argnm))
         names(marg) <- argnm
         ##New supplied arguments override the old ones
         marg[names(oarg)] <- oarg
         ##Already set arguments are left intact
         oldarg <- setdiff(names(argmap),names(oarg))
         marg[oldarg] <- argmap[oldarg]
         argmap <- c(list(Ofunction=Ofunction),marg)
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
    function.opt <- args$Ofunction
    args$Ofunction <- NULL
    if(function.opt=="optim" | function.opt=="spg") {  
        args$p <- x$start.opt
        args$fn <- x$fn0
        if(x$user.gradient) {
            args$gr <- x$gradient
        }
        opt <- do.call(function.opt,args)
        par <- opt$par
        names(par) <- names(coef(x))
    }
    if(function.opt=="nls") {
        rhs <- x$rhs
        if(x$user.gradient) {
            rhs <- function(p) {
                res <- rhs(p)
                attr(res,"gradient") <- x$gr(p)
                res
            }
        }
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

##' Prepare necessary objects for fitting of the MIDAS regression
##'
##' Prepare necessary objects for fitting of the MIDAS regression
##' @param y the response
##' @param X the model matrix
##' @param mt the terms of the formula
##' @param Zenv the environment to evaluate the formula
##' @param cl call of the function
##' @param args additional argument
##' @param start starting values
##' @param Ofunction the optimisation function
##' @param user.gradient see \link{midas_r} documentation
##' @param lagsTable the lagstable from \link{checkARstar}
##' @param unrestricted the unrestricted model
##' @author Vaidotas Zemlys
prepmidas_r <- function(y,X,mt,Zenv,cl,args,start,Ofunction,user.gradient,lagsTable,unrestricted=NULL) {
    
    ##High frequency variables can enter to formula
    ##only within fmls function
    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    term.labels <- attr(mt,"term.labels") 

    rfd <- vector("list",length(terms.lhs))

    wterm <- function(fr,type="fmls") {
         mf <- fr[-4:-5]
         mf[[1]] <- fr[[5]]
         for(j in 3:length(mf)) {
             mf[[j]] <- eval(mf[[j]],Zenv)
         }        
         mf[[3]] <- switch(type,
                           fmls = mf[[3]]+1,
                           dmls = mf[[3]]+1, 
                           mls = length(mf[[3]]))
         rf <- function(p) {
             mf[[2]] <- p
             eval(mf,Zenv)
         }
         gmf <- mf
         ##Make this customizable
         gmf[[1]] <- as.name(paste0(as.character(fr[[5]]),".gradient"))
         grf <- function(p) {
             gmf[[2]] <- p
             eval(gmf,Zenv)
         }
         return(list(weight=rf,
              name=as.character(fr[[2]]),
              gradient=grf,
              start=rep(0,mf[[3]])))
    }
    
    uterm <- function(name,k=1) {
        force(k)
        list(weight=function(p)p,
             name=name,
             gradient=function(p)diag(k),
             start=rep(0,k))
        
    }

    
    for(i in 1:length(rfd)) {
        fr <- terms.lhs[[i]]
        fun <- as.character(fr)[1] 
        rfd[[i]] <- if(fun %in% c("fmls","mls","dmls")){
            if(length(fr)>=5 && fr[[5]] != "*") {
                wterm(fr,fun)
            } else {
                lags <- eval(fr[[3]],Zenv)
                nol <- switch(fun,
                              fmls = lags+1,
                              dmls = lags+1,
                              mls = length(lags)
                              )
                uterm(term.labels[i],nol)
            }            
        }
        else {
            uterm(term.labels[i],1)
        }
    }
   
    if (attr(mt,"intercept")==1)  {
        rfd <- c(list(list(weight=function(p)p,name="(Intercept)",gradient=function(p)return(matrix(1)),start=0)),rfd)
        term.labels <- c("(Intercept)",term.labels)
    }

    rf <- lapply(rfd,with,weight)
    names(rf) <- sapply(rfd,with,name)
    

    weight_names <- setdiff(names(rf), term.labels)

    start_default <- lapply(rfd,with,start)
    names(start_default) <- names(rf)

    if(any(!weight_names%in% names(start)))stop("Starting values for weight hyperparameters must be supplied")
    
    start_default[names(start)] <- start

    restr.no <- sum(sapply(start_default[weight_names], length))
    
    np <- cumsum(sapply(start_default,length))
    
    pinds <- cbind(c(1,np[-length(np)]+1),np)
    pinds <- apply(pinds,1,function(x)list(x[1]:x[2]))
    pinds <- lapply(pinds,function(x)x[[1]])
    names(pinds) <- names(start_default)

    for(i in 1:length(start_default))names(start_default[[i]]) <- NULL
    
    starto <- unlist(start_default)
   
    if(!is.null(lagsTable)) {
      nms <- names(pinds)
      all_coef <- function(p) {
        pp <- lapply(pinds, function(x) p[x])
        cr <- c(1, -p[pinds[[grep("\"*\"", nms)]]])
        res <- mapply(function(fun, param) {
          if(!is.null(lagsTable[[param]])) {
            mltp <- rowSums(lagsTable[[param]] %*% diag(cr))
            mltp <- mltp[rowSums(lagsTable[[param]]) != 0]
          } else {
            mltp <- 1
          }
          fun(pp[[param]]) * mltp
        }, rf, 1:length(pp), SIMPLIFY = FALSE)
        unlist(res)
      }
    } else {    
      all_coef <- function(p) {              
          pp <- lapply(pinds,function(x)p[x])     
          res <- mapply(function(fun,param)fun(param),rf,pp,SIMPLIFY=FALSE)
          unlist(res)
      }
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

    if(!user.gradient) {
        gradD <- function(p)jacobian(all_coef,p)
        gr <- function(p)grad(fn0,p)
    }
    else {
        grf <- sapply(rfd,with,gradient)
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
        gr <- function(p) {
             XD <- X%*%gradD(p)
             resid <- y - X %*% all_coef(p)
             as.vector(-2*apply(as.vector(resid)*XD,2,sum))             
        }
        ##Seems to work
    }
    
    hess <- function(x)numDeriv:::hessian(fn0,x)
      
    if(is.null(unrestricted)) {
        if(ncol(X)<nrow(X)) unrestricted <- lm(y~.-1,data=data.frame(cbind(y,X),check.names=FALSE))
    }

    control <- c(list(Ofunction=Ofunction),args)
    ##Override default method of optim. Use BFGS instead of Nelder-Mead
    if(!("method"%in% names(control)) & Ofunction=="optim") {        
        control$method <- "BFGS"
    }
    
    list(coefficients=starto,
         midas.coefficients=all_coef(starto),
         model=cbind(y,X),
         weights=rf[weight_names],
         unrestricted=unrestricted,
         param.map=pinds,
         fn0=fn0,
         rhs=mdsrhs,
         allcoef=all_coef,
         opt=NULL,
         argmap.opt=control,
         start.opt=starto,
         call=cl,
         terms=mt,
         gradient=gr,
         hessian=hess,
         gradD=gradD,
         Zenv=Zenv,
         user.gradient=user.gradient,
         nobs=nrow(X))   
}

##' Restricted MIDAS regression
##'
##' Fast function for fitting MIDAS regression
##' @param y model response
##' @param X prepared matrix of high frequency variable lags
##' @param z additional low frequency variables
##' @param weight the weight function
##' @param grw the gradient of weight function
##' @param startx the starting values for weigt function
##' @param startz the starting values for additional low frequency variables
##' @param method a method passed to \link{optim}
##' @param ... additional parameters to \link{optim}
##' @return a \code{midas_r} object
##' @author Virmantas Kvedars, Vaidotas Zemlys
##' @import numDeriv
midas_r_fast <- function(y,X,z=NULL,weight,grw=NULL,startx,startz=NULL,method="BFGS",...) {

    d <- ncol(X)
    nw <- length(start)
    if(is.null(z)) {        
        all_coef <- function(p) {
            weight(p,d)
        }
        gradD <- function(p)grw(p,d)
        start <- startx
    }
    else {
        all_coef <- function(p) {
            c(weight(p[1:nw]),p[-nw:-1])
        }
        nz <- ncol(z)
        if(is.null(startz))startz <- rep(0,nz)
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
    model <- na.omit(cbind(y,X,z))
    n <- nrow(model)
    fn0 <- function(p) {
        sum((model[,1]-model[,-1]%*%all_coef(p))^2)
    }
    
    if(is.null(grw)) {
        gradD <- function(p)jacobian(all_coef,p)
        gr <- function(p)grad(fn0,p)
        gr0 <- NULL
    }
    else {
        gr <- function(p) {
             XD <- model[,-1] %*% gradD(p)
             resid <- model[,1] - model[,-1] %*% all_coef(p)
             as.vector(-2*apply(as.vector(resid)*XD,2,sum)) 
        }
        gr0 <- gr
    }
    opt <- optim(start,fn0,gr0,method=method,...)
    call <- match.call()
    fitted.values <- as.vector(model[,-1]%*%all_coef(opt$par))
    list(coefficients=opt$par,
         midas.coefficients=all_coef(opt$par),
         model=model,
         weights=weight,
         fn0=fn0,
         allcoef=all_coef,
         opt=opt,
         call=call,
         gradient=gr,
         hessian=function(x)numDeriv:::hessian(fn0,x),
         gradD=gradD,
         fitted.values=fitted.values,
         residuals=as.vector(model[,1]-fitted.values))
             
}

##' Check whether the MIDAS model is MIDAS-AR* model
##'
##' Checks whether the MIDAS model is MIDAS-AR* model and provides necessary modifications
##' @param trms terms of the model formula
##' @author Julius Vainora
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
        fun <- as.character(z[1])
        
        if(length(z) >= 4 && eval(z[[4]], env) != 1) {
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
