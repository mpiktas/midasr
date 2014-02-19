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
##' theta.h0 <- function(p, dk, ...) {
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
##' mu <- midas_u(y~fmls(x,2,12)-1, list(ldt, hdt))
##'
##' ##Include intercept and trend in regression
##'
##' mu.it <- midas_u(y~fmls(x,2,12)+trend, list(ldt, hdt))
##'
##' ##Pass data as partialy named list
##'
##' mu.it <- midas_u(y~fmls(x,2,12)+trend, list(ldt, x=hdt$x))
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
   
    out <- eval(mf,Zenv)
    out$Zenv <- Zenv
    class(out) <- c("midas_u",class(out))
    out
}

##' Restricted MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares.
##'
##' @param x either formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param data a named list containing data with mixed frequencies
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
##' @references Clements, M. and Galvao, A., \emph{Macroeconomic Forecasting With Mixed-Frequency Data: Forecasting Output Growth in the United States}, Journal of Business and Economic Statistics, Vol.26 (No.4), (2008) 546-554
##' @rdname midas_r
##' @seealso midas_r.midas_r
##' @examples
##' ##The parameter function
##' theta.h0 <- function(p, dk, ...) {
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
##' mr <- midas_r(y~fmls(x,4*12-1,12,theta.h0)-1,
##'               list(y=y,x=x),
##'               start=list(x=c(-0.1,10,-10,-10)))
##'
##' ##Include intercept and trend in regression
##' mr.it <- midas_r(y~fmls(x,4*12-1,12,theta.h0)+trend,
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
##' mr.ar <- midas_r(y.ar ~ trend + mls(y.ar, 1, 1) +
##'                  fmls(xx, 11, 12, nealmon),
##'                  start = list(xx = rep(0, 3)))
##' 
##' ##First order MIDAS-AR* restricted model 
##' mr.arstar <-  midas_r(y.ar ~ trend + mls(y.ar, 1, 1, "*")
##'                      + fmls(xx, 11, 12, nealmon),
##'                      start = list(xx = rep(0, 3)))
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
##' MIDAS-AR* (a model with a common factor, see (Clements and Galvao, 2008)) can be estimated by specifying additional argument, see an example.
##'
##' The restriction function must return the restricted coefficients of
##' the MIDAS regression.
##' 
##' @export
midas_r <- function(x,...)UseMethod("midas_r")


#' @rdname midas_r
#' @method midas_r default
#' @export
midas_r.default <- function(x, data, start, Ofunction="optim", user.gradient=FALSE,...) {

    Zenv <- new.env(parent=environment(x))

    if(missing(data)) {
        ee <- NULL
    }
    else {
        ee <- data_to_env(data)
        parent.env(ee) <- parent.frame()
    }
    
    assign("ee",ee,Zenv)
    x <- as.formula(x)
    cl <- match.call()    
    mf <- match.call(expand.dots = FALSE)
    mf$x <- x
    m <- match(c("x", "data"), names(mf), 0L)
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
        args$par <- x$start.opt
        args$fn <- x$fn0
        if(x$user.gradient) {
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
        args$par <- x$start.opt
        args$fn <- x$fn0
        if(x$user.gradient) {
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
        if(x$user.gradient) {
            orhs <- rhs
            rhs <- function(p) {
                res <- orhs(p)
                attr(res,"gradient") <- x$model[,-1]%*%x$gradD(p)
                res
            }
        }
        y <- x$model[,1]
        args$formula <- formula(y~rhs(p))
        args$start <- list(p=x$start.opt)
        opt <- try(do.call("nls",args),silent=TRUE)
        if(inherits(opt,"try-error")) {
            stop("The optimisation algorithm of MIDAS regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        par <- coef(opt)
        names(par) <- names(coef(x))
        x$convergence <- opt$convInfo$stopCode
    }
    x$opt <- opt
    x$coefficients <- par
    names(par) <- NULL
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
##' @param guess_start if TRUE, get the initial values for non-MIDAS terms via OLS, if FALSE, initialize them with zero.
##' @author Vaidotas Zemlys
prepmidas_r <- function(y,X,mt,Zenv,cl,args,start,Ofunction,user.gradient,lagsTable,unrestricted=NULL,guess_start=TRUE) {
    
    ##High frequency variables can enter to formula
    ##only within fmls function
    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    term.labels <- attr(mt,"term.labels") 

    rfd <- vector("list",length(terms.lhs))

    wterm <- function(fr,type="fmls") {
         mf <- fr[-5]
         mf[[1]] <- fr[[5]]
         noarg <- length(formals(eval(fr[[5]],Zenv)))
         if(noarg<2)stop("The weight function must have at least two arguments")
         mf <- mf[1:min(length(mf),noarg+1)]
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
                     start=rep(0,mf[[3]]),
                     wlabel=as.character(fr[[2]])))
    }
    
    uterm <- function(name,k=1) {
        force(k)
        list(weight=function(p)p,
             name=name,
             gradient=function(p)diag(k),
             start=rep(0,k),
             wlabel="")
        
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
                nm <- as.character(fr[[2]])
                uterm(nm,nol)
            }            
        }
        else {
            uterm(term.labels[i],1)
        }
    }
   
    if (attr(mt,"intercept")==1)  {
        rfd <- c(list(list(weight=function(p)p,name="(Intercept)",gradient=function(p)return(matrix(1)),start=0,wlabel="")),rfd)
        term.labels <- c("(Intercept)",term.labels)
    }
    
    rf <- lapply(rfd,"[[","weight")
    names(rf) <- sapply(rfd,"[[","name")
    

    weight_names <- sapply(rfd,"[[","wlabel")
    weight_inds <- which(weight_names!="")
    weight_names <- weight_names[weight_names!=""]
    
    
    start_default <- lapply(rfd,"[[","start")
    names(start_default) <- names(rf)

    ##If there are no weight functions, we have unrestricted MIDAS model.
    if(length(weight_names)==0)Ofunction <- "lm"
    else {
        if(any(!weight_names%in% names(start)))stop("Starting values for weight hyperparameters must be supplied")
    }
    
    start_default[names(start)] <- start

    
    #restr.no <- sum(sapply(start_default[weight_names], length))
    
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

    if(length(weight_names)>0 && guess_start) {
        npx <- cumsum(sapply(mapply(function(fun,st)fun(st),rf,start_default,SIMPLIFY=FALSE),length))
        xinds <- build_indices(npx,names(start_default))
        wi <- rep(FALSE,length(rf))
        wi[weight_inds] <- TRUE
        Xstart <- mapply(function(fun,st,inds,iswhgt) {        
            if(iswhgt) {
                X[,inds] %*% fun(st)
            }
            else X[,inds]
        }, rf,start_default,xinds,wi,SIMPLIFY=FALSE)

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
    }
    
    starto <- unlist(start_default)
   
    if(!is.null(lagsTable)) {
        yname <- all.vars(mt[[2]])
        nms <- names(pinds)
      all_coef <- function(p) {
        pp <- lapply(pinds, function(x) p[x])
        cr <- c(1, -p[pinds[[yname]]])
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
  
    #aa <- try(mdsrhs(starto))
    
    fn0 <- function(p,...) {
        r <- y - mdsrhs(p)
        sum(r^2)
    }

    if(!user.gradient) {
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
          gradD <- function(p) {
              dind <- which(names(pinds)==yname)
            cr <- c(1, -p[pinds[[dind]]])
            pp <- lapply(pinds, function(x) p[x])
            grmat <- mapply(function(fun, param) fun(param), grf, pp, SIMPLIFY = FALSE)
            if(length(grmat) == 1) {
              res <- grmat[[1]]
            } else {
              res <- matrix(0, nrow = sum(rownos), ncol = sum(colnos))
              mlv <- matrix(0, nrow = sum(rownos), ncol = length(pinds[[dind]]))
              for(j in 1:length(grmat)) {
                ind <- pindm[[j]]
                mltp <- 1
                if(!is.null(lagsTable[[j]])) {
                  mltp <- rowSums(lagsTable[[j]] %*% diag(cr))
                  mltp <- mltp[rowSums(lagsTable[[j]]) != 0]
                } else if(j == dind) {
                  hfv <- !sapply(lagsTable, is.null)
                  for(k in which(hfv)) {
                    for(cl in 1:ncol(mlv) + 1) {
                      ltb <- lagsTable[[k]]
                      mlv[pindm[[k]]$row, cl - 1] <- mlv[ind$row, cl - 1] + rf[[k]](pp[[k]]) * ltb[rowSums(ltb) != 0, cl]
                    }
                  }
                }
                res[ind$row,ind$col] <- grmat[[j]] * mltp
                if(j == dind) 
                  res[, ind$col] <- res[, ind$col] - mlv %*% grmat[[dind]]
              }
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
         weights=rf[weight_inds],
         unrestricted=unrestricted,
         param.map=pinds,
         fn0=fn0,
         rhs=mdsrhs,
         allcoef=all_coef,
         opt=NULL,
         argmap.opt=control,
         start.opt=starto,
         start.list=start,
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
##' Function for fitting MIDAS regression without the formula interface
##' @param y model response
##' @param X prepared matrix of high frequency variable lags
##' @param z additional low frequency variables
##' @param weight the weight function
##' @param grw the gradient of weight function
##' @param startx the starting values for weight function
##' @param startz the starting values for additional low frequency variables
##' @param method a method passed to \link{optim}
##' @param ... additional parameters to \link{optim}
##' @return an object similar to \code{midas_r} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @import numDeriv
##' @import optimx
##' @export
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
         midas.coefficients=all_coef(par),
         model=model,
         weights=weight,
         fn0=fn0,
         allcoef=all_coef,
         opt=opt,
         call=call,
         gradient=gr,
         hessian=function(x)numDeriv::hessian(fn0,x),
         gradD=gradD,
         fitted.values=fitted.values,
         residuals=as.vector(y-fitted.values))
             
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
