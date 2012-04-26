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
##' mu <- midas_u(y~mdslag(x,3)-1, ldt, hdt)
##'
##' ##Include intercept and trend in regression
##'
##' mu.it <- midas_u(y~mdslag(x,3)+trend, ldt, hdt)
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
midas_u <- function(formula, ldata, hdata,...) {
    Zenv <- new.env(parent=environment(formula))
      
    data <- mds(ldata,hdata)

    ee <- as.environment(c(as.list(data$lowfreq),as.list(data$highfreq)))
    parent.env(ee) <- parent.frame()
    assign("ee",ee,Zenv)
    
    mf <- match.call(expand.dots = TRUE)
    mf <- mf[-4]
    mf[[1L]] <- as.name("lm")
    names(mf)[3] <- "data"
    mf[[3L]] <- as.name("ee")   
    eval(mf,Zenv)
}
##' Restricted MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares. Uses \code{optim} for optimisation. Currently only the estimates of the parameters are given, without their standard errors.
##' 
##' @param formula formula for restricted MIDAS regression
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param method the method used for optimisation, see \code{\link{optim}} documentation. All methods are suported except "L-BFGS-B" and "Brent". Default method is "BFGS".
##' @param control.optim a list of control parameters for \code{\link{optim}}.
##' @param ... additional arguments supplied for \code{resfun}
##' @return a list with the following elements:
##' 
##' \item{coefficients}{the estimates of parameters of restrictions}
##' \item{midas.coefficientas}{the estimates of restricted coefficients of MIDAS regression}
##' \item{model}{model data}
##' \item{opt}{the output of call to \code{\link{optim}}}
##' \item{restr.fun}{the restriction function used in estimation.}
##' \item{unrestricted}{unrestricted regression estimated using \code{\link{midas_u}}}
##' \item{restr.no}{the number of parameters used in restriction function}
##' 
##' 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
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
##' mr <- midas_r(y~mdslag(x,3,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)
##'
##' ##Include intercept and trend in regression
##'
##' mr.it <- midas_r(y~mdslag(x,3,theta.h0)+trend,data.frame(y=y,trend=1:500),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)
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
midas_r <- function(formula, ldata, hdata, start, method="BFGS", control.optim=list(), ...) {

    Zenv <- new.env(parent=environment(formula))
      
    data <- mds(ldata,hdata)

    ee <- as.environment(c(as.list(data$lowfreq),as.list(data$highfreq)))
    parent.env(ee) <- parent.frame()
    assign("ee",ee,Zenv)
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "ldata", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    names(mf)[3] <- "data"
    mf[[3L]] <- as.name("ee")   
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
  
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)    
    
    cl <- match.call()

    ##High frequency variables can enter to formula
    ##only within mdslag function
    terms.lhs <- attr(mt,"term.labels") 
    if(attr(mt,"intercept")==1) terms.lhs <- c("(Intercept)",terms.lhs)
    
    rf <- lapply(terms.lhs,function(fr) {     
        if(length(grep("mdslag",fr))>0) {
            fr <- parse(text=fr)[[1]]
            ff <- eval(fr[[4]],parent.frame())
            rf.arg <- formals(ff)
            class(rf.arg) <- "list"
            rf.argnm <-  intersect(names(rf.arg)[-1],names(cl))
            
            for(i in rf.argnm) {
                rf.arg[[i]] <- eval(cl[[i]],parent.frame())
            }
            
            rf.name <- as.character(fr[[4]])
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
        if(length(grep("mdslag",fr))>0) {
            fr <- parse(text=fr)[[1]]
            as.character(fr[[4]])
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
    
    starto <- unlist(start_default)

    all_coef <- function(p) {              
        pp <- lapply(pinds,function(x)p[x])     
        res <- mapply(function(fun,param)fun(param),rf,pp,SIMPLIFY=FALSE)
        unlist(res)
    }   
    
    mdslhs <- function(p) {       
        coefs <- all_coef(p)
        X%*%coefs
    }
  
    aa <- try(mdslhs(starto))
    
    fn0 <- function(p,...) {
        r <- y - mdslhs(p)
        sum(r^2)
    }
    
    opt <- optim(starto,fn0,method=method,control=control.optim,...)    
        
    res <- list(coefficients=opt$par,
         midas.coefficients=all_coef(opt$par),
         model=cbind(y,X),
         opt=opt,
         restrictions=rf[restr.name],
         unrestricted=lm(y~.-1,data=data.frame(cbind(y,X),check.names=FALSE)),
         param.map=pinds)
    class(res) <- "midas_r"
    res
}
##' Return the coefficients of MIDAS regression
##'
##' A helper function for working with output of \code{\link{midas_r}}. Returns the regression coefficients.
##' 
##' @param x an output from \code\{\link{midas_r}}
##' @return vector with coefficients of MIDAS regression
##' @author Vaidotas Zemlys
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
restr_names <- function(x) {
    names(x$restrictions)
}
##' Return the coefficients for mdslag variables
##'
##' Extracts the coefficients and returns those coefficients which name has string \code{mdslag} in it.
##' 
##' @param x an output from \code{\link{midas_u}}
##' @return a vector
##' @author Vaidotas Zemlys
mdslag_coef <- function(x) {
    cf <- coef(x)
    cf[grep("mdslag",names(cf))]
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
##' mr <- midas_r(y~mdslag(x,3,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,0.1,-0.1,-0.001)),dk=4*12)
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
##' hAh.test(mr,grad.h0,dk=4*12)
##'
##' ##Use numerical gradient instead of supplied one 
##' hAh.test(mr)
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

    unrestricted <-  x$unrestricted

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
  
##  restr.no <- sum(sapply(x$param.map[names(x$restrictions)],length))    
    
    D0 <- gr(coef(x),...)

    X <- x$model[,-1]
    
    XtX <- crossprod(X)

    dk <- ncol(XtX)    

    if(nrow(D0) != dk)stop("The gradient dimensions are incorrect. Number of rows does not equal number of unrestricted coefficients")
    
    P <- chol(XtX)

    cfur <- coef(unrestricted)
   
    
    h.0 <- P%*%(cfur-x$midas.coefficients)

    Delta.0 <- D0%*%tcrossprod(ginv(crossprod(D0,XtX)%*%D0),D0)

    A0<-diag(dk)-P%*%tcrossprod(Delta.0,P)

    se2 <- sum(residuals(unrestricted)^2)/(nrow(x$model)-dk)
    STATISTIC <- t(h.0)%*%A0%*%h.0/se2
    
    names(STATISTIC) <- "hAh"
    METHOD <- "hAh restriction test"
    PARAMETER <- dk-length(coef(x))
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD,aa=list(D0=D0,P=P,XtX=XtX,h.0=h.0,Delta.0=Delta.0,A0=A0,se2=se2)), 
        class = "htest")
    
}
