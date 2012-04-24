##' Combine one low and one high frequency time series objects
##' 
##' Given low and high frequency time series, combine them into matrix 
##' suitable for MIDAS regresion estimation.
##' 
##' @param y the response low frequency variable, a \code{ts} object
##' @param x the predictor high frequency variable, a \code{ts} object
##' @param exo the exogenous low frequency variables, a \code{ts} object, vector or matrix. It is forcibly converted to \code{ts} object with exact properties of \code{y}. 
##' @param k the number of lags of predictor variable to include. This can either be a number indicating the highest lag, or the vector with lags which should be included. Note that the here the lag is low frequency lag. Hence effectively high frequency variable is lagged \eqn{k\times m} times, where \eqn{m} is frequency ratio, which is calculated as \code{frequency(x) %/% frequency(y)}
##' @return a matrix
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @details Note that this function does not check whether data is conformable. If \code{exo} variable does not have the same \code{start}, \code{end} and \code{frequency} as response variable \code{y}, error is produced.
##' 
##' @export
##' @import foreach
mmatrix.midas <- function(y, x, exo=NULL, k=0) {
    n.x <- length(x)
    n <- length(y)
    m <- frequency(x) %/% frequency(y)

    
    if(length(k)>1) {
        klags <- k
        k <- max(k)
    }
    else {
        klags <- 0:k
    }
        
    idx <- m*c((n.x/m-n+k+1):(n.x/m))    

        
    X <- foreach(h.x=0:((k+1)*m-1), .combine='cbind') %do% {
        x[idx-h.x]
    }   

    if(is.null(exo)) {
        exo.cn <- character()
    }
    else {
        exo <- try(ts(exo,start=start(y),end=end(y),frequency=frequency(y)))       
        if(inherits(class(exo),"try-error")) stop("Failed to convert exogenous variables to time series object")
        
        if(inherits(exo,"mts")) {
            exo <- exo[(k+1):n,]
            exo.cn <- paste("exo",1:ncol(exo),sep="")
        }
        else {
            exo <- exo[(k+1):n]
            exo.cn <- "exo1"
        }
    }

    y <- y[(k+1):n]

    
    res <- cbind(y,X,exo)
    colnames(res) <- c("y", paste("X", rep(0:k, each=m), ".", rep(m:1, k+1), sep=""),exo.cn)

    ###select only the lags needed
    nmres <- colnames(res)
    lagn <- unlist(lapply(klags,function(no) grep(paste("X",no,"[.]",sep=""),nmres,value=TRUE)))
       
    res[,c("y",lagn,exo.cn)]
}
##' Unrestricted MIDAS regression
##'
##' Estimate unrestricted MIDAS regression using OLS. This function is basically a wrapper for \code{lm}.
##' 
##' @param y the response variable, \code{\link{ts}} object.
##' @param x the predictor variable, \code{\link{ts}} object, such that \code{frequency(x) / frequency(y)} is an integer. 
##' @param exo exogenous variables, \code{\link{ts}} object, having the same properties as \code{y}. The default is NULL, meaning that no exogenous variables are used in the regression.
##' @param k the number of lags to include in MIDAS regression. See \code{\link{mmatrix.midas}} for more details.
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
##' ##Fit unrestricted model
##' mu <- midas.u(y,x,k=3)
##'
##' ##Include intercept and trend in regression
##'
##' mu.it <- midas.u(y,x,exo=cbind(1,1:500),k=3)
##' 
##' @details MIDAS regression has the following form:
##' 
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+u_t}
##'
##' or alternatively
##'
##' \deqn{y_t=\sum_{h=0}^{(k+1)m}\theta_hx_{tm-h}+u_t,}
##' where \eqn{m} is the frequency ratio and
##' \eqn{k} is the number of lags included in the regression. 
##'
##' Given certain assumptions the coefficients can be estimated using usual OLS and they have the familiar properties associated with simple linear regression.
##'
##' MIDAS regression involves times series with different frequencies. In R
##' the frequency property is set when creating time series objects
##' \code{\link{ts}}. Hence the frequency ratio \eqn{m} which figures in MIDAS regression is calculated from frequency property of time series objects supplied.
##' @export
midas.u <- function(y, x, exo = NULL, k) {
    mm <- mmatrix.midas(y, x, exo, k)
    lm(y~.-1,data=data.frame(mm))
}
##' Restricted MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares. Uses \code{optim} for optimisation. Currently only the estimates of the parameters are given, without their standard errors.
##' 
##' @param formula 
##' @param ldata 
##' @param hdata 
##' @param start the starting values for optimisation. Must be a list with named elements \code{resfun} containing vector of starting values for restriction function and \code{exo} containing vector of starting values for exogenous variables.
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
##' \item{unrestricted}{unrestricted regression estimated using \code{\link{midas.u}}
##' \item{restr.no}{the number of parameters used in restriction function}
##' }
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
##' x <- window(x,start=c(1500-500+1,1))
##' 
##' ##Fit restricted model
##' mr <- midas_r(y~mdslag(x,4,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)
##'
##' ##Include intercept and trend in regression
##'
##' mr.it <- midas_r(y~mdslag(x,4,theta0)+trend,data.frame(y=y,trend=1:500),data.frame(x=x),start=list(theta.h0=c(-0.1,10,-10,-10)),dk=4*12)
##' 
##' @details Given MIDAS regression:
##'
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+u_t}
##'
##' estimate the parameters of the restriction
##'
##' \deqn{\theta_h=g(h,\lambda),}
##' where \eqn{h=0,...,(k+1)m}.
##'
##' MIDAS regression involves times series with different frequencies. In R
##' the frequency property is set when creating time series objects
##' \code{\link{ts}}. Hence the frequency ratio \eqn{m} which figures in
##' MIDAS regression is calculated from frequency property of time series
##' objects supplied.
##'
##' The restriction function must return the restricted coefficients of
##' the MIDAS regression. If specific lag structure is needed (see
##' more in \code{\link{mmatrix.midas}}) make sure that the number of
##' coefficients returned by restriction function coincides with number
##' of columns of the lagged high frequency predictor variable matrix.
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

    nm.rf <- sapply(terms.lhs,function(fr) {
        if(length(grep("mdslag",fr))>0) {
            fr <- parse(text=fr)[[1]]
            as.character(fr[[4]])
        }
        else {
            fr
        }
    })
    
    names(rf) <- nm.rf
   
    start_default<- as.list(rep(0,length(rf)))
    names(start_default) <- names(rf)
    
    start_default[names(start)] <- start

    ##Get number of restrictions
    restr.no <- sum(sapply(start_default[setdiff(names(start_default), terms.lhs)], length))
    
    np <- cumsum(sapply(start_default,length))
    
    pinds <- cbind(c(1,np[-length(np)]+1),np)
    colnames(pinds) <- NULL
    rownames(pinds) <- NULL
    
    starto <- unlist(start_default)

    restr_coef <- function(p,...) {
        pp <- apply(pinds,1,function(x)p[x[1]:x[2]])
        res <- mapply(function(fun,param,...)fun(param,...),rf,pp,SIMPLIFY=FALSE)
        unlist(res)
    }
    
    mdslhs <- function(p,...) {       
        coefs <- restr_coef(p,...)
        X%*%coefs
    }
  
    aa <- try(mdslhs(starto))
    
    fn0 <- function(p,...) {
        r <- y - mdslhs(p,...)
        sum(r^2)
    }
    
    opt <- optim(starto,fn0,method=method,control=control.optim,...)    
        
    list(coefficients=opt$par,
         midas.coefficients=restr_coef(opt$par,...),
         model=cbind(y,X),
         opt=opt,
         restr.fun=restr_coef,
         unrestricted=lm(y~.-1,data=data.frame(cbind(y,X),check.names=FALSE)),
         restr.no=restr.no)

}

##' Test restrictions on coefficients of MIDAS regression
##'
##' Perform a test whether the restriction on MIDAS regression coefficients holds.
##' 
##' @param restricted the restricted model, estimated with \code{\link{midas.r}}
##' @param gr the gradient of the restriction function. Must return the matrix with dimensions \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the numbers of coefficients in unrestricted and restricted regressions correspondingly. Default value is \code{NULL}, which means that the numeric approximation of gradient is calculated.
##' @param ... the parameters supplied to gradient function, if \code{gr} is not \code{NULL}
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
##' x <- window(x,start=c(1500-500+1,1))
##' 
##' ##Fit restricted model
##' mr <- midas_r(y~mdslag(x,4,theta.h0),data.frame(y=y),data.frame(x=x),start=list(theta.h0=c(-0.1,0.1,-0.1,-0.001)),dk=4*12)
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
hAh.test <- function(restricted,gr=NULL,...) {

    unrestricted <-  restricted$unrestricted
      
    if(is.null(gr))gr <- function(x,...)jacobian(restricted$restr.fun,x)
    
    D0 <- gr(coef(restricted),...)

    X <- restricted$model[,-1]
    
    XtX <- crossprod(X)

    dk <- ncol(XtX)    

    if(nrow(D0) != dk)stop("The gradient dimensions are incorrect. Number of rows does not equal number of unrestricted coefficients")
    
    P <- chol(XtX)

    cfur <- coef(unrestricted)
   
    
    h.0 <- P%*%(cfur-restricted$midas.coefficients)

    Delta.0 <- D0%*%tcrossprod(ginv(crossprod(D0,XtX)%*%D0),D0)

    A0<-diag(dk)-P%*%tcrossprod(Delta.0,P)

    se2 <- sum(residuals(unrestricted)^2)/(nrow(restricted$model)-dk)
    STATISTIC <- t(h.0)%*%A0%*%h.0/se2
    
    names(STATISTIC) <- "hAh"
    METHOD <- "hAh restriction test"
    PARAMETER <- dk-restricted$restr.no
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")
    
}
