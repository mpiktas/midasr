##' Combine one low and one high frequency time series objects
##' 
##' Given low and high frequency time series, combine them into matrix 
##' suitable for MIDAS regresion estimation.
##' 
##' @param y the response low frequency variable, a \code{ts} object
##' @param x the predictor high frequency variable, a \code{ts} object
##' @param exo the exogenous low frequency variables, a \code{ts} object.
##' @param k the number of lags of predictor variable to include
##' @return a matrix
##' @author Virmantas Kvedaras, Vaidotas Zemlys 
##' @export
##' @import foreach
mmatrix.midas <- function(y, x, exo=NULL, k=0) {
    n.x <- length(x)
    n <- length(y)
    m <- frequency(x) %/% frequency(y)

    idx <- m*c((n.x/m-n+k+1):(n.x/m))    
    y <- y[(k+1):n]
        
    X <- foreach(h.x=0:((k+1)*m-1), .combine='cbind') %do% {
        x[idx-h.x]
    }   

    if(is.null(exo)) {
        exo.cn <- character()
    }
    else {
        if(inherits(exo,"mts")) {
            exo <- exo[(k+1):n,]
            exo.cn <- paste("exo",1:ncol(exo),sep="")
        }
        else {
            exo <- exo[(k+1):n]
            exo.cn <- "exo1"
        }
    }
    res <- cbind(y,X,exo)
    colnames(res) <- c("y", paste("X", rep(0:k, each=m), ".", rep(m:1, k+1), sep=""),exo.cn)
    res
}
##' Unrestricted MIDAS regression
##'
##' Estimate unrestricted MIDAS regression using OLS. This function is basically a wrapper for \code{lm}.
##' 
##' @param y the response variable, \code{\link{ts}} object.
##' @param x the predictor variable, \code{\link{ts}} object, such that \code{frequency(x) / frequency(y)} is an integer. 
##' @param exo exogenous variables, \code{\link{ts}} object, having the same properties as \code{y}. The default is NULL, meaning that no exogenous variables are used in the regression.
##' @param k the number of lags to include in MIDAS regression.
##' @return \code{\link{lm}} object
##' @author Virmantas Kvedaras,Vaidotas Zemlys
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
##' midas.u(y,x,k=3)
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
##' @param y the response variable, \code{\link{ts}} object.
##' @param x the predictor variable, \code{\link{ts}} object, such that \code{frequency(x) / frequency(y)} is an integer.
##' @param exo exogenous variables, \code{\link{ts}} object, having the same properties as \code{y}. The default is NULL, meaning that no exogenous variables are used in the regression.
##' @param resfun the function which returns restricted parameters given
##' the restriction function. The parameters for the restriction function must be the supplied as numeric vector in the first argument of the function. Number of lags of the regression is calculated from the output of this function. 
##' @param start the starting values for optimisation. Must be a list with named elements \code{resfun} containing vector of starting values for restriction function and \code{exo} containing vector of starting values for exogenous variables.
##' @param method the method used for optimisation, see \code{\link{optim}} documentation. All methods are suported except "L-BFGS-B" and "Brent". Default method is "BFGS".
##' @param control.optim a list of control parameters for \code{\link{optim}}.
##' @param ... additional parameters supplied for \code{resfun} and \code{gradfun}
##' @return output suitable for function \code{\link{hAh.test}}
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
##' ##Fit restricted model
##' midas.r(y,x,resfun=theta.h0,start=list(resfun=c(-0.1,10,-10,-10)),dk=4*12)
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
##'MIDAS regression involves times series with different frequencies. In R
##' the frequency property is set when creating time series objects
##' \code{\link{ts}}. Hence the frequency ratio \eqn{m} which figures in
##' MIDAS regression is calculated from frequency property of time series
##' objects supplied. 
##' @export
midas.r <- function(y, x, exo=NULL, resfun, start=list(resfun=NULL,exo=NULL), method="BFGS", control.optim=list(), ...) {
    crstart <- resfun(start$resfun,...)
    if(sum(is.na(crstart))>0) stop("NA coefficients for the starting values")
    m <- frequency(x) %/% frequency(y)
    k <- length(crstart) %/% m - 1

    if((k+1)*m != length(crstart)) stop("Fractional lags are not supported currently")
    
    yx <- mmatrix.midas(y, x, exo, k)

    if(is.null(exo)) {
        X <- yx[,-1]
        y <- yx[,1]
        fn0 <- function(p,...) {
            r <- y-X%*%resfun(p,...)
              sum(r^2)
        }
        starto <- start$resfun
    }
    else {
        
        exonm <- grep("exo",colnames(yx),value=TRUE)
        
        X <- yx[,setdiff(colnames(yx),c("y",exonm))]
        y <- yx[,1]
        exom <- yx[,exonm]
        
        np <- cumsum(sapply(start,length))
        
        fn0 <- function(p,...) {
            r <- y-X%*%resfun(p[1:np[1]],...) - exom%*%p[(np[1]+1):np[2]]
            sum(r^2)
        }
        starto <- unlist(start[c("resfun","exo")])
    }
    opt <- optim(starto,fn0,method=method,control=control.optim,...)
    if(is.null(exo)) {
        resfun.param <- opt$par
        exo.coef <- NULL
    }
    else {
        resfun.param <- opt$par[1:np[1]]
        exo.coef <- opt$par[(np[1]+1):np[2]]
    }
    list(coefficients=resfun(resfun.param,...),parameters=resfun.param,data=yx,opt=opt,exo.coef=exo.coef)
}

##' Test restrictions on coefficients of MIDAS regression
##'
##' Perform a test whether the restriction on MIDAS regression coefficients holds.
##' 
##' @param unrestricted the unrestricted model, estimated with \code{\link{midas.u}}
##' @param restricted the restricted model, estimated with \code{\link{midas.r}}
##' @param gr the gradient of the restriction function. Must return the matrix with dimensions \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the numbers of coefficients in unrestricted and restricted regressions correspondingly.
##' @param ... the parameters supplied to gradient function
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
##' ##Fit restricted model
##' mr <- midas.r(y,x,resfun=theta.h0,start=list(resfun=c(-0.1,0.1,-0.1,-0.001)),dk=4*12)
##' mu <- midas.u(y,x,k=3)
##'
##' ##The gradient function
##' grad.h0<-function(p, dk) {
##'    alpha <- p[1]
##'    beta <- p[2]
##'    lambda <- c(p[3],p[4])
##'    index <- c(1:dk)
##'    i <- (index-1)
##'    pol <- poly(i,2,raw=TRUE) %*%lambda
##'    a <- (alpha+beta*i)*exp(pol)
##'    cbind(a, a*i, a*i*(alpha+beta*i), a*i^2*(alpha+beta*i))
##' }
##'
##' ##Perform test (the expected result should be the acceptance of null)
##'
##' hAh.test(mu,mr,grad.h0,dk=4*12)
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
hAh.test <- function(unrestricted,restricted,gr,...) {

    D0 <- gr(restricted$parameters,...)

    exonm <- grep("exo",colnames(restricted$data),value=TRUE)
    if(length(exonm)>0) {
        X <- restricted$data[,setdiff(colnames(restricted$data),c("y",exonm))]
    }
    else {
        X <- restricted$data[,-1]
    }
    
    XtX <- crossprod(X)

    dk <- ncol(XtX)    

    if(nrow(D0) != dk)stop("The gradient dimensions are incorrect. Number of rows does not equal number of unrestricted coefficients")
    
    P <- chol(XtX)

    cfur <- coef(unrestricted)
    cfur <- cfur[setdiff(names(cfur),exonm)]
    
    h.0 <- P%*%(cfur-coef(restricted))

    Delta.0 <- D0%*%tcrossprod(ginv(crossprod(D0,XtX)%*%D0),D0)

    A0<-diag(dk)-P%*%tcrossprod(Delta.0,P)

    se2 <- sum(residuals(unrestricted)^2)/(nrow(restricted$data)-dk)
    STATISTIC <- t(h.0)%*%A0%*%h.0/se2
    
    names(STATISTIC) <- "hAh"
    METHOD <- "hAh restriction test"
    PARAMETER <- dk-length(restricted$parameters)
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")
    
}
