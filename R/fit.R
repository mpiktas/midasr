##' Combine one low and one high frequency time series objects
##' 
##' Given low and high frequency time series, combine them into matrix 
##' suitable for MIDAS regresion estimation.
##' 
##' @param y the response variable, a \code{ts} object
##' @param x the predictor variable, a \code{ts} object
##' @param k the number of lags of predictor variable to include
##' @return a matrix
##' @author Virmantas Kvedaras, Vaidotas Zemlys 
##' @export
##' @import foreach
mmatrix.midas <- function(y, x, k=0) {
    n.x <- length(x)
    n <- length(y)
    m <- frequency(x) %/% frequency(y)    
    idx <- m*c((n.x/m-n+1):(n.x/m))
    X <- foreach(h.x=0:((k+1)*m-1), .combine='cbind') %do% {
        x[idx-h.x]
    }   
    res <- cbind(y, X)
    colnames(res) <- c("y", paste("X", rep(0:k, each=m), ".", rep(m:1, k+1), sep=""))
    res
}
##' Unrestricted MIDAS regression
##'
##' Estimate unrestricted MIDAS regression using OLS. This function is basically a wrapper for \code{lm}.
##' 
##' @param y the response variable
##' @param x the predictor variable
##' @param k the number of lags to include in MIDAS regression
##' @return \code{lm} object
##' @author Virmantas Kvedaras,Vaidotas Zemlys
##' @references Kvedaras V., Zemlys, V. \emph{Testing the functional constraints on parameters in regressions with variables of different frequency} \url{http://dx.doi.org/10.1016/j.econlet.2012.03.009}
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
##' midas.u(y,x,3)
##' 
##' @details MIDAS regression has the following form:
##' 
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+u_t}
##'
##' or alternatively
##'
##' \deqn{y_t=\sum_{h=0}^{(k+1)m}\theta_hx_{tm-h}+u_t,}
##' where \eqn{m} is the frequency of high-frequency data \eqn{x} and
##' \eqn{k} is the number of lags included in the regression.
##'
##' Given certain assumptions the coefficients can be estimated using usual OLS and they have the familiar properties associated with simple linear regression.
##' @export
midas.u <- function(y, x, k) {
    mm <- mmatrix.midas(y, x, k)
    lm(y~.-1,data=data.frame(mm))
}
##' Restricted MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares. Uses \code{optim} for optimisation. Currently only the estimates of the parameters are given, without their standard errors.
##' 
##' @param y the response variable
##' @param x the predictor variable
##' @param resfun the function which returns restricted parameters given
##' the restriction function. The parameters for the restriction function must be the supplied as numeric vector in the first argument of the function. Number of lags of the regression is calculated from the output of this function. 
##' @param start the starting values for optimisation
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
##' midas.r(y,x,theta.h0,c(-0.1,10,-10,-10),dk=4*12)
##' 
##' @details Given MIDAS regression:
##'
##' \deqn{y_t=\sum_{j=0}^k\sum_{i=0}^{m-1}\theta_{jm+i} x_{(t-j)m-i}+u_t}
##'
##' estimate the parameters of the restriction
##'
##' \deqn{\theta_h=g(h,\lambda),}
##' where \eqn{h=0,...,(k+1)m}.
##' @export
midas.r <- function(y, x, resfun, start, method="BFGS", control.optim=list(), ...) {
    crstart <- resfun(start,...)
    if(sum(is.na(crstart))>0) stop("NA coefficients for the starting values")
    m <- frequency(x) %/% frequency(y)
    k <- length(crstart) %/% m - 1

    if((k+1)*m != length(crstart)) stop("Fractional lags are not supported currently")
    
    yx <- mmatrix.midas(y, x, k)

    X <- yx[,-1]
    y <- yx[,1]
    
    fn0 <- function(p,...) {
        r <- y-X%*%resfun(p,...)
        sum(r^2)
    }
    
    opt <- optim(start,fn0,method=method,control=control.optim,...)
    list(coefficients=resfun(opt$par,...),parameters=opt$par,data=yx,opt=opt)
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
##' @references Kvedaras V., Zemlys, V. \emph{Testing the functional constraints on parameters in regressions with variables of different frequency} \url{http://dx.doi.org/10.1016/j.econlet.2012.03.009}
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
##' ###Change starting values!!!
##' mr <- midas.r(y,x,theta.h0,c(-0.1,0.1,-0.1,-0.001),dk=4*12)
##' mu <- midas.u(y,x,3)
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

    XtX <- crossprod(restricted$data[,-1])
    dk <- ncol(XtX)    

    if(nrow(D0) != dk)stop("The gradient dimensions are incorrect. Number of rows does not equal number of unrestricted coefficients")
    
    P <- chol(XtX)

    h.0 <- P%*%(coef(unrestricted)-coef(restricted))

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
