##' Test restrictions on coefficients of MIDAS regression
##'
##' Perform a test whether the restriction on MIDAS regression coefficients holds.
##' 
##' @param x MIDAS regression model with restricted coefficients, estimated with \code{\link{midas_r}}
##' @return a \code{htest} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Kvedaras V., Zemlys, V. \emph{Testing the functional constraints on parameters in regressions with variables of different frequency} Economics Letters 116 (2012) 250-254
##' @seealso hAhr.test
##' @examples
##' ##The parameter function
##' theta.h0 <- function(p, dk, ...) {
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
##' mr <- midas_r(y~fmls(x,4*12-1,12,theta.h0)-1,list(y=y,x=x),start=list(x=c(-0.1,0.1,-0.1,-0.001)))
##'
##' ##Perform test (the expected result should be the acceptance of null)
##'
##' hAh.test(mr)
##' 
##' ##Fit using gradient function
##'
##' ##The gradient function
##' theta.h0.gradient<-function(p, dk,...) {
##'    i <- (1:dk-1)
##'    a <- exp(p[3]*i + p[4]*i^2)
##'    cbind(a, a*i, a*i*(p[1]+p[2]*i), a*i^2*(p[1]+p[2]*i))
##' }
##'
##' mr <- midas_r(y~fmls(x,4*12-1,12,theta.h0)-1,list(y=y,x=x),start=list(x=c(-0.1,0.1,-0.1,-0.001)),user.gradient=TRUE)
##'
##' ##The test will use an user supplied gradient of weight function. See the
##' ##help of midas_r on how to supply the gradient.
##' 
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
##' @import numDeriv
hAh.test <- function(x) {

    prep <- prep_hAh(x)
    
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
##' @return a \code{htest} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Kvedaras V., Zemlys, V. \emph{The statistical content and empirical testing of the MIDAS restrictions}
##' @seealso hAh.test
##' @examples
##'##The parameter function
##' theta.h0 <- function(p, dk, ...) {
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
##' mr <- midas_r(y~fmls(x,4*12-1,12,theta.h0)-1,list(y=y,x=x),start=list(x=c(-0.1,0.1,-0.1,-0.001)))
##' 
##' ##The gradient function
##' theta.h0.gradient <-function(p, dk,...) {
##'    i <- (1:dk-1)
##'    a <- exp(p[3]*i + p[4]*i^2)
##'    cbind(a, a*i, a*i*(p[1]+p[2]*i), a*i^2*(p[1]+p[2]*i))
##' }
##'
##' ##Perform test (the expected result should be the acceptance of null)
##'
##' hAhr.test(mr)
##' 
##' mr <- midas_r(y~fmls(x,4*12-1,12,theta.h0)-1,list(y=y,x=x),start=list(x=c(-0.1,0.1,-0.1,-0.001)),user.gradient=TRUE)
##'
##' ##Use exact gradient. Note the 
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
##' @importFrom MASS ginv
##' @import sandwich
hAhr.test <- function(x,PHI=vcovHAC(x$unrestricted,sandwich=FALSE)) {
    prep <- prep_hAh(x)
    
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
##' @return a list with necessary matrices
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @seealso hAh.test, hAhr.test
prep_hAh <- function(x) {

    unrestricted <- x$unrestricted
    if(is.null(unrestricted))stop("Unrestricted model cannot be estimated due to the lack of degrees of freedom, testing the restriction is not possible")
    
    D0 <- x$gradD(coef(x))

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

##' Andreou, Ghysels, Kourtellos LM test
##'
##' Perform the test whether hyperparameters of normalized exponential Almon lag weights are zero
##' 
##' @param x MIDAS regression object of class \code{\link{midas_r}}
##' @return a \code{htest} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Andreou E., Ghysels E., Kourtellos A. \emph{Regression models with mixed sampling frequencies} Journal of Econometrics 158 (2010) 246-261 
##' @export
##' @examples
##' ##' ##Load data
##' data("USunempr")
##' data("USrealgdp")
##'
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' t <- 1:length(y)
##'
##' mr <- midas_r(y~t+fmls(x,11,12,nealmon),start=list(x=c(0,0,0)))
##'
##' agk.test(mr)
##'
agk.test <- function(x) {
    tl <- attributes(x$terms)$term.labels
    wn <- grep("nealmon",tl,value=TRUE)   
    if(length(wn)==0)stop("This test can be only used for regressions with normalized Exponential Almon lag weights")    
    X <- x$model[,-1]
    y <- x$model[,1]
    if(attr(x$terms,"intercept")==1) tl <- c("(Intercept)",tl)
    Xa <- lapply(tl,function(nm) {
        if(nm %in% wn) {       
            apply(X[,grep(wn,colnames(X),fixed=TRUE)],1,mean)
        }
        else {
            X[,nm,drop=FALSE]
        }
    })
    Xa <- do.call("cbind",Xa)
   
    ustar <- residuals(lm(y~Xa-1))
    u <- residuals(x)
    w <- weight_param(x)
    if(is.list(w)) {
        r <- sum(sapply(w,length))
    }
    else {
        r <- length(w)
    }
    S.LS <- sum(ustar^2)
    S.M <- sum(u^2)
      
    STATISTIC <- (S.LS-S.M)/S.LS
    names(STATISTIC) <- "agk"
    METHOD <- "Andreou, Ghysels, Kourtellos LM test"
    PARAMETER <- r
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")
}
