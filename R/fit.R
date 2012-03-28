##' Arange data for fitting
##' 
##' Given response and predictor variables arrange them into form for easy
##' fitting.  
##' 
##' @param y the response variable, a \code{ts} object
##' @param x the predictor variable, a \code{ts} object
##' @param k the number of lags of predictor variable to include
##' @return a matrix
##' @author Vaidotas Zemlys
##' @export
##' @import foreach
model.matrix.midas <- function(y,x,k=0) {
    n.x <- length(x)
    n <- length(y)
    m <- frequency(x) %/% frequency(y)    
    idx <- m*c((n.x/m-n+1):(n.x/m))
    X <- foreach(h.x=0:((k+1)*m-1), .combine='cbind')%do%{
        x[idx-h.x]
    }   
    res <- cbind(y,X)
    colnames(res) <- c("y",paste("X.lag",rep(0:k,each=m),".",rep(m:1,k+1),sep=""))
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
##' @author Vaidotas Zemlys
##' @export
midas.u <- function(y,x,k) {
    mm <- model.matrix.midas(y,x,k)
    lm(y~.-1,data=data.frame(mm))
}
##' Restricted MIDAS regression
##'
##' Estimate restricted MIDAS regression using non-linear least squares. Uses \code{optim} for optimisation. 
##' 
##' @param y the response variable
##' @param x the predictor variable
##' @param resfun the function which returns restricted parameters given
##' the restriction function. The parameters for the restriction function must be the supplied as numeric vector in the first argument of the function. Number of lags of the regression is calculated from the output of this function. 
##' @param gradfun the gradient of the restriction function. Must return the matrix with dimensions \code{l x p}, where \code{l} is the number of unrestricted coefficients of the regression and \code{p} is the number of parameters in the restriction function.
##' @param start the starting values for optimisation
##' @param method the method used for optimisation, see \code{optim} documentation. All methods are suported except "L-BFGS-B" and "Brent". Default method is "BFGS".
##' @param control.optim a list of control parameters for \code{optim}.
##' @param ... additional parameters supplied for \code{resfun} and \code{gradfun}
##' @return output suitable for function HAH.test
##' @author Vaidotas Zemlys
##' @export
midas.r <- function(y, x, resfun, start, method="BFGS", control.optim=list(), ...) {
    crstart <- resfun(start,...)
    if(sum(is.na(crstart))>0) stop("NA coefficients for the starting values")
    m <- frequency(x) %/% frequency(y)
    k <- length(crstart) %/% m - 1

    if((k+1)*m != length(crstart)) stop("Fractional lags are not supported currently")
    
    yx <- model.matrix.midas(y, x, k)

    X <- yx[,-1]
    y <- yx[,1]
    
    fn0 <- function(p,...) {
        r <- y-X%*%resfun(p,...)
        sum(r^2)
    }
    
    opt <- optim(fn0,start,method=method,control=control.optim,...)
    list(coefficients=resfun(opt$par),parameters=opt$par,data=yx,opt=opt)
}

##' Test restrictions on coefficients of MIDAS regression
##'
##' Perform a test whether the restriction on MIDAS regression coefficients holds.
##' 
##' @param unrestricted the unrestricted model
##' @param restricted the restricted model
##' @param gr the gradient of the restriction function. Must return the matrix with dimensions \code{l x p}, where \code{l} is the number of unrestricted coefficients of the regression and \code{p} is the number of parameters in the restriction function.
##' @param ... the parameters supplied to gradient function
##' @return a \code{htest} object
##' @author Vaidotas Zemlys
##' @export
##' @import MASS
hAh.test <- function(unrestricted,restricted,gr,...) {

    D0 <- gr(restricted$parameters,...)

    XtX <- crossprod(restricted$data[,-1])
    dk <- ncol(XtX)    
    P <- chol(XtX)

    h.0 <- P%*%(coef(unrestricted)-coef(restricted))

    Delta.0 <- D0%*%tcrossprod(ginv(crossprod(D0,XtX)%*%D0),D0)

    A0<-diag(dk)-P%*%tcrossprod(Delta.0,P)

    se2 <- sum(residuals(unrestricted)^2)/(nrow(restricted$data)-dk)
    STATISTIC <- t(h.0)%*%A0%*%h.0/se2
    
    names(statistic) <- "hAh"
    METHOD <- "hAh restriction test"
    PARAMETER <- dk-length(restricted$parameters)
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")
    
}
