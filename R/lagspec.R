##' Normalized Exponential Almon lag MIDAS coefficients
##'
##' Calculate normalized exponential Almon lag coefficients given the parameters and required number of coefficients.
##' @param p parameters for Almon lag
##' @param d number of the coefficients
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @examples
##'
##' ##Load data
##' data("USunempr")
##' data("USrealgdp")
##'
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' t <- 1:length(y)
##'
##' midas_r(y~t+fmls(x,11,12,nealmon),start=list(x=c(0,0,0)))
##' 
##' @details Given unrestricted MIDAS regression
##'
##' \deqn{y_t=\sum_{h=0}^d\theta_{h}x_{tm-h}+\mathbf{z_t}\beta+u_t}
##' 
##' normalized exponential Almon lag restricts the coefficients \eqn{theta_h} in the following way:
##'
##' \deqn{\theta_{h}=\delta\frac{\exp(\lambda_1(h+1)+\dots+\lambda_r(h+1)^r)}{\sum_{s=0}^d\exp(\lambda_1(s+1)+\dots+\lambda_r(h+1)^r)}}
##'
##' The parameter \eqn{\delta} should be the first element in vector \code{p}. The degree of the polynomial is then decided by the number of the remaining parameters.
##' @export
nealmon <- function(p,d) {
  i <- (1:d)/100
  plc <- poly(i,degree=length(p)-1,raw=TRUE) %*% p[-1]
  as.vector(p[1] * exp(plc)/sum(exp(plc)))
}

##' Gradient function for normalized exponential Almon lag weights
##'
##' Gradient function for normalized exponential Almon lag weights
##' @param p hyperparameters for Almon lag
##' @param d number of coefficients
##' @return the gradient matrix
##' @author Vaidotas Zemlys
nealmon.gradient <- function(p,d) {
    i <- (1:d)/100
    pl <- poly(i,degree=length(p)-1,raw=TRUE)
    eplc <- exp(pl%*%p[-1])[,,drop=TRUE]
    ds <- apply(pl*eplc,2,sum)
    s <- sum(eplc)
    cbind(eplc/s,p[1]*(pl*eplc/s-eplc%*%t(ds)/s^2))
}
