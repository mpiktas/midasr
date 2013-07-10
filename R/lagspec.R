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
##' @export
nealmon.gradient <- function(p,d) {
    i <- (1:d)/100
    pl <- poly(i,degree=length(p)-1,raw=TRUE)
    eplc <- exp(pl%*%p[-1])[,,drop=TRUE]
    ds <- apply(pl*eplc,2,sum)
    s <- sum(eplc)
    cbind(eplc/s,p[1]*(pl*eplc/s-eplc%*%t(ds)/s^2))
}

##' Normalized beta probability density function MIDAS weights specification
##' Calculate MIDAS weights according to normalized beta probability density function specification
##' @param p parameters for normalized beta probability density function
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
nbeta <- function(p,d) {
    eps <- .Machine$double.eps
    xi <- (1:d - 1)/(d - 1)
    xi[1] <- xi[1]+eps
    xi[d] <- xi[d]-eps
    nb <- xi^(p[2] - 1) * (1 - xi)^(p[3] - 1)
    if(sum(nb)<eps) {
        rep(0,length(nb))
    } else {
        p[1] * (nb/sum(nb) + p[4])
    }
}

##' Gradient function for normalized beta probability density function MIDAS weights specification
##' Calculate gradient function for normalized beta probability density function specification of MIDAS weights.
##' @param p parameters for normalized beta probability density function
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
nbeta.gradient <- function(p,d) {
    eps <- .Machine$double.eps
    xi <- (1:d-1)/(d-1)
    xi[1] <- xi[1]+eps
    xi[d] <- xi[d]-eps
    nb <- xi^(p[2]-1)*(1-xi)^(p[3]-1)
    snb <- sum(nb)
    if(snb>eps) {
        nba <- nb*log(xi)
        nbb <- nb*log(1-xi)
        a <- nba/snb-nb*sum(nba)/snb^2
        b <- nbb/snb-nb*sum(nbb)/snb^2
        cbind(nb/snb+p[4],p[1]*a,p[1]*b,p[1])
    }
    else {
       gres <- matrix(0,nrow=d,ncol=4)
       gres[,1] <- p[4]
       gres[,4] <- p[1]
       gres
    }    
}

##' Almon polynomial MIDAS weights specification
##'
##' Calculate Almon polynomial MIDAS weights
##' @param p parameters for Almon polynomial weights
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
almonp <- function(p,d) {
    i <- 1:d/100
    plc <- poly(i,degree=length(p)-1,raw=TRUE) %*%p[-1]+p[1]
    as.vector(plc)
}

##' Gradient function for Almon polynomial MIDAS weights
##'
##' Calculate gradient for Almon polynomial MIDAS weights specification
##' @param p vector of parameters for Almon polynomial specification
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Vaidotas Zemlys
##' @export
almonp.gradient <- function(p,d) {
    i <- 1:d/100
    plc <- poly(i,degree=length(p)-1,raw=TRUE)
    cbind(1,plc)
}

##' Step function specification for MIDAS weights
##'
##' Step function specification for MIDAS weights
##' @param p vector of parameters
##' @param d number of coefficients
##' @param a vector of increasing positive integers indicating the steps
##' @return vector of coefficients
##' @author Vaidotas Zemlys
##' @export
polystep <- function(p,d,a) {
    if(length(a)!=length(p)-1)stop("The number of steps should be number of parameters minus one")
    if(min(a)<=1 | max(a)>=d)stop("The steps are out of bounds")
    a <- c(0,a,d)
    rep(p,times=diff(a))
}
##' Gradient of step function specification for MIDAS weights
##'
##' Gradient of step function specification for MIDAS weights
##' @param p vector of parameters
##' @param d number of coefficients
##' @param a vector of increasing positive integers indicating the steps
##' @return vector of coefficients
##' @author Vaidotas Zemlys
##' @export
polystep.gradient <- function(p,d,a) {
    if(length(a)!=length(p)-1)stop("The number of steps should be number of parameters minus one")
    if(min(a)<=1 | max(a)>=d)stop("The steps are out of bounds")
    a <- c(0,a,d)
    da <- diff(a)
    r <- cbind(a[-length(a)],da)
    apply(r,1,function(x) {
        res <- rep(0,d)
        res[(x[1]+1):(x[1]+x[2])] <- 1
        res
    })
}
    
##' Normalized Gompertz probability density function MIDAS weights specification
##' Calculate MIDAS weights according to normalized Gompertz probability density function specification
##' @param p parameters for normalized Gompertz probability density function
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
gompertzp <- function(p, d) {
  i <- 1:d / d
  gm <- exp(p[3] * i - p[2] * exp(p[3] * i))
  p[1] * gm / sum(gm)
}

##' Gradient function for normalized Gompertz probability density function MIDAS weights specification
##' Calculate gradient function for normalized Gompertz probability density function specification of MIDAS weights.
##' @param p parameters for normalized Gompertz probability density function
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
gompertzp.gradient <- function(p, d) {
  i <- 1:d / d
  gm <- exp(p[3] * i - p[2] * exp(p[3] * i))
  dp2 <- -gm * exp(i * p[3])
  dp3 <- -gm * i * (p[2] * exp(i * p[3]) - 1)
  cbind(gm, p[1] * (dp2 - gm * sum(dp2) / sum(gm)),
        p[1] * (dp3 - gm * sum(dp3) / sum(gm))) / sum(gm) 
}

##' Normalized Nakagami probability density function MIDAS weights specification
##' Calculate MIDAS weights according to normalized Nakagami probability density function specification
##' @param p parameters for normalized Nakagami probability density function
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
nakagamip <- function(p, d) {
  i <- 1:d / d
  ng <- i^(2 * p[2] - 1) * exp(-p[2] / p[3] * i^2)
  p[1] * ng / sum(ng)
}

##' Gradient function for normalized Nakagami probability density function MIDAS weights specification
##' Calculate gradient function for normalized Nakagami probability density function specification of MIDAS weights.
##' @param p parameters for normalized Nakagami probability density function
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
nakagamip.gradient <- function(p, d) {
  i <- 1:d / d
  ng <- i^(2 * p[2] - 1) * exp(-p[2] / p[3] * i^2)
  dp2 <- ((2 * log(i) * p[3] - i^2) / p[3]) * ng
  dp3 <- (p[2] * i^2 / p[3]^2) * ng
  cbind(ng, (dp2 - ng * sum(dp2) / sum(ng)) * p[1],
        (dp3 - ng * sum(dp3) / sum(ng)) * p[1]) / sum(ng)
}

##' Normalized log-Cauchy probability density function MIDAS weights specification
##' Calculate MIDAS weights according to normalized log-Cauchy probability density function specification
##' @param p parameters for normalized log-Cauchy probability density function
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
lcauchyp <- function(p, d) {
  i <- 1:d / d
  lc <- 1 / (i * ((log(i) - p[2])^2 + p[3]^2))
  p[1] * lc / sum(lc)
}

##' Gradient function for normalized log-Cauchy probability density function MIDAS weights specification
##' Calculate gradient function for normalized log-Cauchy probability density function specification of MIDAS weights.
##' @param p parameters for normalized log-Cauchy probability density function
##' @param d number of coefficients
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
lcauchyp.gradient <- function(p, d) {
  i <- 1:d / d
  lc <- 1 / (i * ((log(i) - p[2])^2 + p[3]^2))
  dp2 <- 2 * lc^2 * i * (log(i) - p[2])
  dp3 <- -2 * lc^2 * p[3] * i
  cbind(lc, p[1] * (dp2 - lc * sum(dp2) / sum(lc)),
        p[1] * (dp3 - lc * sum(dp3) / sum(lc))) / sum(lc)
}
