##' Normalized Exponential Almon lag MIDAS coefficients
##'
##' Calculate normalized exponential Almon lag coefficients given the parameters and required number of coefficients.
##' @param p parameters for Almon lag
##' @param d number of the coefficients
##' @param m the frequency, currently ignored.
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
nealmon <- function(p,d,m) {
  i <- 1:d
  plc <- poly(i,degree=length(p)-1,raw=TRUE) %*% p[-1]
  as.vector(p[1] * exp(plc)/sum(exp(plc)))
}


##Improve the following documentation                                      

##' Produces weights for aggregates based MIDAS regression
##'
##' Suppose a weight function \eqn{w(\beta,\theta)} satisfies the following equation:
##' \deqn{w(\beta,\theta)=\beta g(\theta)}
##' 
##' The following combinations are defined, corresponding to structure types \code{A}, \code{B} and \code{C} respectively:
##' \deqn{(w(\beta_1,\theta_1),...,w(\beta_k,\theta_k))}
##' \deqn{(w(\beta_1,\theta),...,w(\beta_k,\theta))}
##' \deqn{\beta(w(1,\theta_1),...,w(1,\theta_k))}
##'
##' The starting values \eqn{p} should be supplied then as follows:
##' \deqn{(\beta_1,\theta_1,...,\beta_k,\theta_k)}
##' \deqn{(\beta_1,...,\beta_k,\theta)}
##' \deqn{(\beta,\theta_1,...,\theta_k)}
##'
##' 
##' @title Weights for aggregates based MIDAS regressions
##' @param p parameters for weight functions, see details.
##' @param d number of high frequency lags
##' @param m the frequency
##' @param weight the weight function
##' @param type type of structure, a string, one of A, B or C.
##' @return a vector of weights
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
amweights <- function(p,d,m,weight=nealmon,type=c("A","B","C")) {
    type <- match.arg(type)
    hf <- d%/%m
    
    if(d%%m!=0)stop("Number of high frequency lags should be a multiple of frequency")
    if(type=="A") {
        return(as.vector(apply(matrix(1:length(p),ncol=hf),2,function(x)weight(p[x],m,m))))
    }
    if(type=="B") {
        theta <- p[-1:-hf]
        return(as.vector(sapply(p[1:hf],function(beta)weight(c(beta,theta),m))))
    }
    if(type=="C") {
        return(p[1]*as.vector(apply(matrix(2:length(p),ncol=hf),2,function(x)weight(c(1,p[x]),m,m))))
    }
}


##' Gradient function for normalized exponential Almon lag weights
##'
##' Gradient function for normalized exponential Almon lag weights
##' @param p hyperparameters for Almon lag
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return the gradient matrix
##' @author Vaidotas Zemlys
##' @export
nealmon_gradient <- function(p,d,m) {
    i <- 1:d
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
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
nbeta <- function(p,d,m) {
    nbetaMT(c(p,0),d,m)
}

##' Gradient function for normalized beta probability density function MIDAS weights specification
##' Calculate gradient function for normalized beta probability density function specification of MIDAS weights.
##' @param p parameters for normalized beta probability density function
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
nbeta_gradient <- function(p,d,m) {
    nbetaMT_gradient(c(p,0),d,m)[,1:3]
}


##' Normalized beta probability density function MIDAS weights specification (MATLAB toolbox compatible)
##' Calculate MIDAS weights according to normalized beta probability density function specification. Compatible with the specification in MATLAB toolbox.
##' @param p parameters for normalized beta probability density function
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
nbetaMT <- function(p,d,m) {
    eps <- .Machine$double.eps
    xi <- (1:d - 1)/(d - 1)
    xi[1] <- xi[1]+eps
    xi[d] <- xi[d]-eps
    nb <- xi^(p[2] - 1) * (1 - xi)^(p[3] - 1)
    if(sum(nb)<eps) {
        if(abs(p[4])<eps) rep(0,d)
        else p[1]*rep(1/d,length(nb))
    } else {
        w <- (nb/sum(nb) + p[4])
        p[1]*w/sum(w)
    }
}

##' Gradient function for normalized beta probability density function MIDAS weights specification (MATLAB toolbox compatible)
##' Calculate gradient function for normalized beta probability density function specification of MIDAS weights.
##' @param p parameters for normalized beta probability density function
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
nbetaMT_gradient <- function(p,d,m) {
    eps <- .Machine$double.eps
    xi <- (1:d-1)/(d-1)
    xi[1] <- xi[1]+eps
    xi[d] <- xi[d]-eps
    nb <- xi^(p[2]-1)*(1-xi)^(p[3]-1)
    snb <- sum(nb)
    wnb <- 1+d*p[4]
    if(snb>eps) {
        nba <- nb*log(xi)
        nbb <- nb*log(1-xi)
        a <- nba/snb-nb*sum(nba)/snb^2
        b <- nbb/snb-nb*sum(nbb)/snb^2        
        cbind((nb/snb+p[4])/wnb,p[1]*a/wnb,p[1]*b/wnb,p[1]*(1-nb/snb*d)/wnb^2)
    }
    else {
       gres <- matrix(0,nrow=d,ncol=4)
       if(abs(p[4])>eps)gres[,1] <- 1/d
       gres
    }    
}

##' Almon polynomial MIDAS weights specification
##'
##' Calculate Almon polynomial MIDAS weights
##' @param p parameters for Almon polynomial weights
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
almonp <- function(p,d,m) {
    i <- 1:d
    plc <- poly(i,degree=length(p)-1,raw=TRUE) %*%p[-1]+p[1]
    as.vector(plc)
}

##' Gradient function for Almon polynomial MIDAS weights
##'
##' Calculate gradient for Almon polynomial MIDAS weights specification
##' @param p vector of parameters for Almon polynomial specification
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Vaidotas Zemlys
##' @export
almonp_gradient <- function(p,d,m) {
    i <- 1:d
    plc <- poly(i,degree=length(p)-1,raw=TRUE)
    cbind(1,plc)
}

##' Step function specification for MIDAS weights
##'
##' Step function specification for MIDAS weights
##' @param p vector of parameters
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @param a vector of increasing positive integers indicating the steps
##' @return vector of coefficients
##' @author Vaidotas Zemlys
##' @export
polystep <- function(p,d,m,a) {
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
##' @param m the frequency ratio, currently ignored
##' @param a vector of increasing positive integers indicating the steps
##' @return vector of coefficients
##' @author Vaidotas Zemlys
##' @export
polystep_gradient <- function(p,d,m,a) {
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
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
gompertzp <- function(p, d, m) {
  i <- 1:d / d
  gm <- exp(p[3] * i - p[2] * exp(p[3] * i))
  p[1] * gm / sum(gm)
}

##' Gradient function for normalized Gompertz probability density function MIDAS weights specification
##' Calculate gradient function for normalized Gompertz probability density function specification of MIDAS weights.
##' @param p parameters for normalized Gompertz probability density function
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
gompertzp_gradient <- function(p, d, m) {
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
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
nakagamip <- function(p, d, m) {
  i <- 1:d / d
  ng <- i^(2 * p[2] - 1) * exp(-p[2] / p[3] * i^2)
  p[1] * ng / sum(ng)
}

##' Gradient function for normalized Nakagami probability density function MIDAS weights specification
##' Calculate gradient function for normalized Nakagami probability density function specification of MIDAS weights.
##' @param p parameters for normalized Nakagami probability density function
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
nakagamip_gradient <- function(p, d, m) {
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
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
lcauchyp <- function(p, d, m) {
  i <- 1:d / d
  lc <- 1 / (i * ((log(i) - p[2])^2 + p[3]^2))
  p[1] * lc / sum(lc)
}

##' Gradient function for normalized log-Cauchy probability density function MIDAS weights specification
##' Calculate gradient function for normalized log-Cauchy probability density function specification of MIDAS weights.
##' @param p parameters for normalized log-Cauchy probability density function
##' @param d number of coefficients
##' @param m the frequency ratio, currently ignored
##' @return vector of coefficients
##' @author Julius Vainora
##' @export
lcauchyp_gradient <- function(p, d, m) {
  i <- 1:d / d
  lc <- 1 / (i * ((log(i) - p[2])^2 + p[3]^2))
  dp2 <- 2 * lc^2 * i * (log(i) - p[2])
  dp3 <- -2 * lc^2 * p[3] * i
  cbind(lc, p[1] * (dp2 - lc * sum(dp2) / sum(lc)),
        p[1] * (dp3 - lc * sum(dp3) / sum(lc))) / sum(lc)
}

##' HAR(3)-RV model MIDAS weights specification
##'
##' MIDAS weights for Heterogeneous Autoregressive model of Realized Volatilty (HAR-RV). It is assumed that month has 20 days.
##' @title HAR(3)-RV model MIDAS weights specification
##' @param p parameters for Almon lag
##' @param d number of the coefficients
##' @param m the frequency, currently ignored.
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Corsi, F., \emph{A Simple Approximate Long-Memory Model of Realized Volatility}, Journal of Financial Econometrics Vol. 7 No. 2 (2009) 174-196 
##' @export
harstep <- function(p,d,m) {
    if(d!=20) stop("HAR(3)-RV process requires 20 lags")
    out <- rep(0,20)
    out[1] <- p[1]+p[2]/5+p[3]/20
    out[2:5] <- p[2]/5+p[3]/20
    out[6:20] <- p[3]/20
    out
}
##' Gradient function for HAR(3)-RV model MIDAS weights specification
##'
##' MIDAS weights for Heterogeneous Autoregressive model of Realized Volatilty (HAR-RV). It is assumed that month has 20 days.
##' @title Gradient function for HAR(3)-RV model MIDAS weights specification
##' @param p parameters for Almon lag
##' @param d number of the coefficients
##' @param m the frequency, currently ignored.
##' @return vector of coefficients
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @references Corsi, F., \emph{A Simple Approximate Long-Memory Model of Realized Volatility}, Journal of Financial Econometrics Vol. 7 No. 2 (2009) 174-196 
##' @export
harstep_gradient <- function(p,d,m) {
   if(d!=20) stop("HAR(3)-RV process requires 20 lags")
   out <- matrix(0,ncol=3,nrow=d)
   out[1,1] <- 1
   out[1:5,2] <- 1/5
   out[1:20,3] <- 1/20
   out
}
