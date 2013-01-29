##' Restricted MIDAS regression with I(1) regressors
##'
##' Estimate restricted MIDAS regression using non-linear least squares, when the regressor is I(1)
##'
##' @param x either formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are thearguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param gradient the gradient of the restriction function. Must return the matrix with dimensions \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the numbers of coefficients in unrestricted and restricted regressions correspondingly. Default value is \code{NULL}, which means that the numeric approximation of gradient is calculated.
##' @param ... additional arguments supplied to optimisation function
##' @return a \code{midas_r} object which is the list with the following elements:
##' 
##' \item{coefficients}{the estimates of parameters of restrictions}
##' \item{midas.coefficientas}{the estimates of restricted coefficients of MIDAS regression}
##' \item{model}{model data}
##' \item{restrictions}{the restriction function(s) used in estimation.}
##' \item{unrestricted}{unrestricted regression estimated using \code{\link{midas_u}}}
##' \item{param.map}{parameter map for optimisation function}
##' \item{fn0}{optimisation function for non-linear least squares problem solved in restricted MIDAS regression}
##' \item{rhs}{the function which evaluates the right-hand side of the MIDAS regression}
##' \item{allcoef}{the function which evaluates the restricted coefficientsof MIDAS regression}
##' \item{opt}{the output of optimisation procedure}
##' \item{argmap.opt}{the list containing the name of optimisation function together with arguments for optimisation function}
##' \item{start.opt}{the starting values used in optimisation}
##' 
##' 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname midas_r
##' @seealso midas_r.midas_r
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
##' x <- simplearma.sim(list(ar=1),1500*12,1,12)
##'
##' ##Simulate the response variable
##' y <- midas.sim(500,theta0,x,1)
##'
##' ##Remove unnecessary history of x
##' x <- window(x,start=start(y))
##' 
##' ##Fit restricted model
##' mr <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(x=c(-0.1,10,-10,-10)))
##'
##' ##Include intercept and trend in regression
##'
##' mr.it <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)+trend,data.frame(y=y,trend=1:500),data.frame(x=x),start=list(x=c(-0.1,10,-10,-10)))
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
##' It is assumed that \eqn{x} is a I(1) process, hence the special transformation is made. After the transformation \link{midas_r} is used for estimation.
##'
##' MIDAS regression involves times series with different frequencies. 
##'
##' The restriction function must return the restricted coefficients of
##' the MIDAS regression.
##' 
##' @export
imidas_r <- function(x,...)UseMethod("imidas_r")

is.imidas_r <- function(x) inherits(x,"imidas_r")

#' @rdname imidas_r
#' @method imidas_r default
#' @export
imidas_r.default <- function(x, ldata=NULL, hdata=NULL, start, Ofunction="optim", gradient=NULL,...) {

    Zenv <- new.env(parent=environment(x))

    mt <- terms(formula(x),specials="fmls")
    vl <- as.list(attr(mt,"variables"))
    vl <- vl[-1]
    
    fr <- vl[[attr(mt,"specials")$fmls]]
    mf <- fr[ -4:-5]
    mf[[1]] <- fr[[5]]
    for(j in 3:length(mf)) {
             mf[[j]] <- eval(mf[[j]],Zenv)
         }        
    
    pp <- function(p,d) {
        mf[[2]] <- p
        mf[[3]] <- d
        r <- eval(mf,Zenv)
        cumsum(r)
    }
    mfg <- mf
    mfg[[1]] <- as.name(paste(as.character(mf[[1]]),"gradient",sep="."))
    pp.gradient <- function(p,d) {
        mfg[[2]] <- p
        mfg[[3]] <- d
        r <- eval(mfg,Zenv)
        apply(r,2,cumsum)
    }

    formula <- expandfmls(formula(x),"pp",Zenv)
    cl <- match.call(expand.dots=TRUE)
    cl[[2]] <- formula
    cl[[1]] <- as.name("midas_r")
    eval(cl,Zenv)
}


##Function for expanding the formula in I(1) case
expandfmls <- function(expr,wfun,Zenv) {
    if(length(expr)==3) {
        expr[[2]] <- expandfmls(expr[[2]],wfun,Zenv)
        expr[[3]] <- expandfmls(expr[[3]],wfun,Zenv)
    }
    if(length(expr)==5) {
        if(expr[[1]]==as.name("fmls")) {
            res <- expression(a+b)[[1]]
            t2 <- expression(mls(x,a,b))[[1]]
            nol <- eval(expr[[3]],Zenv)
            m <- eval(expr[[4]],Zenv)
            expr[[1]] <- as.name("dmls")
            t2[[3]] <- nol+1
            t2[[4]] <- m
            expr[[5]] <- as.name(wfun)
            res[[2]] <-expr
            res[[3]] <- t2
            return(res)
        }
        else return(expr)
    }
    return(expr)
}
