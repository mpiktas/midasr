##' Restricted MIDAS regression with I(1) regressors
##'
##' Estimate restricted MIDAS regression using non-linear least squares, when the regressor is I(1)
##'
##' @param x either formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are thearguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param user.gradient the default value is \code{FALSE}, which means that the numeric approximation of weight function gradient is calculated. If \code{TRUE} it is assumed that the R function for weight function gradient has the name of the weight function appended with \code{.gradient}. This function must return the matrix with dimensions \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the numbers of coefficients in unrestricted and restricted regressions correspondingly.
##' @param ... additional arguments supplied to optimisation function
##' @return a \code{midas_r} object which is the list with the following elements:
##' 
##' \item{coefficients}{the estimates of parameters of restrictions}
##' \item{midas.coefficientas}{the estimates of restricted coefficients of MIDAS regression}
##' \item{model}{model data}
##' \item{weights}{the restriction function(s) used in estimation.}
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
##' @rdname imidas_r
##' @seealso midas_r.midas_r
##' @examples
##' theta.h0 <- function(p, dk) {
##'   i <- (1:dk-1)/100
##'   pol <- p[3]*i + p[4]*i^2
##'  (p[1] + p[2]*i)*exp(pol)
##' }
##'
##' theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)
##' 
##' xx <- simplearma.sim(list(ar=1),1500*12,1,12)
##' y <- midas.sim(500,theta0,xx,1)
##' x <- window(xx,start=start(y))
##'
##' imr <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,start=list(x=c(-0.1,10,-10,-10)))
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
imidas_r.default <- function(x, ldata=NULL, hdata=NULL, start, Ofunction="optim", user.gradient=FALSE,...) {

    Zenv <- new.env(parent=environment(x))
        
    mt <- terms(formula(x),specials="fmls")
    vl <- as.list(attr(mt,"variables"))
    vl <- vl[-1]
    
    pl <- attr(mt,"specials")$fmls
    if(length(pl)>1) stop("Only one high frequency term is supported currently")
    fr <- vl[[pl]]
    wf <- fr[ -4:-5]
    wf[[1]] <- fr[[5]]
    for(j in 3:length(wf)) {
             wf[[j]] <- eval(wf[[j]],Zenv)
         }
    wf[[3]] <- wf[[3]]+1
    pp <- function(p,d) {
        wf[[2]] <- p        
        r <- eval(wf,Zenv)
        cumsum(r)[1:d]
    }
    mfg <- wf
    mfg[[1]] <- as.name(paste(as.character(wf[[1]]),"gradient",sep="."))
    pp.gradient <- function(p,d) {
        mfg[[2]] <- p
        r <- eval(mfg,Zenv)
        apply(r,2,cumsum)[1:d,]
    }
                              
    formula <- expandfmls(formula(x),"pp",Zenv,0)    
    cl <- match.call(expand.dots=TRUE)
    cl <- cl[names(cl)!="model"]
 
    assign("pp",pp,Zenv)
    assign("pp.gradient",pp.gradient,Zenv)
    environment(formula) <- Zenv
    cl[[2]] <- formula    
    cl[[1]] <- as.name("midas_r")
    res <- eval(cl,Zenv)
    class(res) <- c(class(res),"imidas_r")
    return(res)
}
##' Restricted MIDAS regression with I(1) regressors
##'
##' Reestimate the MIDAS regression with I(1) regressors with different starting values
##' 
##' @param x \code{imidas_r} object 
##' @param start the starting values
##' @param Ofunction a character string of the optimisation function to use. The default value is to use the function of previous optimisation.
##' @param ... further arguments to optimisation function. If none are supplied, the arguments of previous optimisation are used.
##' @return \code{imidas_r} object
##' @method imidas_r imidas_r
##' @seealso imidas_r
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
imidas_r.imidas_r <- function(x,start=coef(x),Ofunction=x$argmap.opt$Ofunction,...) {
    midas_r.midas_r(x,start=start,Ofunction=Ofunction,...)
}

##Function for expanding the formula in I(1) case
expandfmls <- function(expr,wfun,Zenv,diff=0,truncate=FALSE) {
    if(length(expr)==3) {
        expr[[2]] <- expandfmls(expr[[2]],wfun,Zenv,diff,truncate)
        expr[[3]] <- expandfmls(expr[[3]],wfun,Zenv,diff,truncate)
    }
    if(length(expr)==5) {
        if(expr[[1]]==as.name("fmls")) {
            rr <- modifyfmls(expr,wfun,Zenv,diff)
            if(truncate) {
                return(rr[[2]])
            }
            else {
                return(rr)
            }
        }
        else return(expr)
    }
    return(expr)
}

##Substitute fmls with dmls
modifyfmls <- function(expr,wfun,Zenv,diff) {
     res <- expression(a+b)[[1]]
     t2 <- expression(mls(x,a,b))[[1]]
     nol <- eval(expr[[3]],Zenv)
     expr[[3]] <- nol+diff
     m <- eval(expr[[4]],Zenv)
     expr[[1]] <- as.name("dmls")
     t2[[3]] <- nol+1+diff
     t2[[4]] <- m
     expr[[5]] <- as.name(wfun)
     res[[2]] <-expr
     res[[3]] <- t2
     res
}
