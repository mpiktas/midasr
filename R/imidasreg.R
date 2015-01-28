##' Restricted MIDAS regression with I(1) regressors
##'
##' Estimate restricted MIDAS regression using non-linear least squares, when the regressor is I(1)
##'
##' @param formula formula for restricted MIDAS regression. Formula must include \code{\link{fmls}} function
##' @param data a named list containing data with mixed frequencies
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are the arguments passed to this function. The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param weight_gradients a named list containing gradient functions of weights. The weight gradient function must return the matrix with dimensions
##' \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the number of coefficients in unrestricted and restricted regressions correspondingly.
##' The names of the list should coincide with the names of weights used in formula.
##' The default value is NULL, which means that the numeric approximation of weight function gradient is calculated. If the argument is not NULL, but the
##' name of the weight used in formula is not present, it is assumed that there exists an R function which has  
##' the name of the weight function appended with \code{.gradient}. 
##' @param ... additional arguments supplied to optimisation function
##' @return a \code{midas_r} object which is the list with the following elements:
##' 
##' \item{coefficients}{the estimates of parameters of restrictions}
##' \item{midas_coefficients}{the estimates of MIDAS coefficients of MIDAS regression}
##' \item{model}{model data}
##' \item{unrestricted}{unrestricted regression estimated using \code{\link{midas_u}}}
##' \item{term_info}{the named list. Each element is a list with the information about the term, such as its frequency, function for weights, gradient function of weights, etc.}
##' \item{fn0}{optimisation function for non-linear least squares problem solved in restricted MIDAS regression}
##' \item{rhs}{the function which evaluates the right-hand side of the MIDAS regression}
##' \item{gen_midas_coef}{the function which generates the MIDAS coefficients of MIDAS regression}
##' \item{opt}{the output of optimisation procedure}
##' \item{argmap_opt}{the list containing the name of optimisation function together with arguments for optimisation function}
##' \item{start_opt}{the starting values used in optimisation}
##' \item{start_list}{the starting values as a list}
##' \item{call}{the call to the function}
##' \item{terms}{terms object}
##' \item{gradient}{gradient of NLS objective function}
##' \item{hessian}{hessian of NLS objective function}
##' \item{gradD}{gradient function of MIDAS weight functions} 
##' \item{Zenv}{the environment in which data is placed}
##' \item{use_gradient}{TRUE if user supplied gradient is used, FALSE otherwise}
##' \item{nobs}{the number of effective observations}
##' \item{convergence}{the convergence message}
##' \item{fitted.values}{the fitted values of MIDAS regression}
##' \item{residuals}{the residuals of MIDAS regression}
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
##' xx <- ts(cumsum(rnorm(600*12)), frequency = 12)
##'
##' ##Simulate the response variable
##' y <- midas_sim(500, xx, theta0)
##'
##' x <- window(xx, start=start(y))
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
imidas_r <- function(formula, data, start, Ofunction = "optim", weight_gradients = NULL,...) {

    Zenv <- new.env(parent=environment(formula))
    
    mt <- terms(formula(formula),specials="fmls")


    if(missing(data)) {
        nms <- all.vars(mt)        
        insample <- lapply(nms,function(nm)eval(as.name(nm),Zenv))
        names(insample) <- nms
        insample <- insample[!sapply(insample,is.function)]
    }
    else {
        insample <- data_to_list(data)
    }
    
    vl <- as.list(attr(mt,"variables"))
    vl <- vl[-1]
    
    pl <- attr(mt,"specials")$fmls
    if(length(pl)>1) stop("Only one high frequency term is supported currently")
    fr <- vl[[pl]]

    ##Do an ugly hack
    xname <- as.character(fr[[2]])
    insample <- c(list(insample[[xname]]),insample)
    names(insample)[1] <- ".Level"
   
    assign(".IMIDASdata",insample,envir=Zenv)
    
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
                              
    formula <- expandfmls(formula(formula),"pp",Zenv,0)    
    cl <- match.call(expand.dots=TRUE)
    cl <- cl[names(cl)!="model"]
 
    assign("pp",pp,Zenv)
    assign("pp.gradient",pp.gradient,Zenv)
    environment(formula) <- Zenv
    cl[[2]] <- formula    
    cl[[1]] <- as.name("midas_r")
    cl$data <- as.name(".IMIDASdata")
    res <- eval(cl,Zenv)
    class(res) <- c(class(res),"imidas_r")
    return(res)
}

## @method update imidas_r
## @export
#update.imidas_r <- update.midas_r

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
     t2[[2]] <- as.name(".Level")
     t2[[3]] <- nol+1+diff
     t2[[4]] <- m
     expr[[5]] <- as.name(wfun)
     res[[2]] <-expr
     res[[3]] <- t2
     res
}
