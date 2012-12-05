##' Check whether non-linear least squares restricted MIDAS regression problem has converged
##'
##' Computes the gradient and hessian of the optimisation function of restricted MIDAS regression and checks whether the conditions of local optimum are met. Numerical estimates are used.
##' @param x \code{\link{midas_r}} object
##' @param tol a tolerance, values below the tolerance are considered zero
##' @return a list with gradient, hessian of optimisation function and convergence message
##' @rdname deriv_tests
##' @seealso midas_r
##' @export
##' @author Vaidotas Zemlys
deriv_tests<- function(x,tol=1e-6) UseMethod("deriv_tests")

#' @rdname deriv_tests
#' @method deriv_tests midas_r
#' @export
deriv_tests.midas_r <- function(x,tol=1e-6) {
    gr <- x$gradient(coef(x))
    hess <- x$hessian(coef(x))
    egh <- eigen(hess)$values

    first <- sum(abs(gr))<tol
    second <- !any(egh<0) & !any(abs(egh)<tol)
    list(first=first,second=second,gradient=gr,eigenval=egh)
}

##' Gradient function of coefficient function in MIDAS regression
##'
##' Returns gradient function of coefficient function
##' 
##' @param x \code{\link{midas_r}} object
##' @param gr a function or a list of functions for each restriction.
##' @param ... arguments supplied to gradient function(s)
##' @return a function
##' @author Vaidotas Zemlys
##' @details MIDAS regression can be written in the following form:
##'
##' \deqn{y=Xf(\eta)+u}
##'
##' This function calculates
##'
##' \deqn{\frac{\partial f}{\partial \eta}}
##' 
gradD <- function(x,gr=NULL,...) {
    rf <- lapply(x$param.map,function(x)return(function(p,...)p)) 
    names(rf) <- names(x$param.map)

    restr.name <- names(x$restrictions)   
    
    if(is.null(gr)) {
        rf[names(x$restrictions)] <- x$restrictions
        all_coef <- function(p) {
             pp <- lapply(x$param.map,function(xx)p[xx])     
             res <- mapply(function(fun,param)fun(param),rf,pp,SIMPLIFY=FALSE)
             unlist(res)
         }

        gr <- function(p,...)jacobian(all_coef,p)
    }
    else {
        if(is.function(gr) & length(x$restrictions)==1) {
            grf <- list(gr)
            names(grf) <- names(x$restrictions)
        }
        else {
            if(is.null(names(gr)) | !is.list(gr)) stop("Argument gr must be a named list")
            if(any(!names(gr) %in% names(x$restrictions))) stop("The names of the list must coincide with the names of the restriction functions")
            grf <- lapply(x$param.map,function(x)return(function(p)return(matrix(1))))
            grf[names(x$restrictions)] <- gr[names(x$restrictions)]
        }
        pp0 <- lapply(x$param.map,function(xx)coef(x)[xx])            
        grmat0 <- mapply(function(fun,param,MoreArgs=...)fun(param,...),grf,pp0,SIMPLIFY=FALSE)
        colnos <- sapply(grmat0,ncol)
        rownos <- sapply(grmat0,nrow)
        np <- length(colnos)
        ccol <- cumsum(colnos)
        rrow <- cumsum(rownos)
        pinds <- cbind(c(1,rrow[-np]+1),rrow,
                       c(1,ccol[-np]+1),ccol)
        pinds <- apply(pinds,1,function(x)list(row=x[1]:x[2],col=x[3]:x[4]))
                
        gr <- function(p,...) {
            pp <- lapply(x$param.map,function(x)p[x])
            grmat <- mapply(function(fun,param,MoreArgs=...)fun(param,...),grf,pp,SIMPLIFY=FALSE)
            if(length(grmat)==1) {
                res <- grmat[[1]]
            }
            else {
                res <- matrix(0,nrow=sum(rownos),ncol=sum(colnos))
                mapply(function(m,ind){res[ind$row,ind$col] <- m},
                       grmat,pinds)
                   
            }
            res
        }
    }
    function(p)gr(p,...)
}
