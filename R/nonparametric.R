##' Estimates non-parametric MIDAS regression
##'
##' Estimates non-parametric MIDAS regression accodring Breitung et al.
##' 
##' @title Estimate non-parametric MIDAS regression
##' @param x formula specifying MIDAS regression 
##' @param data a named list containing data with mixed frequencies
##' @param lambda smoothing parameter, defaults to \code{NULL}, which means that it is chosen by minimising AIC.
##' @return a \code{midas_r_np} object
##' @author Vaidotas Zemlys
##' @references Breitung J, Roling C, Elengikal S (2013). \emph{The statistical content and empirical testing of the MIDAS restrictions} Working paper, URL http://www.ect.uni-bonn.de/mitarbeiter/ joerg-breitung/npmidas.
##' @export
##' @import Matrix
##' @examples
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##' midas_r_np(y~fmls(x,12,12)-1)
midas_r_np <- function(x,data,lambda=NULL) {
    Zenv <- new.env(parent=environment(x))

    if(missing(data)) {
        ee <- NULL
    }
    else {
        ee <- data_to_env(data)
        parent.env(ee) <- parent.frame()
    }
    
    assign("ee",ee,Zenv)
    x <- as.formula(x)
    cl <- match.call()    
    mf <- match.call(expand.dots = FALSE)
    mf$x <- x
    m <- match(c("x", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- as.name("na.omit")
    names(mf)[c(2,3,4)] <- c("formula","data","na.action")
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")

    if(attr(mt,"intercept")==1)stop("Models with intercept are not supported yet")

   # args <- list(...)
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)    

    k <- ncol(X)
           
    D <- bandSparse(k-2,k,c(0,1,2),diagonals=list(rep(1,k-2),rep(-2,k-2),rep(1,k-2)))
    DD <- crossprod(D)
    ol <- opt_lambda(y,X,DD,lambda)

    fit <- X%*%ol$beta
    res <- y-fit
    
    out <- list(coefficients=as.numeric(ol$beta),
                midas.coefficients=as.numeric(ol$beta),
                model=cbind(y,X),
                call=cl,        
                terms=mt,
                fitted.values=as.numeric(fit),
                residuals=as.numeric(res),
                lambda=ol$lambda,
                klambda=ol$klambda,
                AIC=ol$AIC,
                Zenv=Zenv
                )
    
    class(out) <- "midas_r_np"
    out
}

##' @export
##' @method AIC midas_r_np
AIC.midas_r_np <- function(object, ..., k) {
    object$AIC(object$lambda)
}
##' @export
##' @method forecast midas_r_np
forecast.midas_r_np <- forecast.midas_r

##' @export
##' @method predict midas_r_np
predict.midas_r_np <- predict.midas_r

##' @export
##' @method print midas_r_np
print.midas_r_np <- function(x,...) {
    cat("Nonparametric MIDAS regression model")
    cat("\nThe smoothing parameter: ", x$lambda)
    cat("\nThe effective number of parameters:", x$klambda)
    cat("\nAIC of the model: ",AIC(x))
    cat("\nResidual standard error: ", sqrt(mean(residuals(x)^2)),"\n")
}

##' @export
##' @method summary midas_r_np
summary.midas_r_np <- function(object,...) {
   print(object)
}

opt_lambda <- function(y,X,DD,lambda) {
    n <- length(y)
    XX <- crossprod(X)
    Xy <- crossprod(X,y)
    tX <- t(X)
    AIC <- function(lambda) {
            Qlambda <- XX+lambda*n*DD   
            klambda <- sum(diag(X%*%solve(Qlambda,tX)))
            beta <- solve(Qlambda,Xy)
            res <- y-X%*%beta
            log(sum(res^2))+2*(klambda+1)/(n-klambda-2)
        }
    if(is.null(lambda)) {       
        opt <- optim(1,AIC,method="BFGS")
        lambda <- opt$par
    }
    Qlambda <- XX+lambda*n*DD   
    klambda <- sum(diag(X%*%solve(Qlambda,tX)))
    beta <- solve(Qlambda,Xy)    
    AICl <- AIC(lambda)
    list(beta=beta,klambda=klambda,lambda=lambda,AIC=AIC)
    
}
