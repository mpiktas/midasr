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

    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    term.labels <- attr(mt,"term.labels") 

    rfd <- vector("list",length(terms.lhs))
    yname <- all.vars(mt[[2]])

    for(i in 1:length(rfd)) {
        fr <- terms.lhs[[i]]
        #This a behaviour of R we rely on. It might be non-standard one.
        fun <- as.character(fr)[1]
        rfd[[i]] <- if(fun %in% c("fmls","mls","dmls")) {
            lags <- eval(fr[[3]],Zenv)
            nm <-as.character(fr[[2]])            
            nol <- switch(fun,
                              fmls = lags+1,
                              dmls = lags+1,
                              mls = length(lags)
                              )
            if(nol<3 & nm!=yname)stop("For nonparametric MIDAS you need at least 3 high frequency lags")          
           
            wlab <- ifelse(nm==yname,"",nm)
            list(length=nol,name=nm,wlabel=wlab,weight=function(p)p)
        } else {
            list(length=1,name=term.labels[i],wlabel="",weight=function(p)p)
        }
    }
    
    if(attr(mt,"intercept")==1) {
        rfd <- c(list(list(length=1,name="(Intercept)",wlabel="",weight=function(p)p)),rfd)
    }


    names(rfd) <- sapply(rfd,"[[","name")
    rf <- lapply(rfd,"[[","weight")
    names(rf) <- sapply(rfd,"[[","name")
    
    weight_names <- sapply(rfd,"[[","wlabel")
    weight_inds <- which(weight_names!="")
    weight_names <- weight_names[weight_names!=""]
    lengths <- sapply(rfd,"[[","length")
    
    build_indices <- function(ci,nm) {
        inds <- cbind(c(1,ci[-length(ci)]+1),ci)
        inds <- apply(inds,1,function(x)list(x[1]:x[2]))
        inds <- lapply(inds,function(x)x[[1]])
        names(inds) <- nm
        inds
    }
    
    pinds <- build_indices(cumsum(lengths),names(rfd))
        
    if(length(weight_names)>1)stop("Only one non-autoregressive mixed frequency term is currently supported")
    
    
    resplace <- pinds[[weight_names]][1]     
    rno <- rfd[[weight_names]]$length
    
   # args <- list(...)
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)    
    
    k <- ncol(X)
    
    D <- bandSparse(rno-2,k,resplace-1+c(0,1,2),diagonals=list(rep(1,rno-2),rep(-2,rno-2),rep(1,rno-2)))
  
    DD <- crossprod(D)
    ol <- opt_lambda(y,X,DD,lambda)

    fit <- X%*%ol$beta
    res <- y-fit
    
    cf <- as.numeric(ol$beta)
    names(cf) <- names(unlist(pinds))
    
    out <- list(coefficients=cf,
                midas.coefficients=cf,
                model=cbind(y,X),
                call=cl,        
                terms=mt,
                fitted.values=as.numeric(fit),
                residuals=as.numeric(res),
                param.map=pinds,
                weights=rf[weight_inds],
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
    cat("Nonparametric MIDAS regression model",paste0("(",nrow(x$model)), "low frequency observations)")
    cat("\nFormula: ", deparse(formula(terms(x))))    
    cat("\nThe smoothing parameter: ", x$lambda)
    cat("\nThe effective number of parameters:", x$klambda)
    cat("\nAIC of the model: ",AIC(x))
    cat("\nRoot mean squared error: ", sqrt(mean(residuals(x)^2)),"\n")
}

##' @export
##' @method summary midas_r_np
summary.midas_r_np <- function(object,...) {
  print(object,...)
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
