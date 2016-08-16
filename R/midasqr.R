
midas_qr <- function(formula, data, tau = 0.5, start, Ofunction="nlrq", weight_gradients=NULL,...) {
    Zenv <- new.env(parent=environment(formula))
    
    if(missing(data)) {
        ee <- NULL
    }
    else {
        ee <- data_to_env(data)
        parent.env(ee) <- parent.frame()
    }
    
    if(missing(start)) {
        stop("Please supply starting values.")
    } 
    
    assign("ee",ee,Zenv)
    formula <- as.formula(formula)
    cl <- match.call()    
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- formula
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- as.name("na.omit")
    names(mf)[c(2,3,4)] <- c("formula","data","na.action")
    
    ##Add check for dynamic terms. 
    ##They are not supported as they do not make sense in quantile regression
    itr <- checkARstar(terms(eval(mf[[2]], Zenv)))
    if(!is.null(itr$lagsTable)) mf[[2]] <- itr$x
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
    args <- list(...)
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)
    
    #Save ts information
    if(is.null(ee)) { 
        yy <- eval(formula[[2]], Zenv)
    }else {
        yy <- eval(formula[[2]], ee)
    }
    
    y_index <- 1:length(yy) 
    if(!is.null(attr(mf, "na.action"))) {
        y_index <- y_index[-attr(mf, "na.action")]
    }
    if(length(y_index)>1) {
        if(sum(abs(diff(y_index) - 1))>0) warning("There are NAs in the middle of the time series")                
    }
    ysave <- yy[y_index]
    
    if(inherits(yy, "ts")) {
        class(ysave) <- class(yy)
        attr(ysave, "tsp") <- c(time(yy)[range(y_index)], frequency(yy))
    }
    
    prepmd <- prepmidas_r(y,X,mt,Zenv,cl,args,start,Ofunction,weight_gradients,itr$lagsTable, guess_start = TRUE, tau = tau)
    
    prepmd <- c(prepmd, list(lhs = ysave, lhs_start= y_start, lhs_end = y_end, ))
    
    class(prepmd) <- "midas_qr"
    
    midas_qr.fit(prepmd)    
    
}

midas_qr.fit <- function(x) {
    args <- x$argmap_opt
    function.opt <- args$Ofunction
    args$Ofunction <- NULL
    if(function.opt=="nlrq") {
        rhs <- x$rhs
        if(x$use_gradient) {
            orhs <- rhs
            rhs <- function(p) {
                res <- orhs(p)
                attr(res,"gradient") <- x$model[,-1]%*%x$gradD(p)
                res
            }
        }
        z <- x$model[,1]
        args$formula <- formula(z~rhs(p))
        args$start <- list(p=x$start_opt)
        args$tau <- x$tau
        nlrq(z~rhs(p),start=list(p=x$start_opt))
        opt <- try(do.call("nlrq",args),silent = TRUE)
        if(inherits(opt,"try-error")) {
            stop("The optimisation algorithm of MIDAS regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        par <- coef(opt)
        names(par) <- names(coef(x))
        x$convergence <- opt$convInfo$stopCode
    }
    x$opt <- opt
    x$coefficients <- par
    names(par) <- NULL
    x$midas_coefficients <- x$gen_midas_coef(par)
    x$fitted.values <- as.vector(x$model[,-1]%*%x$midas_coefficients)
    x$residuals <- as.vector(x$model[,1]-x$fitted.values)
    x
    
}