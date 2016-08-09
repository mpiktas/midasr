
midas_qr <- function(formula, data, tau = 0.5, start, Ofunction="nlrq", weight_gradients=NULL,...) {
    mf <- match.call(expand.dots = TRUE)
    
    mff <- mf
    mff[[1]] <- as.name("midas_r")
    mff[[match("tau",names(mff))]] <- NULL
    mff$Ofunction <- "dry_run"
    ##Do a nicer argument matching
    
    x <- eval(mff, envir=parent.frame())
    x$argmap_opt$Ofunction <- Ofunction
    x$tau <- tau
    x$argmap_opt$tau <- NULL
    class(x) <- "midas_qr"
    midas_qr.fit(x)
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
        y <- x$model[,1]
        args$formula <- formula(y~rhs(p))
        args$start <- list(p=x$start_opt)
        args$tau <- x$tau
        args$method <- NULL
        browser()
        opt <- try(do.call("nlrq",args),silent=TRUE)
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