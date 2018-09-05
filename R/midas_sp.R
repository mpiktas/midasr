##' Semi-parametric MIDAS regression
##'
##' Estimate semi-parametric MIDAS regression using non-linear least squares.
##'
##' @param formula formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param data a named list containing data with mixed frequencies
##' @param bws a bandwith specification, starting values for bandwith, if unset the default will be used.
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param degree the degree of local polynomial. 0 corresponds to local-constant, 1 local-linear. For univariate models higher values can be provided.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are the arguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param ... additional arguments supplied to optimisation function
##' @return a \code{midas_sp} object which is the list with the following elements:
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
##' \item{nobs}{the number of effective observations}
##' \item{convergence}{the convergence message}
##' \item{fitted.values}{the fitted values of MIDAS regression}
##' \item{residuals}{the residuals of MIDAS regression}
##' 
##' @author Virmantas Kvedaras, Vaidotas Zemlys-Balevičius
##' @rdname midas_sp
##' @details Given MIDAS regression:
##'
##' \deqn{y_t = \sum_{j=1}^p\alpha_jy_{t-j} +\sum_{i=0}^{k}\sum_{j=0}^{l_i}\beta_{j}^{(i)}x_{tm_i-j}^{(i)} + u_t,}
##' 
##' estimate the parameters of the restriction
##'
##' \deqn{\beta_j^{(i)}=g^{(i)}(j,\lambda).}
##'
##' Such model is a generalisation of so called ADL-MIDAS regression. It is not required that all the coefficients should be restricted, i.e the function \eqn{g^{(i)}}
##' might be an identity function. The regressors \eqn{x_\tau^{(i)}} must be of higher
##' (or of the same) frequency as the dependent variable \eqn{y_t}. 
##'
##' @importFrom stats as.formula formula model.matrix model.response terms lsfit time
##' @importFrom zoo index index2char
##' @importFrom Formula Formula
##' @export
midas_sp <- function(formula, data, bws, start, degree = 1,  Ofunction="optim", ...) {
    
    Zenv <- new.env(parent = environment(formula))
    
    if (missing(data)) {
        ee <- NULL
    }
    else {
        ee <- data_to_env(data)
        parent.env(ee) <- parent.frame()
    }
    
    if (missing(start)) {
        stop("Please supply starting values.")
    } 
    
    assign("ee",ee,Zenv)
    formula <- Formula(formula)
    environment(formula) <- Zenv
    cl <- match.call()    
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- formula
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- as.name("na.omit")
    names(mf)[c(2,3,4)] <- c("formula","data","na.action")
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
    args <- list(...)
    y <- as.numeric(model.response(mf, "numeric"))
    X <- model.matrix(formula, data = mf, rhs = 1)
    if(length(attr(formula,"rhs")) > 1) {
        Z <- model.matrix(formula, data = mf, rhs = 2)
    } else Z <- NULL
    
    #Save ts/zoo information
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
    
    if(inherits(yy,c("zoo","ts"))) {
        y_start <- index2char(index(ysave)[1], frequency(ysave))
        y_end <- index2char(index(ysave)[length(ysave)], frequency(ysave))
    } else {
        y_start <- y_index[1]
        y_end <- y_index[length(y_index)]
    }
    
    prepmd <- prep_midas_sp(y, X, Z, bws, degree, formula, Zenv,cl,args,start,Ofunction)
    
    prepmd <- c(prepmd, list(lhs = ysave, lhs_start = y_start, lhs_end = y_end))
    
    class(prepmd) <- "midas_sp"
    
    midas_nlpr.fit(prepmd)    
}

## Prepare necessary objects for fitting of the non-linear parametric MIDAS regression
##
## y the response
## X the model matrix
## mt the terms of the formula
## Zenv the environment to evaluate the formula
## cl call of the function
## args additional argument
## start starting values
## Ofunction the optimisation function
## weight_gradients a list of gradient functions for weights
## lagsTable the lagstable from checkARstar
## unrestricted the unrestricted model
## guess_start if TRUE, get the initial values for non-MIDAS terms via OLS, if FALSE, initialize them with zero.
prep_midas_sp <- function(y, X, Z, bws, degree, f, Zenv, cl, args, start, Ofunction,  guess_start = TRUE) {
    
    
    start <- start[!sapply(start,is.null)]
    
    if (!is.null(args$guess_start)) {
        guess_start <- args$guess_start
        args$guess_start <- NULL
    }  
    if(is.null(Z)) {
        ##We are having pure SI model.
        Z <- X
        X <- 0
        mt1 <- NULL
        mt2 <- terms(f, rhs = 1)
    } else {
        mt1 <- terms(f, rhs = 1)
        mt2 <- terms(f, rhs = 2)
        
        if (attr(mt1,"intercept") == 1)  {
            X <- X[, -1, drop = FALSE]
        }
        
    }

    if (attr(mt2,"intercept") == 1)  {
        Z <- Z[, -1, drop = FALSE]
    }
    
    bws_length <- length(bws)
    
    if(!is.null(mt1)) {
        terms.rhs1 <- as.list(attr(mt1,"variables"))[-2:-1]
        rfd1 <- lapply(terms.rhs1, dterm_nlpr,  Zenv = Zenv)
        term_names1 <- sapply(rfd1,"[[","term_name")
        names(rfd1) <- term_names1
        rf1 <- lapply(rfd1,"[[","weight")
   
        start_default1 <- lapply(rfd1,"[[","start")
        start_default1[intersect(names(start), term_names1)] <- start[intersect(names(start), term_names1)]
        np1 <- cumsum(sapply(start_default1,length))
        pinds1 <- build_indices(np1,names(start_default1))
        fake_x_coef <- all_coef_full(unlist(start_default1), pinds1, rf1) 
        npx <- cumsum(sapply(fake_x_coef,length))
        xinds <- build_indices(npx,names(start_default1))
        pinds1 <- lapply(pinds1, function(x) x + bws_length)
        cf1 <- function(p) unlist(all_coef_full(p, pinds1, rf1))
        p2_offset <- max(pinds1[[length(pinds1)]])
    } else {
        cf1 <- function(p) return(0)
        p2_offset <- bws_length
        start_default1 <- NULL
    }
    
    terms.rhs2 <- as.list(attr(mt2,"variables"))[-2:-1]
    rfd2 <- lapply(terms.rhs2, dterm_nlpr,  Zenv = Zenv)
    term_names2 <- sapply(rfd2,"[[","term_name")
    names(rfd2) <- term_names2
    rf2 <- lapply(rfd2,"[[","weight")
    weight_names2 <- sapply(rfd2,"[[","weight_name")
    weight_inds2 <- which(weight_names2 != "")
    weight_names2 <- names(rf2)[weight_names2 != ""]
    
    start_fake2 <- lapply(rfd2,"[[","start")
    start_fake2[intersect(names(start), weight_names2)] <- start[intersect(names(start), weight_names2)]
    np_fake2 <- cumsum(sapply(start_fake2,length))
    pinds_fake2 <- build_indices(np_fake2,names(start_fake2)) 
    fake_z_coef <- all_coef_full(unlist(start_fake2), pinds_fake2, rf2) 
    npz <- cumsum(sapply(fake_z_coef,length))
    zinds <- build_indices(npz,names(start_fake2))
   
    start_default2 <- start_fake2[weight_names2]
    if (length(start_default2) > 0) {
        np2 <- cumsum(sapply(start_default2,length))
        pinds2 <- build_indices(np2,names(start_default2))
        pinds2 <- lapply(pinds2, function(x) x + p2_offset)
    } else start_default2 <- NULL
    
    cf0 <- function(p) p[1:bws_length]
   
    
    rhs_cv <- function(p) {
        h <- cf0(p)
        xi <- as.numeric(X %*% cf1(p))
        zi <- do.call("rbind", lapply(term_names2, function(t2) {
            if (t2 %in% weight_names2) Z[, zinds[[t2]]] %*% rf2[[t2]](p[pinds2[[t2]]])
            else Z[, zinds[[t2]], drop = FALSE]
        }))
        if (ncol(zi) == 1) zi <- as.numeric(zi)
        u <- y - xi
        xi + cv_np(u, zi, h, degree) 
    }
    
    fn0 <- function(p,...) {
        r <- y - rhs_cv(p)
        sum(r^2)
    }
    
    hess <- function(x)numDeriv::hessian(fn0,x)
    
    control <- c(list(Ofunction = Ofunction),args)
    ##Override default method of optim. Use BFGS instead of Nelder-Mead
    if (!("method" %in% names(control)) & Ofunction == "optim") {        
        control$method <- "BFGS"
    }
    if(bws_length > 1) names(bws) <- paste0("bw", 1:bws_length)
    else names(bws) <- "bw"
    
    if(is.null(mt1)) {
        starto <- c(bws, unlist(start_default2)) 
        term_info1 <- NULL   
        model_matrix_map = list(y = 1,
                                X = NULL,
                                Z = 1 + 1:ncol(Z))
        model <- cbind(y, Z)
    } else {
        starto <- c(bws, unlist(start_default1), unlist(start_default2)) 
        term_info1 <- mapply(function(l, pind, xind) {
            l$coef_index <- pind
            l$midas_coef_index <- xind
            l$term_type <- "X"
            l
        },rfd1, pinds1, xinds, SIMPLIFY = FALSE)
        model_matrix_map = list(y = 1,
                                X = 1 + 1:ncol(X),
                                Z = ncol(X) + 1 + 1:ncol(Z))
        model <- cbind(y, X, Z)
    }
    
    term_info2 <- lapply(term_names2, function(t2) {
        res <- rfd2[[t2]]
        if(t2 %in% weight_names2) {
            res$coef_index <- pinds2[[t2]]
            res$midas_coef_index <- zinds[[t2]]
        } else {
            res$midas_coef_index <- zinds[[t2]]
        }
        res$term_type <- "Z"
        res
    })
    
    names(term_info2) <- term_names2
    term_info <- c(term_info1, term_info2)
    
    
    list(coefficients = starto,
         model = model, 
         fn0 = fn0,
         rhs = rhs_cv,
         rhs_cv = rhs_cv,
         coef_X  = cf1,
         model_matrix_map = model_matrix_map,
         opt = NULL,
         argmap_opt = control,
         start_opt = starto,
         start_list = c(list(bws = bws), start1 = start_default1, start2 = start_default2),
         call = cl,
         terms = mt1,
         terms2 = mt2,
         bws = bws,
         term_info = term_info,
         hessian = hess,
         Zenv = Zenv,
         nobs = nrow(Z),
         degree = degree)   
}


##' MIDAS Single index regression
##'
##' Function for fitting SI MIDAS regression without the formula interface
##' @param y model response
##' @param X prepared matrix of high frequency variable lags for MMM term
##' @param p.ar length of AR part
##' @param weight the weight function
##' @param degree the degree of local polynomial
##' @param start_bws the starting values bandwith
##' @param start_x the starting values for weight function
##' @param start_ar the starting values for AR part. Should be the same length as \code{p}
##' @param method a method passed to \link{optimx}
##' @param ... additional parameters to \link{optimx}
##' @return an object similar to \code{midas_r} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' 
##' @import optimx
##' @importFrom stats na.omit 
##'
##' @export
##' 
midas_si_plain <- function(y, X, p.ar = NULL, weight, degree = 1, start_bws, start_x, start_ar = NULL, method = c("Nelder-Mead"), ...) {
    d <- ncol(X)
    
    yy <- NULL 
    if(!is.null(p.ar)) {
        p.ar <- as.integer(p.ar)
        if (length(start_ar)!=p.ar) stop("The length of starting values vector for AR terms must the same as the number of AR terms")
        yy <- as.matrix(mls(y, 1:p.ar, 1))
    } 
    
    model <- na.omit(cbind(y,X,yy))
    
    y <- model[,1]
    X <- model[, 2:(ncol(X) + 1)]
    
    if (is.null(yy)) { 
        yy <- 0
    } else {
        yy <- model[, (ncol(X) + 2):ncol(model)]
    }
    n <- nrow(model)
    
    sx <- length(start_x)
    
    rhs_cv <- function(p) {
        h <- p[1]
        pr <- p[2:(1 + sx)]
        if (is.null(yy)) pyy <- 0 else {
            pyy <- p[(2 + sx):length(p)]
        }
        yar <- yy %*% pyy
        u <- y - yar
        yar + cv_np(u, as.vector(X %*% weight(pr, ncol(X))), h, degree)
    }
    
    fn0 <- function(p) {
        sum((y - rhs_cv(p))^2)
    }
    
    start <- c(start_bws, start_x, start_ar)
    opt <- optimx(start, fn0, method = method,...)
    bmet <- which.min(opt$value)
    par <- as.numeric(opt[bmet, 1:length(start)])   
    call <- match.call()
    
    rhs <- function(p) {
        h <- p[1]
        pr <- p[2:(1 + sx)]
        if (is.null(yy)) pyy <- 0 else {
            pyy <- p[(2 + sx):length(p)]
        }
        yar <- yy %*% pyy
        u <- y - yar
        xi <-  as.vector(X %*% weight(pr, ncol(X)))
        np <- g_np(u,xi,  xeval = xi, h, degree)
        list(fitted = yar + np, xi = xi, np = np)
    }
    fitted.values <- as.vector(y - rhs(par)$fitted)
    names(par) <- c("h", paste0("x", 1:length(start_x)), paste0("y", 1:length(start_ar)))

    list(coefficients = par,
         midas_coefficients = weight(par[2:(1 + sx)],ncol(X)),
         bws = par[1],
         model = model,
         weights = weight,
         fn0 = fn0,
         rhs_cv = rhs_cv,
         rhs = rhs,
         opt = opt,
         call = call,
         hessian = function(x)numDeriv::hessian(fn0,x),
         fitted.values = fitted.values,
         residuals = as.vector(y - fitted.values),
         start = start)
    
}

##' MIDAS Partialy linear non-parametric regression
##'
##' Function for fitting PL MIDAS regression without the formula interface
##' @param y model response
##' @param X prepared matrix of high frequency variable lags for MMM term
##' @param z a vector, data for the non-parametric part
##' @param p.ar length of AR part
##' @param weight the weight function
##' @param degree the degree of local polynomial
##' @param start_bws the starting values bandwith
##' @param start_x the starting values for weight function
##' @param start_ar the starting values for AR part. Should be the same length as \code{p}
##' @param method a method passed to \link{optimx}
##' @param ... additional parameters to \link{optimx}
##' @return an object similar to \code{midas_r} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' 
##' @import optimx
##' @importFrom stats na.omit 
##'
##' @export
##' 
midas_pl_plain <- function(y, X, z, p.ar = NULL, weight, degree = 1, start_bws, start_x, start_ar = NULL, method = c("Nelder-Mead"), ...) {
    d <- ncol(X)
    
    yy <- NULL 
    
    if(!is.null(p.ar)) {
        p.ar <- as.integer(p.ar)
        if (length(start_ar)!=p.ar) stop("The length of starting values vector for AR terms must the same as the number of AR terms")
        yy <- as.matrix(mls(y, 1:p.ar, 1))
    } 
    z <- as.numeric(z)
    y <- as.numeric(y)
    if ((length(z) != length(y)) &  (length(y) != nrow(X))) stop("The dimensions of the data do not match, z and y must be vectors of equal length, which must coincide with number of rows of X")
    
    model <- na.omit(cbind(y, X, z, yy ))
    
    y <- model[,1]
    X <- model[, 2:(ncol(X) + 1)]
    z <- model[, ncol(X) + 2]
   
    if (is.null(yy)) { 
        yy <- 0
    } else {
        yy <- model[, (ncol(X) + 3):ncol(model)]
    }
    n <- nrow(model)
    
    sx <- length(start_x)
    
    rhs_cv <- function(p) {
        h <- p[1]
        pr <- p[2:(1 + sx)]
        if (is.null(yy)) pyy <- 0 else {
            pyy <- p[(2 + sx):length(p)]
        }
        yar <- yy %*% pyy
        xi <- as.vector(X %*% weight(pr, ncol(X)))
        u <- y - yar - xi
        yar + xi + cv_np(u, z, h, degree)
    }
    
    fn0 <- function(p) {
        sum((y - rhs_cv(p))^2)
    }
    start <- c(start_bws, start_x, start_ar)
    opt <- optimx(start, fn0, method = method,...)
    bmet <- which.min(opt$value)
    par <- as.numeric(opt[bmet, 1:length(start)])   
    call <- match.call()
    
    rhs <- function(p) {
        h <- p[1]
        pr <- p[2:(1 + sx)]
        if (is.null(yy)) pyy <- 0 else {
            pyy <- p[(2 + sx):length(p)]
        }
        yar <- yy %*% pyy
        xi <-  as.vector(X %*% weight(pr, ncol(X)))
        u <- y - yar  - xi
        
        np <- g_np(u, z, xeval = z, h, degree)
        list(fitted = yar + xi + np, xi = xi, np = np, z = z)
    }
    fitted.values <- as.vector(y - rhs(par)$fitted)
    names(par) <- c("h", paste0("x", 1:length(start_x)), paste0("y", 1:length(start_ar)))
    
    list(coefficients = par,
         midas_coefficients = weight(par[2:(1 + sx)],ncol(X)),
         bws = par[1],
         model = model,
         weights = weight,
         fn0 = fn0,
         rhs_cv = rhs_cv,
         rhs = rhs,
         opt = opt,
         call = call,
         hessian = function(x)numDeriv::hessian(fn0,x),
         fitted.values = fitted.values,
         residuals = as.vector(y - fitted.values),
         start = start)
    
}


cv_np <- function(y, x, h, degree = 1) {
    cvg <- rep(NA, length(y))
    for (i in 1:length(y)) {
        if (is.null(dim(x))) w <- kfun(x[i], x, h)
        else  w <- kfun(x[i, ], x, h)
        if (degree > 0) {
            if (degree == 1) {
                if (is.null(dim(x))) { 
                    xc <- matrix(x - x[i], ncol = 1)
                } else xc <- cbind(1, sweep(x, 2, x[i, ], "-"))
            } else {
                xc <- poly(x - x[i], degree = degree, raw = TRUE, simple = TRUE)
            }
            cc <- lsfit(xc[-i,], y[-i], wt = w[-i], intercept = TRUE)
            cvg[i] <- coef(cc)[1]
        } else {
            cvg[i] <- sum(w[-i]*y[-i])/sum(w[-i])
        }
    }
    cvg
}

g_np <- function(y, x, xeval, h, degree = 1) {
    res <- rep(NA, length(xeval))
    for (i in 1:length(xeval)) {
        w <- kfun(xeval[i], x, h)
        if(degree > 0) {
            xc <- poly(x - xeval[i], degree = degree, raw = TRUE, simple = TRUE)
            cc <- lsfit(xc, y, wt = w, intercept = TRUE)
            res[i] <- coef(cc)[1]
        } else {
            res[i] <- sum(w*y)/sum(w)
        }
    }
    res
}

#' @importFrom stats dnorm
kfun <- function(z0, z, h) {
    if (is.null(ncol(z))) {
        dnorm((z - z0)/exp(h))
    } else {
        z0 <- as.numeric(z0)
        h <- as.numeric(h)
        apply(dnorm(scale(z, center = as.numeric(z0), scale = exp(h))), 1, prod)
    }
}

all_coef_full <- function(p, pinds, rf) {              
    pp <- lapply(pinds,function(x)p[x])     
    res <- mapply(function(fun,param)fun(param),rf ,pp,SIMPLIFY=FALSE)
    return(res)
}