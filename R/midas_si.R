##' MIDAS Single index regression
##'
##' Function for fitting MMM MIDAS regression without the formula interface
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
    
    z <- NULL 
    if(!is.null(p.ar)) {
        p.ar <- as.integer(p.ar)
        if (length(start_ar)!=p.ar) stop("The length of starting values vector for AR terms must the same as the number of AR terms")
        z <- as.matrix(mls(y, 1:p.ar, 1))
    } 
    
    model <- na.omit(cbind(y,X,z))
    
    y <- model[,1]
    X <- model[, 2:(ncol(X) + 1)]
    if (is.null(z)) { 
        z <- 0
    } else {
        z <- model[, (ncol(X) + 2):ncol(model)]
    }
    n <- nrow(model)
    
    sx <- length(start_x)
    
    rhs_cv <- function(p) {
        h <- p[1]
        pr <- p[2:(1 + sx)]
        if (is.null(z)) pz <- 0 else {
            pz <- p[(2 + sx):length(p)]
        }
        yar <- z %*% pz
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
        if (is.null(z)) pz <- 0 else {
            pz <- p[(2 + sx):length(p)]
        }
        yar <- z %*% pz
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

cv_np <- function(y, x, h, degree = 1) {
    cvg <- rep(NA, length(y))
    for(i in 1:length(y)) {
        w <- kfun(x[i], x, h)
        if(degree > 0) {
            xc <- poly(x - x[i], degree = degree, raw = TRUE, simple = TRUE)
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

kfun <- function(z0, z, h) {
    dnorm((z - z0)/exp(h))
}