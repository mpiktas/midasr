##' Restricted MIDAS regression with I(1) regressors
##'
##' Estimate restricted MIDAS regression using non-linear least squares, when the regressor is I(1)
##'
##' @param x either formula for restricted MIDAS regression or \code{midas_r} object. Formula must include \code{\link{fmls}} function
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param model one of \code{"onestep"}, \code{"twosteps"} or \code{"reduced"}, see the details.
##' @param start the starting values for optimisation. Must be a list with named elements.
##' @param Ofunction the list with information which R function to use for optimisation. The list must have element named \code{Ofunction} which contains character string of chosen R function. Other elements of the list are thearguments passed to this function.  The default optimisation function is \code{\link{optim}} with argument \code{method="BFGS"}. Other supported functions are \code{\link{nls}}
##' @param gradient the default value is \code{NULL}, which means that the numeric approximation of weight function gradient is calculated. For any other value it is assumed that the R function for weight function gradient has the name of the weight function appended with \code{.gradient}. This function must return the matrix with dimensions \eqn{d_k \times q}, where \eqn{d_k} and \eqn{q} are the numbers of coefficients in unrestricted and restricted regressions correspondingly.
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
##' imr.t0 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="reduced",start=list(x=c(-0.1,10,-10,-10)))
##' 
##' imr.t2 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="twosteps",start=list(x=c(-0.1,10,-10,-10)))
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
imidas_r.default <- function(x, ldata=NULL, hdata=NULL, model=c("onestep","twosteps","twostep","reduced"), start, Ofunction="optim", gradient=NULL,...) {

    Zenv <- new.env(parent=environment(x))

    step1 <- vector("list",5)
    names(step1) <- c("weights","gradD","unrestricted","mcf","betad")

    model <- match.arg(model)
    
    
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
#        mfg[[3]] <- d
        r <- eval(mfg,Zenv)
        apply(r,2,cumsum)[1:d,]
    }

    step1$weights <- pp
    
    if(is.null(gradient)) {
        step1$gradD <- function(p,d)jacobian(pp,p,d=d)        
    }
    else {
        step1$gradD <- pp.gradient
    }
    if(model=="reduced") {
        diff <- -1
    }
    else {
        diff <- 0
    }
                          
    formula <- expandfmls(formula(x),"pp",Zenv,diff)    
    cl <- match.call(expand.dots=TRUE)
    cl <- cl[names(cl)!="model"]
 
    assign("pp",pp,Zenv)
    assign("pp.gradient",pp.gradient,Zenv)
    environment(formula) <- Zenv
    cl[[2]] <- formula
    if(model=="onestep") {
        cl[[1]] <- as.name("midas_r")
        res <- eval(cl,Zenv)
    }
    else {
        m <- match(c("x","ldata","hdata"),names(cl),0L)
      
        mf.mu <- cl[c(1,m)]        
        mf.mu[[1]] <- as.name("midas_u")
        names(mf.mu)[2] <- "formula"
        mu <- eval(mf.mu,Zenv)
        trform <- expandfmls(formula(x),"pp",Zenv,diff,truncate=TRUE)
        mttr <- terms(trform)
        cf <- setdiff(attr(mu$terms,"term.labels"),attr(mttr,"term.labels"))
        step1$unrestricted <- mu
        step1$mcf <- coef(mu)[names(coef(mu))!=cf]
        step1$betad <- coef(mu)[cf]
        
        u <- mu$model[,1]-mu$model[,cf]*coef(mu)[cf]        
        
        if(missing(ldata)|missing(hdata)) {
             ee <- NULL
         } else {
             data <- check_mixfreq(ldata,hdata)
                 
             ee <- as.environment(c(as.list(data$lowfreq),as.list(data$highfreq)))
             parent.env(ee) <- parent.frame()
        }
        
        assign("ee",ee,Zenv)
        cl <- match.call()
        mf <- match.call(expand.dots = FALSE)
        ##Fix this!!
        m <- match(c("x", "ldata"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf[[1L]] <- as.name("model.frame")
        mf[[2L]] <- trform
        mf[[3L]] <- as.name("ee")   
        mf[[4L]] <- as.name("na.omit")
        names(mf)[c(2,3,4)] <- c("formula","data","na.action")
        
        mf <- eval(mf,Zenv)
        mt <- attr(mf, "terms")
        args <- list(...)

        X <- model.matrix(mt,mf)

        prepmd <- prepmidas_r(u,X,mt,Zenv,cl,args,start,Ofunction,gradient)        
        class(prepmd) <- "midas_r"
        res <- midas_r.fit(prepmd)
    }
    res$imodel <- model
    res$step1 <- step1
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
##' Test restrictions on coefficients of MIDAS regression with I(1) regressors
##'
##' Perform a test whether the restriction on MIDAS regression coefficients with I(1) regressors holds.
##' @param x a MIDAS regression model with restricted coefficients, estimated with \code{\link{imidas_r}}
##' @return a \code{htest} object
##' @author Benediktas Bilinskas, Virmantas Kvedaras, Vaidotas Zemlys
##' @seealso hAh.test, hAhr.test
##' @examples
##' theta.h0 <- function(p, dk) {
##'   i <- (1:dk-1)/100
##'   pol <- p[3]*i + p[4]*i^2
##'  (p[1] + p[2]*i)*exp(pol)
##' }
##' theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)
##' 
##' xx <- simplearma.sim(list(ar=1),1500*12,1,12)
##' y <- midas.sim(500,theta0,xx,1)
##' x <- window(xx,start=start(y))
##'
##' imr <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,start=list(x=c(-0.1,10,-10,-10)))
##' 
##' imr.t0 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="reduced",start=list(x=c(-0.1,10,-10,-10)))
##' 
##' imr.t2 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="twosteps",start=list(x=c(-0.1,10,-10,-10)))
##'
##' ihAh.test(imr)
##' ihAh.test(imr.t0)
##' ihAh.test(imr.t2)
##' 
##' @details The test is chosen depending on the model. For \code{"onestep"} MIDAS regression test  \code{hAh.test} is used. For \code{"twosteps"} a modified version of the test \code{ihAh.nls.test} is used. For \code{"reduced"} special test \code{ihAh.Td.test} is calculated.
##' @export
ihAh.test <- function(x) {
    switch(x$imodel,
           onestep=hAh.test(x)
#           twosteps=ihAh.nls.test(x),
#           reduced=ihAh.Td.test(x)
           )
}

##' Test restrictions on coefficients of MIDAS regression with I(1) regressors
##'
##' Perform a test whether the restriction on MIDAS regression coefficients with I(1) regressors holds, in case where twostep estimation is used.
##' @param x a MIDAS regression model with restricted coefficients, estimated with \code{\link{imidas_r}} 
##' @return a \code{htest} object
##' @author Benediktas Bilinskas, Virmantas Kvedaras, Vaidotas Zemlys
##' @examples
##' theta.h0 <- function(p, dk) {
##'   i <- (1:dk-1)/100
##'   pol <- p[3]*i + p[4]*i^2
##'  (p[1] + p[2]*i)*exp(pol)
##' }
##' theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)
##' xx <- simplearma.sim(list(ar=1),1500*12,1,12)
##' y <- midas.sim(500,theta0,xx,1)
##' x <- window(xx,start=start(y))
##'
##' imr.t2 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="twosteps",start=list(x=c(-0.1,10,-10,-10)))
##'
##' ihAh.nls.test(imr.t2)
ihAh.nls.test <- function(x,se.type=c("ols","nls")) {
    se.type <- match.arg(se.type)
    X <- x$model[,-1]
    XtX <- crossprod(X)

    n_d <- nrow(X)
    dk <- ncol(XtX)

    if(se.type=="ols") {
        mu <- x$step1$unrestricted
        se2 <- sum(residuals(mu)^2)/(n_d-length(coef(mu)))
    }
    else {
        se2 <- sum(residuals(x)^2)/(n_d-length(coef(x)))
    }
        
    D0 <- x$gradD(coef(x))    
    Delta.0 <- D0%*%tcrossprod(ginv(1/n_d*crossprod(D0,XtX)%*%D0),D0)
    Sigma <- ginv(XtX/n_d)-Delta.0

    cfur <- x$step1$mcf
    h.0 <- sqrt(n_d)*(cfur-x$midas.coefficients)/sqrt(se2)

    STATISTIC <- t(h.0)%*%ginv(Sigma)%*%h.0

    names(STATISTIC) <- "ihAh"
    METHOD <- "ihAh restriction test"
    PARAMETER <- dk-length(coef(x))
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")        
}

##' Test restrictions on coefficients of MIDAS regression with I(1) regressors
##'
##' Perform a test whether the restriction on MIDAS regression coefficients with I(1) regressors holds, in case where twostep estimation is used.
##' @param x a MIDAS regression model with restricted coefficients, estimated with \code{\link{imidas_r}}
##' @return a \code{htest} object
##' @author Benediktas Bilinskas, Virmantas Kvedaras, Vaidotas Zemlys
##' @examples
##' theta.h0 <- function(p, dk) {
##'   i <- (1:dk-1)/100
##'   pol <- p[3]*i + p[4]*i^2
##'  (p[1] + p[2]*i)*exp(pol)
##' }
##' theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)
##' 
##' xx <- simplearma.sim(list(ar=1),1500*12,1,12)
##' y <- midas.sim(500,theta0,xx,1)
##' x <- window(xx,start=start(y))
##'
##' imr.t0 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="reduced",start=list(x=c(-0.1,10,-10,-10)))
##'
##' ihAh.Td.test(imr.t0)
ihAh.Td.test <- function(x,se.type=c("ols","nls")) {
    X <- x$model[,-1]
    XtX <- crossprod(X)

    n_d <- nrow(X)
    dk <- ncol(XtX)
    
    se.type <- match.arg(se.type)
    if(se.type=="ols") {
        mu <- x$step1$unrestricted
        se2 <- sum(residuals(mu)^2)/(n_d-length(coef(mu)))
    }
    else {
        se2 <- sum(residuals(x)^2)/(n_d-length(coef(x)))
    }

    D0.all <- x$step1$gradD(coef(x),dk+1)

    D0.start <- D0.all[1:dk,]
    d.end <- D0.all[dk+1,]
    
    Sigma <- t(d.end)%*%ginv(1/n_d*t(D0.start)%*%XtX%*%D0.start)%*%d.end
    print(Sigma)
    h <- sqrt(n_d)*(x$step1$betad-x$step1$weights(coef(x),dk+1)[dk+1])/sqrt(se2)
    print(h)
    print(se2)
    STATISTIC <- h^2/Sigma
    names(STATISTIC) <- "ihAh T_d"
    METHOD <- "ihAh T_d restriction test"
    PARAMETER <- 1
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")            
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

imidas_r_fast <- function(y,x,dk,weight,start,imodel=c("twosteps","twosteps2","reduced","reduced2"),model.matrix=NULL) {

    imodel <- match.arg(imodel)
    
    if(imodel=="twosteps"| imodel=="twosteps2") {
        wt <- function(p) {
            cumsum(weight(p,dk+1))
        }
        wt1 <- NULL
        gradD1 <- NULL
    }
    else {
        wt1 <- function(p,d) {
            cumsum(weight(p,dk+1))[1:d]
        }
        wt <- function(p) {
            wt1(p,dk)
        }
        gradD1 <- function(p,d)jacobian(wt1,p,d=d)
    }    
    
    if(is.null(model.matrix)) {
        dd <- switch(imodel,twosteps=dk,twosteps2=dk,reduced=dk-1,reduced2=dk-1)
        m <- frequency(x)        
        V <- dmls(x,dd,m)
        xmd <- mls(x,dd+1,m)    
        y <- as.numeric(y)
        model <- na.omit(cbind(y,xmd,V))
    }
    else {
        model <- model.matrix
    }

    mu <- lsfit(model[,-1],model[,1],intercept=FALSE)

    if(imodel=="reduced2" | imodel=="twosteps2") {
        mu.alt <- lsfit(model[,2],model[,1],intercept=FALSE)
        betad <- coef(mu.alt)[1]
        mcf <- NULL
        if(imodel=="twosteps2") {
            u.alt <- model[,1]-betad*model[,2]
            mu.mcf <- lsfit(model[,-2:-1],u.alt,intercept=FALSE)
            mcf <- coef(mu.mcf)            
        }
    }
    else {
        mcf <- coef(mu)[-1]
        betad <- coef(mu)[1]
    }
    u <- model[,1]-betad*model[,2]

    model <- model[,-2]
    fun <- function(p) {
        r <- wt(p)
        sum((u-model[,-1]%*%r)^2)
    }
    opt <- optim(start,fun,method="BFGS",control=list(maxit=1000,reltol=sqrt(.Machine$double.eps)/10),hessian=T)
    step1 <- list(unrestricted=mu,
                  mcf=mcf,
                  betad=betad,
                  weights=wt1,
                  gradD=gradD1)
    
    opt$gr0 <- grad(fun,opt$par)

    resid <- u-model[,-1]%*%wt(opt$par)
    
    list(coefficients=opt$par,
         midas.coefficients=wt(opt$par),
         gradD=function(p)jacobian(wt,p),
         residuals=resid,
         opt=opt,
         model=model,
         step1=step1,
         fn0=fun,
         imodel=imodel)
}

