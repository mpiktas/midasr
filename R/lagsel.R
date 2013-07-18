##' Select lags of MIDAS regression using information criteria
##'
##' Select lags of MIDAS regression using information criteria
##' @param x the formula for MIDAS regression, the lag selection is performed for the first mls term in the formula
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation
##' @param kmin the minimum high frequency lag, defaults to zero.
##' @param kmax the highest high frequency lag, defaults to square root of number of low fequency observations
##' @param IC the information criterias on which to base selection
##' @param Ofunction see \link{midasr}
##' @param user.gradient see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a table with 
##' @author Vaidotas Zemlys
lagsel <- function(x,ldata=NULL,hdata=NULL,start,kmin=NULL,kmax=NULL,IC=list(AIC,BIC),Ofunction="optim",user.gradient=FALSE,...) {
    
    Zenv <- new.env(parent=environment(x))
      
    if(missing(ldata)|missing(hdata)) {
        ee <- NULL
    }
    else {
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
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- as.name("na.omit")
    names(mf)[c(2,3,4)] <- c("formula","data","na.action")

    mff <- mf
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
        
    args <- list(...)
    y <- model.response(mf, "numeric")

    if(is.null(kmax))kmax <- sqrt(length(y))

    formulas <- lagformula(x,Zenv,kmin=kmin,kmax=kmax)
    
    modellist <- lapply(formulas, function(f) {    
        mff[[2L]] <- f
        mmf <- eval(mff,Zenv)
        mmt <- attr(mmf, "terms")
        y <- model.response(mmf, "numeric")
        X <- model.matrix(mmt, mmf)
        list(mt=mmt,y=y,X=X)
    })

    lrn <- rownames(modellist[[length(modellist)]]$X)
    modellist <- lapply(modellist,function(mm) {
        rn <- rownames(mm$X)
        ind <- match(lrn,rn)
        mm$y <- mm$y[ind]
        mm$X <- mm$X[ind,]
        mm
    })
    mrm <- lapply(modellist,function(mm) {
        res <- prepmidas_r(mm$y,mm$X,mm$mt,Zenv,cl,args,start,Ofunction,user.gradient,NULL,FALSE)
        class(res) <- "midas_r"
        res
    })
    candlist <- lapply(mrm,midas_r)

    candlist
}

lagformula <- function(x,Zenv,kmin=NULL,kmax) {
    last.term <- x[[3]]
    if(length(last.term)==3) {
        last.term <- x[[c(3,3)]]
        ind <- c(3,3)
    }
    else {
        ind <- 3
    }

    mtype <- as.character(last.term[[1]])
    if(!(mtype%in%c("fmls","dmls","mls")))stop("The last term in the formula must be a MIDAS lag term")

    if(mtype=="mls") {
        mkmin <- min(as.numeric(eval(last.term[[3]],Zenv)))
        if(is.null(kmin)) {
            kmin <- min(mkmin)
        }
    }
    else {
        if(is.null(kmin)) kmin <- 0
        mkmin <- 0
    }
    last.term[[1]] <- as.name("mls")

    lags <- lapply(kmin:kmax,function(x)mkmin:x)

    formulas <- vector("list",length(lags))
    for(i in 1: length(formulas)) {
        res <- x
        lt <- last.term
        lt[[3]] <- lags[[i]]
        res[[ind]] <- lt
        formulas[[i]] <- res
    }
    formulas    
}
