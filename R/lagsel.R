##' Create a lag selection table for MIDAS regression model
##'
##' Creates a lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param x the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation
##' @param kmin the minimum high frequency lag, defaults to zero.
##' @param kmax the highest high frequency lag, defaults to square root of number of low fequency observations
##' @param IC the information criterias which to compute
##' @param Ofunction see \link{midasr}
##' @param user.gradient see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a \code{midas_r_iclagtab} object which is the list with the following elements:
##'
##' \item{table}{the table where each row contains calculated information criteria for both restricted and unrestricted MIDAS regression model with given lag structure}
##' \item{candlist}{the list containing fitted models}
##' \item{IC}{the argument IC}
##' @examples
##'
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##' 
##' mlr <- iclagtab(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=4,kmax=6)
##' mlr
##'
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} of the last term of the given formula
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
iclagtab <- function(x,ldata=NULL,hdata=NULL,start,kmin=NULL,kmax=NULL,IC=c("AIC","BIC"),Ofunction="optim",user.gradient=FALSE,...) {
    
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

    laginfo <- lagformula(x,Zenv,kmin=kmin,kmax=kmax)
    
    modellist <- lapply(laginfo$formulas, function(f) {    
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

    ICfun <- lapply(IC,function(ic)eval(as.name(ic)))
    ictab <- lapply(candlist,function(mm) {
        c(sapply(ICfun,function(ic)ic(mm)),
          sapply(ICfun,function(ic)ic(mm$unrestricted)))        
    })
    ictab <- do.call("rbind",ictab)

    colnames(ictab) <- c(paste(IC,"restricted",sep="."),paste(IC,"unrestricted",sep="."))
    ictab <- data.frame(lags=sapply(laginfo$lags,deparse),ictab)
    res <- list(table=ictab,candlist=candlist,IC=IC)
    class(res) <- "midas_r_iclagtab"
    res
}

##' Select the model based on given information criteria
##'
##' Selects the model with minimum of given information criteria and model type
##' @param x and output from iclagtab function
##' @param IC the name of information criteria to choose
##' @param type the type of MIDAS model, either restricted or unrestricted
##' @return the best model based on information criteria, \link{midas_r} object
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
##' @examples
##'
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##' 
##' mlr <- iclagtab(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=4,kmax=6)
##' 
##' mwr <- icwtab(y~trend+fmls(x,12,12,nealmon),weights=c("nealmon","nbeta"),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)))
##'
##' mwlr <- icwlagtab(y~trend+fmls(x,12,12,nealmon),weights=c("nealmon","nbeta"),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)),kmin=4,kmax=6)
##' 
##' icsel(mlr,"BIC","unrestricted")
##'
##' icsel(mwr,"BIC","unrestricted")
##'
##' icsel(mwlr,"BIC","unrestricted")
##' 
##' @details This function selects the model from the model selection table for which the chosen information criteria achieves the smallest value. The function works with model tables produced by functions \link{iclagtab}, \link{icwtab} and \link{icwlagtab}.
 
icsel <- function(x,IC=x$IC[1],type=c("restricted","unrestricted")) {
    if(!(IC%in%x$IC))stop("The supplied information criteria was not used in creating lag selection table")
    type <- match.arg(type)
    coln <- paste(IC,type,sep=".")
    i <- which.min(x$table[,coln])
    cat("\n Selected model with ")
    #if("lags"%in% names(x$table))cat("lag structure: ", as.character(x$table$lags[i]))
    #if("weights"%in% names(x$table))cat(" weight function: ", as.character(x$table$weights[i]))
    
    cat(IC," = ",x$table[i,coln],"")
    cat("\n Based on ",type, " MIDAS regression model\n")
    print(summary(x$candlist[[i]]))
    invisible(x)
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
    list(lags=lags,formulas=formulas)    
}

weightformula <- function(x,Zenv,weights) {
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

    formulas <- vector("list",length(weights))
    for(i in 1: length(formulas)) {
        res <- x
        lt <- last.term
        lt[[5]] <- as.name(weights[i])
        res[[ind]] <- lt
        formulas[[i]] <- res
    }
    list(formulas=formulas,varname=as.character(last.term[[2]]))
}
##' Create a lag selection table for MIDAS regression model
##'
##' Creates a lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param x the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation
##' @param weights the names of weight function
##' @param wstart the starting values for weight functions, a named list.
##' @param IC the information criterias which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param user.gradient see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a \code{midas_r_iclagtab} object which is the list with the following elements:
##'
##' \item{table}{the table where each row contains calculated information criteria for both restricted and unrestricted MIDAS regression model with given lag structure}
##' \item{candlist}{the list containing fitted models}
##' \item{IC}{the argument IC}
##' @examples
##'
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##' mwr <- icwtab(y~trend+fmls(x,12,12,nealmon),weights=c("nealmon","nbeta"),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)))
##'
##' mwr
##'
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} of the last term of the given formula
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
icwtab <- function(x,ldata=NULL,hdata=NULL,start=NULL,weights,wstart,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
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


    winfo <- weightformula(x,Zenv,weights)
    
    modellist <- lapply(winfo$formulas, function(f) {    
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
    
    mrm <- mapply(function(mm,wst) {
        st <- c(list(wst),start)
        names(st)[1] <- winfo$varname
        res <- prepmidas_r(mm$y,mm$X,mm$mt,Zenv,cl,args,st,Ofunction,user.gradient,NULL,FALSE)
        class(res) <- "midas_r"
        res
    },modellist,wstart,SIMPLIFY=FALSE)
    browser()
    candlist <- lapply(mrm,midas_r)
    
    ICfun <- lapply(IC,function(ic)eval(as.name(ic)))
    tfun <- lapply(test,function(ic)eval(as.name(ic)))
    tab <- lapply(candlist,function(mm) {
        c(sapply(ICfun,function(ic)ic(mm)),
          sapply(ICfun,function(ic)ic(mm$unrestricted)),
          sapply(tfun,function(tt)tt(mm)$p.value)
          )        
    })
    tab <- do.call("rbind",tab)

    colnames(tab) <- c(paste(IC,"restricted",sep="."),paste(IC,"unrestricted",sep="."),paste(test,"p.value",sep="."))
    tab <- data.frame(weights=weights,tab)
    res <- list(table=tab,candlist=candlist,IC=IC,weights=weights)
    class(res) <- "midas_r_icwtab"
    res
}

##' Create a weight and lag selection table for MIDAS regression model
##'
##' Creates a weight and lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param x the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param ldata low frequency data, a \code{data.frame} object
##' @param hdata high frequency data, a \code{data.frame} object
##' @param start the starting values for optimisation excluding the starting values for the last term
##' @param weights the names of weight function
##' @param wstart the starting values for weight functions, a named list.
##' @param kmin the minimum high frequency lag, defaults to zero.
##' @param kmax the highest high frequency lag, defaults to square root of number of low fequency observations
##' @param IC the information criterias which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param user.gradient see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a \code{midas_r_icwlagtab} object which is the list with the following elements:
##'
##' \item{table}{the table where each row contains calculated information criteria for both restricted and unrestricted MIDAS regression model with given lag structure}
##' \item{candlist}{the list containing fitted models}
##' \item{IC}{the argument IC}
##' @examples
##'
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##' 
##'  
##' mwlr <- icwlagtab(y~trend+fmls(x,12,12,nealmon),weights=c("nealmon","nbeta"),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)),kmin=4,kmax=6)
##'
##' mwlr
##'
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} and varying the weights of the last term of the given formula 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
icwlagtab <- function(x,ldata=NULL,hdata=NULL,start=NULL,weights,wstart,kmin,kmax,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
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


    laginfo <- lagformula(x,Zenv,kmin=kmin,kmax=kmax)
    winfo <- weightformula(x,Zenv,weights)
    
    starts <- lapply(wstart,function(wst) {
        st <- c(list(wst),start)
        names(st)[1] <- winfo$varname
        st
    })
    
    
    wlinfo <- mapply(function(f,ll) {
        wf <- weightformula(f,Zenv,weights)
        res <- mapply(list,wf$formulas,starts,rep(list(ll),length(starts)),as.list(weights),SIMPLIFY=FALSE)
        res <- lapply(res,function(x){
            names(x) <- c("formula","start","lags","weight")
            x
        })
        res
    },laginfo$formulas,laginfo$lags,SIMPLIFY=FALSE)

    wlinfo <- do.call("c",wlinfo)

    modellist <- lapply(wlinfo, function(f) {    
        mff[[2L]] <- f$formula
        mmf <- eval(mff,Zenv)
        mmt <- attr(mmf, "terms")
        y <- model.response(mmf, "numeric")
        X <- model.matrix(mmt, mmf)
        list(mt=mmt,y=y,X=X,start=f$start)
    })

    maxlag <- which.max(sapply(wlinfo,with,max(lags)))
    
    lrn <- rownames(modellist[[maxlag]]$X)

    modellist <- lapply(modellist,function(mm) {
        rn <- rownames(mm$X)
        ind <- match(lrn,rn)
        mm$y <- mm$y[ind]
        mm$X <- mm$X[ind,]
        mm
    })
    
    mrm <- lapply(modellist,function(mm) {
        res <- prepmidas_r(mm$y,mm$X,mm$mt,Zenv,cl,args,mm$start,Ofunction,user.gradient,NULL,FALSE)
        class(res) <- "midas_r"
        res
    })
    
    candlist <- lapply(mrm,midas_r)
    
    ICfun <- lapply(IC,function(ic)eval(as.name(ic)))
    tfun <- lapply(test,function(ic)eval(as.name(ic)))
    tab <- lapply(candlist,function(mm) {
        c(sapply(ICfun,function(ic)ic(mm)),
          sapply(ICfun,function(ic)ic(mm$unrestricted)),
          sapply(tfun,function(tt)tt(mm)$p.value)
          )        
    })
    tab <- do.call("rbind",tab)

    colnames(tab) <- c(paste(IC,"restricted",sep="."),paste(IC,"unrestricted",sep="."),paste(test,"p.value",sep="."))
    tab <- data.frame(weights=sapply(wlinfo,with,weight),
                      lags=sapply(wlinfo,with,deparse(lags)),
                      tab)
    res <- list(table=tab,candlist=candlist,IC=IC,weights=weights)
    class(res) <- "midas_r_icwlagtab"
    res
}


##' @export
##' @method print midas_r_iclagtab
print.midas_r_iclagtab <- function(x,...) {
    print(x$table,...)
}

##' @export
##' @method print midas_r_icwtab
print.midas_r_icwtab <- function(x,...) {
    print(x$table,...)
}

##' @export
##' @method print midas_r_icwlagtab
print.midas_r_icwlagtab <- function(x,...) {
    print(x$table,...)
}
