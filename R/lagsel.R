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

##' @export
##' @method print midas_r_iclagtab
print.midas_r_iclagtab <- function(x,...) {
    print(x$table,...)
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
##' lagsel(mlr,"BIC","unrestricted")
##'
##' @details This function selects the model from the lag selection table for which the chosen information criteria achieves the smallest value
lagsel <- function(x,IC=x$IC[1],type=c("restricted","unrestricted")) {
    if(!(IC%in%x$IC))stop("The supplied information criteria was not used in creating lag selection table")
    type <- match.arg(type)
    coln <- paste(IC,type,sep=".")
    i <- which.min(x$table[,coln])
    cat("\n Selected model with lag structure: ", as.character(x$table$lags[i]))
    cat("\n ",IC," = ",x$table[i,coln],"")
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
