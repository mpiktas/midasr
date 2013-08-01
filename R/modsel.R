##' Create a high frequency lag selection table for MIDAS regression model
##'
##' Creates a high frequency lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation
##' @param kmin the minimum high frequency lag, defaults to zero.
##' @param kmax the highest high frequency lag, defaults to square root of number of low fequency observations
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
##' 
##' mlr <- hf_lags_table(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=4,kmax=6)
##' mlr
##'
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} of the last term of the given formula
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
hf_lags_table<- function(formula,data,start,kmin=NULL,kmax=NULL,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepenv(data,Zenv,cl,mf,parent.frame())

    if(is.null(kmax))kmax <- round(10*log(length(prep$y),base=10))
    lti <- lastterminfo(formula,prep$Zenv)

    nop <- length(start[[lti$varname]])

    mkmin <- min(lti$lags)
    if(is.null(kmin)) {
        kmin <- mkmin+nop
    }
    else {
        kmin <- min(kmin,mkmin+nop)
    }

    if(kmax<kmin) {
        stop("Maximum number of lags must be larger than the sum of minimum lag and number of parameters")
    }
    lags <- lapply(kmin:kmax,function(x)mkmin:x)
    wstart <- start[lti$varname]
    names(wstart) <- lti$weight
    start <- start[names(start)!=lti$varname]
 
    make_ic_table(formula,data,start=start,wstart=wstart,table=list(lags,lti$weight),IC=IC,test=test,Ofunction=Ofunction,user.gradient=FALSE,...)
}

##' Create a low frequency lag selection table for MIDAS regression model
##'
##' Creates a low frequency lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation
##' @param kmin the minimum low frequency lag, defaults to zero.
##' @param kmax the highest low frequency lag, defaults to square root of number of low fequency observations
##' @param IC the information criterias which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param user.gradient see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a \code{midas_r_ic_table} object which is the list with the following elements:
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
##' mlr <- lf_lags_table(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=2,kmax=3)
##' mlr
##'
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} of the last term of the given formula
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
lf_lags_table<- function(formula,data,start,kmin=NULL,kmax=NULL,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepenv(data,Zenv,cl,mf,parent.frame())
    lti <- lastterminfo(formula,prep$Zenv)
    nop <- length(start[[lti$varname]])
    m <- lti$frequency

    if(is.null(kmax))kmax <- round(10*log(length(prep$y),base=10))%/%m
    if((kmax*m)>length(prep$y))stop("Maximum lag is too big, not enough data points")
       
    mkmin <- min(lti$lags)
    if(is.null(kmin)) {
        kmin <- max((mkmin+nop)%/%m+1,1)
    }
    else {
        kmin <- min(kmin,(mkmin+nop)%/%m+1)
    }

    if(kmax<kmin) {
        stop("Maximum number of lags must be larger than the sum of minimum lag and number of parameters")
    }
    lags <- lapply(kmin:kmax,function(x)(mkmin*m):(x*m-1))
    wstart <- start[lti$varname]
    names(wstart) <- lti$weight
    start <- start[names(start)!=lti$varname]
 
    make_ic_table(formula,data,start=start,wstart=wstart,table=list(lags,lti$weight),IC=IC,test=test,Ofunction=Ofunction,user.gradient=FALSE,...)
}


lastterminfo <- function(x,Zenv) {
    last.term <- x[[3]]
    if(length(last.term)==3) {
        last.term <- x[[c(3,3)]]
        ind <- c(3,3)
    }
    else {
        ind <- 3
    }
    lags <- as.numeric(eval(last.term[[3]],Zenv))
    freq <- as.numeric(eval(last.term[[4]],Zenv))
    weightname <- as.character(last.term[[5]])
    mtype <- as.character(last.term[[1]])
  
    if(!(mtype%in%c("fmls","dmls","mls")))stop("The last term in the formula must be a MIDAS lag term")
    if(mtype=="fmls")lags <- 0:lags

    list(lags=lags,weight=weightname,varname=as.character(last.term[[2]]),frequency=freq)  
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
##' mhfr <- hf_lags_table(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=4,kmax=6)
##' 
##' mlfr <- lf_lags_table(y~trend+fmls(x,12,12,nealmon),start=list(x=rep(0,3)),kmin=1,kmax=2)
##'
##' mwlr <- make_ic_table(y~trend+fmls(x,12,12,nealmon),table=list(list(0:4,0:5,0:6),c("nealmon","nealmon","nbeta")),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)))
##' 
##' modsel(mhfr,"BIC","unrestricted")
##'
##' modsel(mlfr,"BIC","unrestricted")
##'
##' modsel(mwlr,"BIC","unrestricted")
##' 
##' @details This function selects the model from the model selection table for which the chosen information criteria achieves the smallest value. The function works with model tables produced by functions \link{lf_lags_table}, \link{hf_lags_table} and \link{make_ic_table}.
modsel <- function(x,IC=x$IC[1],type=c("restricted","unrestricted")) {
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

##' Create a weight function selection table for MIDAS regression model
##'
##' Creates a weight function selection table for MIDAS regression model with given information criteria and weight functions.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation
##' @param weights the names of weight functions. Weight functions must have only 3 default parameters.
##' @param wstart the starting values for weight functions, a named list.
##' @param IC the information criterias which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param user.gradient see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a \code{midas_r_ic_table} object which is the list with the following elements:
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
##' mwr <- weights_table(y~trend+fmls(x,12,12,nealmon),weights=c("nealmon","nbeta"),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)))
##'
##' mwr
##'
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} of the last term of the given formula
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
weights_table <- function(formula,data,start=NULL,weights,wstart,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
    Zenv <- new.env(parent=environment(formula))
      
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepenv(data,Zenv,cl,mf,parent.frame())

    lti <- lastterminfo(formula,prep$Zenv)
    lags <- lti$lags
 
    make_ic_table(formula,data,start=start,wstart=wstart,table=list(list(lags),weights),IC=IC,test=test,Ofunction=Ofunction,user.gradient=FALSE,...)
}


##' Create a weight and lag selection table for MIDAS regression model
##'
##' Creates a weight and lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation excluding the starting values for the last term
##' @param table a list with first element a list of lag structures and the second element the vector of weight function names. The second element is recycled if its length does not coincide with the length of the first element.
##' @param wstart the starting values for weight functions, a named list.
##' @param IC the information criterias which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param user.gradient see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a \code{midas_r_ICtable} object which is the list with the following elements:
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
##' mwlr <- make_ic_table(y~trend+fmls(x,12,12,nealmon),table=list(list(0:4,0:5,0:6),c("nealmon","nealmon","nbeta")),wstart=list(nealmon=rep(0,3),nbeta=c(1,1,1,0)))
##'
##' mwlr
##'
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} and varying the weights of the last term of the given formula 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
make_ic_table <- function(formula,data=NULL,start=NULL,table,wstart,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
    Zenv <- new.env(parent=environment(formula))
    formula <- as.formula(formula)
    args <- list(...)
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- formula

    prep <- prepenv(data,Zenv,cl,mf,parent.frame())

    Zenv <- prep$Zenv
    mff <- prep$mf
    
    wlinfo <- formula_table(formula,Zenv,table,start,wstart)
    varname <- wlinfo$varname
    wlinfo$varname <- NULL

    ##Remove those formulas for which the number of parameters is less or equal than number of lags.
    
    cond <- mapply(
        function(lags,start){
            ifelse(length(lags)>length(start),TRUE,FALSE)
        },
        wlinfo$lags,
        lapply(wlinfo$starts,function(x)x[[varname]]),
        SIMPLIFY=TRUE)
    
    wlinfo <- lapply(wlinfo,function(x)x[cond])
    
    modellist <- mapply(function(f,st) {    
        mff[[2L]] <- f
        mmf <- eval(mff,Zenv)
        mmt <- attr(mmf, "terms")
        y <- model.response(mmf, "numeric")
        X <- model.matrix(mmt, mmf)
        list(mt=mmt,y=y,X=X,start=st)
    },wlinfo$formulas,wlinfo$starts,SIMPLIFY=FALSE)

    maxlag <- which.max(sapply(wlinfo$lags,max))
    
    lrn <- rownames(modellist[[maxlag]]$X)

    modellist <- lapply(modellist,function(mm) {
        rn <- rownames(mm$X)
        ind <- match(lrn,rn)
        mm$y <- mm$y[ind]
        mm$X <- mm$X[ind,]
        mm
    })
    
    mrm <- lapply(modellist,function(mm) {
        res <- prepmidas_r(mm$y,mm$X,mm$mt,Zenv,cl,args,mm$start,Ofunction,user.gradient,NULL)
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
    tab <- data.frame(weights=names(wlinfo$weights),
                      lags=sapply(wlinfo$lags,deparse),
                      tab)
    res <- list(table=tab,candlist=candlist,IC=IC,weights=weights)
    class(res) <- "midas_r_ic_table"
    res
}

##' @export
##' @method print midas_r_ic_table
print.midas_r_ic_table <- function(x,...) {
    print(x$table,...)
}

formula_table <- function(x,Zenv,table,start,wstart) {
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

    last.term[[1]] <- as.name("mls")
       
    if(length(table[[2]])<length(table[[1]])) {
        table[[2]] <- rep(table[[2]],length.out=length(table[[1]]))
    }else {
        if(length(table[[1]])<length(table[[2]])) {
            table[[1]] <- rep(table[[1]],length.out=length(table[[2]]))
        }
    }
    
    formulas <- vector("list",length(table[[1]]))
    starts <- formulas
    varname <- as.character(last.term[[2]])
    table[[2]] <- as.list(table[[2]])
    for(i in 1:length(formulas)) {
        res <- x
        lt <- last.term
        wght <- table[[c(2,i)]]
        if(is.character(wght)) {
            lt[[5]] <- as.name(wght)
            nmwght <- wght
            names(table[[2]])[i] <- nmwght
        }
        else {
            if(!is.function(wght))stop("Supply either function name or a function")
            nmwght <- names(table[[2]])[i]
            assign(nmwght,wght,environment(res))
            lt[[5]] <- as.name(nmwght)
        }
        lt[[3]] <- table[[c(1,i)]]
        res[[ind]] <- lt
        formulas[[i]] <- res
        wst <- list(wstart[[nmwght]])
        names(wst) <- varname
        starts[[i]] <- c(start,wst)
    }
    ##FIX!!!
    ##Adapt the code so it correctly assumes that wieghts is the list now
    list(formulas=formulas,
         lags=table[[1]],
         weights=table[[2]],
         starts=starts,
         varname=varname)
}

prepenv <- function(data,Zenv,cl,mf,pf) {
##Get the response of the model to get the number of observations
##Get the model.frame object, not evaluated!
##Prepare data if necessary    
    if(is.null(data)) {
        ee <- NULL
    }
    else {
        if(is.matrix(data)) data <- data.frame(data)
        if(is.data.frame(data)) {
            ee <- as.enviroment(as.list(data))
        }
        else {
            if(is.list(data)) {
                data <- mapply(function(x,nm){
                    if(is.null(dim(x))) {
                        x <- list(x)
                        names(x) <- nm
                        x
                    } else {
                        as.list(x)
                    }
                },data,names(data),SIMPLIFY=FALSE)
                names(data) <- NULL
                ee <- as.environment(do.call("c",data))
            } else {
                stop("Argument data must be a matrix, data.frame or a list")
            }
        }
        parent.env(ee) <- pf
    }
    assign("ee",ee,Zenv)
    
    ##Fix this!!
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf[[3L]] <- as.name("ee")   
    mf[[4L]] <- as.name("na.omit")
    names(mf)[c(2,3,4)] <- c("formula","data","na.action")

    
    mff <- mf
    
    mf <- eval(mf,Zenv)
    mt <- attr(mf, "terms")
            
    y <- model.response(mf, "numeric")
    list(Zenv=Zenv,y=y,mf=mff)
}

