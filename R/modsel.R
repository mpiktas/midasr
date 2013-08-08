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
hf_lags_table<- function(formula,data,start,lagsinfo,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())

#    if(is.null(kmax))kmax <- round(10*log(length(prep$y),base=10))

    lti <- sapply(names(lagsinfo),function(nm)term_info(prep$mt,nm,prep$Zenv))

    nop <-sapply(start[sapply(lti,with,varname)],length)
    m <- sapply(lti,with,frequency)
    
    mkmin <- sapply(lapply(lti,with,lags),min)
    ###Fix to the end
    if(is.null(kmin)) {
        kmin <- mkmin+nop
    }
    else {
        kmin <- min(kmin,mkmin+nop)
    }

    if(kmax<kmin) {
        stop("The maximum number of lags is not large enough, no degrees of freedom for the weight function")
    }

    wstart <- start[lti$varname]
    names(wstart) <- lti$weight
    start <- start[names(start)!=lti$varname]
    
    table <- expand_weights_lags(lti$weight,kmin,kmax,1,mkmin,start=wstart)
    
    midas_r_ic_table(formula,data,start=start,table=table,IC=IC,test=test,Ofunction=Ofunction,user.gradient=FALSE,...)
}

##' Create a low frequency lag selection table for MIDAS regression model
##'
##' Creates a low frequency lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation
##' @param kmin the minimum low frequency lag, defaults to 1.
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
lf_lags_table <- function(formula,data,start,kmin=NULL,kmax=NULL,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())
    lti <- last_term_info(formula,prep$Zenv)  
    m <- lti$frequency

    if(is.null(kmax))kmax <- round(10*log(length(prep$y),base=10))%/%m
       
    mkmin <- min(lti$lags)
    nop <- length(start[[lti$varname]])

    ##Reuse high frequency code
    if(is.null(kmin)) {
        kmin <- mkmin+nop
    }
    else {
        kmin <- min(kmin*m,mkmin+nop)
    }
    kmin <- min(kmin%/%m,1)
    
    if(kmax<kmin) {
        stop("Maximum number of lags must be larger than the sum of minimum lag and number of parameters")
    }
       
    wstart <- start[lti$varname]
    names(wstart) <- lti$weight
    start <- start[names(start)!=lti$varname]
    
    table <- expand_weights_lags(lti$weight,kmin,kmax,m,mkmin,start=wstart)
        
    midas_r_ic_table(formula,data,start=start,table=table,IC=IC,test=test,Ofunction=Ofunction,user.gradient=FALSE,...)
}

last_term_info <- function(x,Zenv) {
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

term_info <- function(mt,term.name,Zenv) {
    vars <- as.list(attr(mt,"variables"))[-1]
    term.no <- find_mls_terms(varname,vars)
    
    last.term <- vars[[term.no]]
    
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
##' @param IC the name of information criteria to base the choosing of the model
##' @param test the name of the test for which to print out the p-value
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
##' modsel(mhfr,"BIC","unrestricted")
##'
##' modsel(mlfr,"BIC","unrestricted")
##' 
##' @details This function selects the model from the model selection table for which the chosen information criteria achieves the smallest value. The function works with model tables produced by functions \link{lf_lags_table}, \link{hf_lags_table}, \link{ghysels_table} and \link{midas_r_ic_table}.
modsel <- function(x,IC=x$IC[1],test=x$test[1],type=c("restricted","unrestricted")) {
    if(!(IC%in%x$IC))stop("The supplied information criteria was not used in creating lag selection table")
    type <- match.arg(type)
    coln <- paste(IC,type,sep=".")
    i <- which.min(x$table[,coln])
    cat("\n Selected model with ")
    #if("lags"%in% names(x$table))cat("lag structure: ", as.character(x$table$lags[i]))
    #if("weights"%in% names(x$table))cat(" weight function: ", as.character(x$table$weights[i]))
    
    cat(IC," = ",x$table[i,coln],"")
    cat("\n Based on",type, "MIDAS regression model\n")
    cat(" The p-value for the null hypothesis of the test", test, "is", x$table[i,paste(test,"p.value",sep=".")],"\n")
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
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())

    lti <- last_term_info(formula,prep$Zenv)
    lags <- lti$lags

    table <- expand_weights_lags(weights,1,1,1,0,start=wstart)
    table$lags <- rep(list(lags),length.out=length(table$weights))
    
    midas_r_ic_table(formula,data,start=start,table=table,IC=IC,test=test,Ofunction=Ofunction,user.gradient=FALSE,...)
}


##' Create a weight and lag selection table for MIDAS regression model
##'
##' Creates a weight and lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation excluding the starting values for the last term
##' @param table an wls_table object, see \link{expand_weights_lags}
##' @param IC the names of information criterias which to compute
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
##'  
##' mwlr <- midas_r_ic_table(y~trend+fmls(x,12,12,nealmon), table=list(weights=as.list(c("nealmon","nealmon","nbeta")), lags=list(0:4,0:5,0:6),starts=list(rep(0,3),rep(0,3,),c(1,1,1,0))))
##'
##' mwlr
##'
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} and varying the weights of the last term of the given formula 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname midas_r_ic_table
##' @export
midas_r_ic_table <- function(formula,...) UseMethod("midas_r_ic_table")

is.midas_r_ic_table <- function(x) inherits(x,"midas_r_ic_table")

#' @rdname midas_r_ic_table
#' @method midas_r_ic_table default
#' @export
midas_r_ic_table.default <- function(formula,data=NULL,start=NULL,table,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    
    Zenv <- new.env(parent=environment(formula))
    formula <- as.formula(formula)
    args <- list(...)
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- formula

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())

    Zenv <- prep$Zenv
    mff <- prep$mf

    ##Remove those formulas for which the number of parameters is less or equal than number of lags.
    remove_incomplete <- function(info,nm) {
        cond <- mapply(
            function(lags,start){
                ifelse(length(lags)>length(start),TRUE,FALSE)
            },
            info$lags,
            lapply(info$starts,function(x)x[[nm]]),
            SIMPLIFY=TRUE)    
        lapply(info,function(x)x[cond])
    }
    
    wlinfo <- remove_incomplete(formula_table(prep$mt,names(table)[1],Zenv,table[[1]],start),names(table)[1])

    combine <- function(l) {
        nms <- names(l[[1]])
        out <- lapply(nms,function(nm)do.call("c",lapply(l,function(x)x[[nm]])))
        names(out) <- nms
        out
    }
    
    if(length(table)>1) {
        for(i in 2:length(table)) {
            res <- mapply(function(f,s) {
                formula_table(terms(f),names(table)[i],Zenv,table[[i]],s)
            },wlinfo$formulas,wlinfo$starts,SIMPLIFY=FALSE)
            wlinfo <- combine(res)
            wlinfo <- remove_incomplete(wlinfo,names(table)[i])
        }        
    }
    
    
    
    modellist <- mapply(function(f,st) {    
        mff[[2L]] <- f
        if(!is.null(prep$itr$lagsTable)) {
            itr <-  checkARstar(terms(eval(mff[[2]], Zenv)))
            mff[[2]] <- itr$x
        } else itr <- NULL
        mmf <- eval(mff,Zenv)
        mmt <- attr(mmf, "terms")
        y <- model.response(mmf, "numeric")
        X <- model.matrix(mmt, mmf)
       # st <- st[!sapply(st,is.null)] #this is AR*, not working currently
        list(mt=mmt,y=y,X=X,start=st,itr=itr)
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
        res <- prepmidas_r(mm$y,mm$X,mm$mt,Zenv,cl,args,mm$start,Ofunction,user.gradient,mm$itr$lagsTable)
        class(res) <- "midas_r"
        res
    })
        
    candlist <- lapply(mrm,midas_r)

    make_ic_table(candlist,IC,test)
}

##' @method midas_r_ic_table midas_r_ic_table
##' @export
midas_r_ic_table.midas_r_ic_table <- function(formula,...) {
    do.call("make_ic_table",formula[-1])
}

make_ic_table <- function(candlist,IC,test) {

  #  if(is.null(weights)|is.null(lags)) {
  #      wl <- get_wl_from_cl(candlist)
  #      if(is.null(weights)) weights <- wl$weights
  #      if(is.null(lags)) lags <- wl$lags
  #  }
    
    ICfun <- lapply(IC,function(ic)eval(as.name(ic)))
    tfun <- lapply(test,function(ic)eval(as.name(ic)))

    tab <- lapply(candlist,function(mm) {
        c(sapply(ICfun,function(ic)ic(mm)),
          sapply(ICfun,function(ic)ifelse(is.null(mm$unrestricted),NA,ic(mm$unrestricted))),
          sapply(tfun,function(tt){
              tst <- try(tt(mm))
              ifelse(class(tst)=="try-error",NA,
              tst$p.value)
          })
          )        
    })
    tab <- do.call("rbind",tab)

    colnames(tab) <- c(paste(IC,"restricted",sep="."),paste(IC,"unrestricted",sep="."),paste(test,"p.value",sep="."))
    tab <- data.frame(model=sapply(candlist,with,deparse(terms)),
                      tab)
    res <- list(table=tab,candlist=candlist,IC=IC,test=test,weights=tab[,1],lags=tab[,2])
    class(res) <- "midas_r_ic_table"
    res
}

get_wl_from_cl <- function(candlist) {
    info <- lapply(candlist,function(l) {
        last_term_info(formula(l),l$Zenv)
    })   
    list(weights=sapply(info,with,weight),lags=sapply(info,with,deparse(lags)))
}


##' @export
##' @method print midas_r_ic_table
print.midas_r_ic_table <- function(x,...) {
    print(x$table,...)
}

formula_table <- function(mt,varname,Zenv,table,start) {
    if(is.null(names(table)))names(table) <- c("weights","lags","names")

    vars <- as.list(attr(mt,"variables"))[-1]
    term.no <- find_mls_terms(varname,vars)
    
#    last.term <- x[[3]]
#    if(length(last.term)==3) {
#        last.term <- x[[c(3,3)]]
#        ind <- c(3,3)
#    }
#    else {
#        ind <- 3
#    }
    
#    mtype <- as.character(last.term[[1]])
#    if(!(mtype%in%c("fmls","dmls","mls")))stop("The last term in the formula must be a MIDAS lag term")
    vars[[term.no]][[1]] <- as.name("mls")
#    last.term[[1]] <- as.name("mls")

    
    formulas <- vector("list",length(table$lags))
    starts <- formulas
    for(i in 1:length(formulas)) {
        res <- vars
        lt <- res[[term.no]]
        wght <- table$weights[[i]]
        if(is.character(wght)) {
            lt[[5]] <- as.name(wght)
        }
        else {
            if(!is.function(wght))stop("Supply either function name or a function")
            nmwght <- names(table$weights)[i]
            assign(nmwght,wght,environment(res))
            lt[[5]] <- as.name(nmwght)
        }
        lt[[3]] <- table$lags[[i]]
        res[[term.no]] <- lt
        formulas[[i]] <- variables_to_formula(res)
        wst <- list(table$starts[[i]])
        names(wst) <- varname
        starts[[i]] <- c(start,wst)
    }
    list(formulas=formulas,
         lags=table$lags,
         weights=table$weights,
         starts=starts)
}

prepare_model_frame <- function(data,Zenv,cl,mf,pf) {
##Get the response of the model to get the number of observations
##Get the model.frame object, not evaluated!
##Prepare data if necessary    
    if(is.null(data)) {
        ee <- NULL
    }
    else {
        if(is.matrix(data)) data <- data.frame(data)
        if(is.data.frame(data)) {
            ee <- as.environment(as.list(data))
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

    itr <- checkARstar(terms(eval(mf[[2]], Zenv)))
    
    mff <- mf
    mtf <- eval(mf,Zenv)
    mtt <- attr(mtf,"terms")

    #We need only response to get the number of low frequency observations.
    resf <- mf$formula
    resf[[3]] <- 1
    mf$formula <- resf

    mf <- eval(mf,Zenv)
    
    y <- model.response(mf, "numeric")
    list(Zenv=Zenv,y=y,mf=mff,itr=itr,mt=mtt)
}

##' Creates table of weights, lags and starting values 
##'
##' For each weight function creates lags starting from \code{kmin} to \code{kmax}. This is a convenience function for easier work with the function \link{midas_r_ic_table}. 
##' @title Create table of weights, lags and starting values
##' @param weights either a vector with names of the weight functions or a named list of weight functions
##' @param from the high frequency lags from which to start the fitting
##' @param to a vector of length two, containing minimum and maxmimum lags, high frequency if \code{m=1}, low frequency otherwise.
##' @param m the frequency ratio
##' @param start a named list with the starting values for weight functions
##' @return a \code{lws_table} object, a list with elements \code{weights}, \code{lags} and \code{starts}.
##' @examples
##'
##' expand_weights_lags(c("nealmon","nbeta"),0,c(4,8),1,start=list(nealmon=rep(0,3),nbeta=rep(0,4)))
##' nlmn <- expand_weights_lags("nealmon",0,c(4,8),1,start=list(nealmon=rep(0,3)))
##' nbt <- expand_weights_lags("nbeta",0,c(4,8),1,start=list(nbeta=rep(0,4)))
##'
##' nlmn+nbt
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
##' 
expand_weights_lags <- function(weights,from=0,to,m=1,start) {
    weights <- as.list(weights)
    kmin <- min(to)
    kmax <- max(to)
    mkmin <- from
    chnm <- sapply(weights,is.character)
    names(weights)[chnm] <- unlist(weights[chnm])
    if(!identical(names(weights),names(start)))stop("Mismatch between the  weight function names and the names of starting values")
    
    if(m>1) {
        lags <- lapply(kmin:kmax,function(x)(mkmin):(x*m-1))
    }
    else {
        lags <- lapply(kmin:kmax,function(x)mkmin:x)
    }
    
    weights <- rep(weights,each=length(lags))
    starts <- rep(start,each=length(lags))
    lags <- rep(lags,length.out=length(weights))
    out <- list(weights=weights,lags=lags,starts=starts)
    class(out) <- "lws_table"
    out
}

## @export
is.lws_table <- function(x) {
    length(unique(sapply(x,length)))==1
}

## @export
## @method print lws_table
print.lws_table <- function(x) {
    if(is.null(names(x)))names(x) <- c("weights","lags","starts")
    print(data.frame(weights=names(x$weights),lags=sapply(x$lags,deparse),starts=sapply(x$starts,deparse)))
}

##' Create table of weights, lags and starting values for Ghysels weight schema, see \link{ghyselslag}
##'
##' Given weight function creates lags starting from \code{kmin} to \code{kmax} and replicates starting values for each low frequency lag.
##' @title Create table of weights, lags and starting values for Ghysels weight schema
##' @param weight the names of weight functions
##' @param type the type of Ghysels schema, \code{"A"}, \code{"B"} or \code{"C"}
##' @param kmin the minimum low frequency lag
##' @param kmax the maximum high frequency lag
##' @param m the frequency ratio
##' @param start the starting values for the weights of the one low frequency lag
##' @param mkmin the minimum low frequency lag from which to start the fitting
##' @return a \code{lws_table} object, a list with elements \code{weights}, \code{lags} and \code{starts}
##' @examples
##' expand_ghysels("nealmon","A",1,2,12,c(0,0,0))
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
expand_ghysels <- function(weight,type=c("A","B","C"),kmin,kmax,m,start,mkmin=0) {
    mkmin <- mkmin%/%m
    lags <- lapply(kmin:kmax,function(x)(mkmin*m):(x*m-1))
    d <- sapply(lags,length)
    nm <- paste(weight,type,d,sep="_")
    type <- match.arg(type)
    if(type=="A") {
        starts <- lapply(d%/%m,function(lf)rep(start,times=lf))
    }
    if(type=="B") {
        starts <- lapply(d%/%m,function(lf)c(rep(start[1],lf),start[-1]))
    }
    if(type=="C") {
        starts <- lapply(d%/%m,function(lf)c(start[1],rep(start[-1],times=lf)))
    }    
    names(starts) <- nm
    ff <- expression(ghyselslag(p,d,m,weight,type))
    ff[[c(1,4)]] <- m
    ff[[c(1,5)]] <- as.name(weight)
    ff[[c(1,6)]] <- as.character(type)
    weights <- mapply(function(fun,e,dk){
        e[[c(1,3)]] <- as.numeric(dk)
        fun[[c(1,3,3,2)]] <- e[[1]]
        eval(fun)
    },rep(list(expression(f<-function(p,d,m){p})),length(d)),rep(list(ff),length(d)),as.list(d),SIMPLIFY=FALSE)
    names(weights) <- nm
    out <- list(weights=weights,lags=lags,starts=starts)
    class(out) <- "lws_table"
    out
}

##' Combines \code{lws_table} objects
##'
##' The \code{lws_table} objects have similar structure to table, i.e. it is a list with 3 elements which are the lists with the same number of elements. The base function \code{c} would \code{cbind} such tables. This function \code{rbind}s them.
##'  
##' @title Combine \code{lws_table} objects
##' @param ... \code{lws_table} object
##' @param check logical, if TRUE checks that the each \code{lws_table} object is named a list with names \code{c("weights","lags","starts")}
##' @return \code{lws_table} object
##' @rdname lws_table-add
##' @method + lws_table
##' @examples
##' nlmn <- expand_weights_lags("nealmon",4,8,1,0,start=list(nealmon=rep(0,3)))
##' nbt <- expand_weights_lags("nbeta",4,8,1,0,start=list(nbeta=rep(0,4)))
##'
##' nlmn+nbt
##' @export
##' @author Virmantas Kvedaras, Vaidotas Zemlys
"+.lws_table" <- function(...,check=TRUE)  {

    l <- list(...)
    if(check) {
        nms <- c("weights","lags","starts")
        for(i in 1:length(l)) {
            if(is.null(names(l[[i]])))names(l[[i]]) <- nms
        }
    }
    out <- lapply(nms,function(nm)do.call("c",lapply(l,function(x)x[[nm]])))
    names(out) <- nms
    class(out) <- "lws_table"
    out
}
##' Create weight and lag selection table for the MIDAS regression model with Ghysels schema
##'
##' This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} and varying the weights of the last term of the given formula 
##' @title Weight and lag selection table for MIDAS regression model with Ghysels schema
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param weights the names of weights used in Ghysels schema
##' @param wstart the starting values for the weights of the firs low frequency lag
##' @param type the type of Ghysels schema see \link{ghyselslag}, can be a vector of types
##' @param start the starting values for optimisation excluding the starting values for the last term
##' @param kmin the minimum low frequency lag, defaults to 1.
##' @param kmax the highest low frequency lag, defaults to square root of number of low fequency observations
##' @param IC the names of information criteria which should be calculated
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction  see \link{midasr}
##' @param user.gradient see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a \code{midas_r_ic_table} object which is the list with the following elements:
##'
##' \item{table}{the table where each row contains calculated information criteria for both restricted and unrestricted MIDAS regression model with given lag structure}
##' \item{candlist}{the list containing fitted models}
##' \item{IC}{the argument IC}
##' \item{test}{the argument test}
##' \item{weights}{the names of weight functions}
##' \item{lags}{the lags used in models}
##' ##' @examples
##'
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##' 
##' tb <- ghysels_table(y~trend+fmls(x,12,12,nealmon),data=list(y=y,x=x,trend=trend),weights=c("nealmon"),wstart=list(nealmon=c(0,0,0)),start=list(trend=1),type=c("A"),kmax=4)
##'
##' 
##' @details This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} and varying the weights of the last term of the given formula 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
ghysels_table <- function(formula,data,weights,wstart,type,start=NULL,kmin=NULL,kmax=NULL,IC=c("AIC","BIC"),test=c("hAh.test"),Ofunction="optim",user.gradient=FALSE,...) {
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())

    lti <- last_term_info(formula,prep$Zenv)
    lags <- lti$lags
 
    m <- lti$frequency

    if(is.null(kmax))kmax <- round(10*log(length(prep$y),base=10))%/%m
    
    mkmin <- min(lti$lags)
    nop <- length(start[[lti$varname]])
    ##Reuse high frequency code
    if(is.null(kmin)) {
        kmin <- mkmin+nop
    }
    else {
        kmin <- min(kmin*m,mkmin+nop)
    }
    kmin <- max(kmin%/%m,1)
    
    if(kmax<kmin) {
        stop("Maximum number of lags must be larger than the sum of minimum lag and number of parameters")
    }
    
    if(!is.list(weights)) {
        weights <- as.list(weights)
        ###Ugly code fix this
        names(weights) <- sapply(weights,function(x)x[[1]])
    }
    if(!is.list(wstart)) {
        wstart <- list(wstart)
        names(wstart) <- names(weights)
    }
    
    
    if(!identical(names(weights),names(wstart))) stop("Mismatch between the  weight function names and the names of starting values")

    if(length(weights)<length(type)) {
        weights <- rep(weights,length.out=length(type))
        wstart <- rep(wstart,length.out=length(type))
    }else {
        if(length(type)<length(weights)) {
            type <- rep(type,length.out=length(weights))
       }
    }
    
    tb <- mapply(function(w,t,s){
        expand_ghysels(weight=w,type=t,kmin=kmin,kmax=kmax,m=m,start=s,mkmin=mkmin)
    },as.list(weights),as.list(type),as.list(wstart),SIMPLIFY=FALSE)
    names(tb) <- NULL
    
    table <- do.call("+",tb)
    midas_r_ic_table(formula,data,start=start,table=table,IC=IC,test=test,Ofunction=Ofunction,user.gradient=FALSE,...)
}
    

add_expressions <- function(l) {
    if(length(l)<2) stop("You need 2 elements for addition")
    base <- expression(a+b)
    base[[c(1,2)]] <- l[[1]]
    base[[c(1,3)]] <- l[[2]]
    if(length(l)>2) {
        l <- l[-2:-1]
        for(i in 1:length(l)) {
            tmp <- expression(a+b)
            tmp[[c(1,2)]] <- base[[1]]
            tmp[[c(1,3)]] <- l[[i]]
            base[[1]] <- tmp[[1]]
        }
    }
    base[[1]]    
}

variables_to_formula <- function(vars,intercept=0) {
    rhs <- add_expressions(vars[-1])
    if(intercept==1) {
        res <- formula(a~b-1)
        res[[2]] <- vars[[1]]
        res[[c(3,2)]] <- rhs        
    }
    else {
        res <- formula(a~b)
        res[[2]] <- vars[[1]]
        res[[3]] <- rhs
    }
    res
}

find_mls_terms <- function(term.name,vars) {
    res<-sapply(vars, function(l) {
        if(length(l)>1) {
            if(as.character(l[[1]])%in%c("mls","fmls","dmls")) {
                ifelse(as.character(l[[2]])==term.name, TRUE,FALSE)
            }
            else FALSE
        }
        else FALSE
    })
    which(res)
}

