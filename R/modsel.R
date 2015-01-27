if(getRversion() >= "2.15.1")  utils::globalVariables("X")

##' Create a high frequency lag selection table for MIDAS regression model
##'
##' Creates a high frequency lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation
##' @param from a named list, or named vector with lag numbers which are the beginings of MIDAS lag structures. The names should correspond to the MIDAS lag terms in the formula for which to do the lag selection. Value NA indicates lag start at zero
##' @param to a named list where each element is a vector with two elements. The first element is the lag number from which the lag selection starts, the second is the lag number at which the lag selection ends. NA indicates lowest (highest) lag numbers possible.
##' @param IC the information criteria which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param weight_gradients see \link{midas_r}
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
##' mlr <- hf_lags_table(y~trend+fmls(x,12,12,nealmon),
##'                      start=list(x=rep(0,3)),
##'                      from=c(x=0),to=list(x=c(4,4)))
##' mlr
##'
##' @details This function estimates models sequentially increasing the midas lag from \code{kmin} to \code{kmax} of the last term of the given formula
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
hf_lags_table<- function(formula,data,start,from,to,IC=c("AIC","BIC"),test=c("hAh_test"),Ofunction="optim",weight_gradients=NULL,...) {

    if(!identical(names(from),names(to)))stop("The names of lag structure start and end should be identical")
    from <- as.list(from)
    varnames <- names(from)
    
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())

    kmax <- round(10*log(length(prep$y),base=10))

    lti <- lapply(varnames,function(nm)term_info(prep$mt,nm,prep$Zenv))
    names(lti) <- varnames
   
    table <- mapply(function(ti,st,end) {
        wstart <- start[ti$varname]
        names(wstart) <- ti$weight
        if(is.na(st))st <- 0
        nop <- length(wstart[[1]])
        if(is.na(end[1])) {
            end[1] <- st+nop
        }
        else {
            end[1] <- max(end[1],st+nop)
        }
        if(is.na(end[2])) {
            end[2] <- kmax
        }
        else {
            if(end[1]>end[2])stop("The lag number which ends the selection should be larger or equal to the lag number which starts the selction")
        }        
        expand_weights_lags(ti$weight,st,end,1,wstart)
    },lti,from,to,SIMPLIFY=FALSE)
    names(table) <- varnames
    
    start <- start[!(names(start)%in% varnames)]
    if(length(start)==0)start <- NULL
  
    midas_r_ic_table(formula,data,start=start,table=table,IC=IC,test=test,Ofunction=Ofunction,weight_gradients=weight_gradients,...)
}

##' Create a low frequency lag selection table for MIDAS regression model
##'
##' Creates a low frequency lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation
##' @param from a named list, or named vector with high frequency (NB!) lag numbers which are the beginnings of MIDAS lag structures. The names should correspond to the MIDAS lag terms in the formula for which to do the lag selection. Value NA indicates lag start at zero
##' @param to a named list where each element is a vector with two elements. The first element is the low frequency lag number from which the lag selection starts, the second is the low frequency lag number at which the lag selection ends. NA indicates lowest (highest) lag numbers possible.
##' @param IC the information criteria which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param weight_gradients see \link{midas_r}
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
##' mlr <- lf_lags_table(y~trend+fmls(x,12,12,nealmon),
##'                      start=list(x=rep(0,3)),
##'                      from=c(x=0),to=list(x=c(3,4)))
##' mlr
##'
##' @details This function estimates models sequentially increasing the midas lag from \code{kmin} to \code{kmax} of the last term of the given formula
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
lf_lags_table <- function(formula,data,start,from,to,IC=c("AIC","BIC"),test=c("hAh_test"),Ofunction="optim",weight_gradients=NULL,...) {

    if(!identical(names(from),names(to)))stop("The names of lag structure start and end should be identical")
    from <- as.list(from)
    varnames <- names(from)
    
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())

    kmax <- round(10*log(length(prep$y),base=10))

    lti <- lapply(varnames,function(nm)term_info(prep$mt,nm,prep$Zenv))
    names(lti) <- varnames
    
    table <- mapply(function(ti,st,end) {
        wstart <- start[ti$varname]
        names(wstart) <- ti$weight
        m <- ti$frequency
        if(is.na(st))st <- 0
        nop <- length(wstart[[1]])
        if(is.na(end[1])) {
            end[1] <- max((st+nop)%/%m,1)
        }
        else {
            end[1] <- max(end[1],max((st+nop)%/%m,1))
        }
        if(is.na(end[2])) {
            end[2] <- kmax %/% m
        }
        else {
            if(end[1]>end[2])stop("The lag number which ends the selection should be larger or equal to the lag number which starts the selction")
        }        
        expand_weights_lags(ti$weight,st,end,m,wstart)
    },lti,from,to,SIMPLIFY=FALSE)
    names(table) <- varnames

    start <- start[!(names(start)%in% varnames)]
    if(length(start)==0)start <- NULL
        
    midas_r_ic_table(formula,data,start=start,table=table,IC=IC,test=test,Ofunction=Ofunction,weight_gradients=weight_gradients,...)
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
    term.no <- find_mls_terms(term.name,vars)

    if(length(term.no)==0) stop("No mls terms for variable ", term.name)
    if(length(term.no)>1) stop("There can be only one mls term for variable ", term.name)
    
    mls_term <- vars[[term.no]]
    
    lags <- as.numeric(eval(mls_term[[3]],Zenv))
    freq <- as.numeric(eval(mls_term[[4]],Zenv))
    weightname <- as.character(mls_term[[5]])
    mtype <- as.character(mls_term[[1]])
      
    if(mtype=="fmls")lags <- 0:lags

    list(lags=lags,weight=weightname,varname=as.character(mls_term[[2]]),frequency=freq)  
}


##' Select the model based on given information criteria
##'
##' Selects the model with minimum of given information criteria and model type
##' @param x and output from iclagtab function
##' @param IC the name of information criteria to base the choosing of the model
##' @param test the name of the test for which to print out the p-value
##' @param type the type of MIDAS model, either restricted or unrestricted
##' @param print logical, if TRUE, prints the summary of the best model.
##' @return (invisibly) the best model based on information criteria, \link{midas_r} object
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
##' mhfr <- hf_lags_table(y~trend+fmls(x,12,12,nealmon),
##'                       start=list(x=rep(0,3)),
##'                       from=list(x=0),to=list(x=c(4,6)))
##' 
##' mlfr <- lf_lags_table(y~trend+fmls(x,12,12,nealmon),
##'                       start=list(x=rep(0,3)),
##'                       from=list(x=0),to=list(x=c(2,3)))
##'
##' modsel(mhfr,"BIC","unrestricted")
##'
##' modsel(mlfr,"BIC","unrestricted")
##' 
##' @details This function selects the model from the model selection table for which the chosen information criteria achieves the smallest value. The function works with model tables produced by functions \link{lf_lags_table}, \link{hf_lags_table}, \link{amidas_table} and \link{midas_r_ic_table}.
modsel <- function(x,IC=x$IC[1],test=x$test[1],type=c("restricted","unrestricted"),print=TRUE) {
    if(!(IC%in%x$IC))stop("The supplied information criteria was not used in creating lag selection table")
    type <- match.arg(type)
    coln <- paste(IC,type,sep=".")
    i <- which.min(x$table[,coln])
    
    if(print) {
        cat("\n Selected model with ")
    
        cat(IC," = ",x$table[i,coln],"")
        cat("\n Based on",type, "MIDAS regression model\n")
        cat(" The p-value for the null hypothesis of the test", test, "is", x$table[i,paste(test,"p.value",sep=".")],"\n")
        print(summary(x$candlist[[i]]))
    }
    invisible(x$candlist[[i]])
}

##' Create a weight function selection table for MIDAS regression model
##'
##' Creates a weight function selection table for MIDAS regression model with given information criteria and weight functions.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation
##' @param IC the information criteria which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param weight_gradients see \link{midas_r}
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
##' mwr <- weights_table(y~trend+fmls(x,12,12,nealmon),
##'                      start=list(x=list(nealmon=rep(0,3),
##'                      nbeta=c(1,1,1,0))))
##'
##' mwr
##'
##' @details This function estimates models sequentially increasing the midas lag from \code{kmin} to \code{kmax} of the last term of the given formula
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
weights_table <- function(formula,data,start=NULL,IC=c("AIC","BIC"),test=c("hAh_test"),Ofunction="optim",weight_gradients=NULL,...) {
    
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())
    
    varnames <- names(start)[sapply(start,is.list)]
    
    table <- lapply(varnames,function(nm) {
        ti <- term_info(terms(formula),nm,prep$Zenv)
        out <- expand_weights_lags(names(start[[nm]]),from=1,to=c(1,1),0,start=start[[nm]])
        out$lags <- rep(list(ti$lags),length.out=length(out$weights))
        out
    })
    names(table) <- varnames
    start <- start[!(names(start) %in% varnames)]
    if(length(start)==0)start <- NULL
    midas_r_ic_table(formula,data,start=start,table=table,IC=IC,test=test,Ofunction=Ofunction,weight_gradients=weight_gradients,...)
}


##' Create a weight and lag selection table for MIDAS regression model
##'
##' Creates a weight and lag selection table for MIDAS regression model with given information criteria and minimum and maximum lags.
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param start the starting values for optimisation excluding the starting values for the last term
##' @param table an wls_table object, see \link{expand_weights_lags}
##' @param IC the names of information criteria which to compute
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param weight_gradients see \link{midas_r}
##' @param show_progress logical, TRUE to show progress bar, FALSE for silent evaluation
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
##' mwlr <- midas_r_ic_table(y~trend+fmls(x,12,12,nealmon),
##'                    table=list(x=list(weights=
##'                    as.list(c("nealmon","nealmon","nbeta")),
##'                    lags=list(0:4,0:5,0:6),
##'                    starts=list(rep(0,3),rep(0,3,),c(1,1,1,0)))))
##'
##' mwlr
##'
##' @details This function estimates models sequentially increasing the midas lag from \code{kmin} to \code{kmax} and varying the weights of the last term of the given formula 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
midas_r_ic_table <- function(formula,data=NULL,start=NULL,table,IC=c("AIC","BIC"),test=c("hAh_test"),Ofunction="optim",weight_gradients=NULL,show_progress=TRUE,...) {
    
    Zenv <- new.env(parent=environment(formula))
    formula <- as.formula(formula)
    args <- list(...)
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- formula

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())

    Zenv <- prep$Zenv
    mff <- prep$mf
    isstar <- any(sapply(table,with,any(names(weights)=="*")))
    
    ##Remove those formulas for which the number of parameters is less or equal than number of lags.
    remove_incomplete <- function(info,nm) {
        cond <- mapply(
            function(lags,start){
                ifelse(length(lags)>=length(start),TRUE,FALSE)
            },
            lapply(info$lags,function(x)x[[nm]]),
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
            res <- mapply(function(f,s,lg) {
                out <- formula_table(terms(f),names(table)[i],Zenv,table[[i]],s)
                out$lags <- lapply(out$lags,function(ll)c(ll,lg))
                out
            },wlinfo$formulas,wlinfo$starts,wlinfo$lags,SIMPLIFY=FALSE)

            wlinfo <- combine(res)
            wlinfo <- remove_incomplete(wlinfo,names(table)[i])
        }        
    }
        
    modellist <- mapply(function(f,st) {    
        mff[[2L]] <- f
        ###Add condition for catching the star in the table
        if(isstar) {
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
  
    #maxlag <- which.max(sapply(wlinfo$lags,function(ll)max(sapply(ll,max))))
    maxlag <- which.min(sapply(modellist,with,nrow(X)))
    lrn <- rownames(modellist[[maxlag]]$X)

    modellist <- lapply(modellist,function(mm) {
        rn <- rownames(mm$X)
        ind <- match(lrn,rn)
        mm$y <- mm$y[ind]
        mm$X <- mm$X[ind,]
        mm
    })

    new_cl <- cl
    m <- match(c("table","IC","test","show_progress"),names(new_cl))
    new_cl[m] <- NULL
    
    if("Ofunction" %in% names(new_cl)) new_cl$Ofunction<- eval(new_cl$Ofunction)
    if("weight_gradients" %in% names(new_cl)) new_cl$weight_gradients <- eval(new_cl$weight_gradients)
    if(is.null(eval(new_cl$data))) new_cl$data <- NULL
    
    mrm <- lapply(modellist,function(mm) {
        cll <- new_cl
        cll[[1]] <- as.name("midas_r")
        cll$formula <- formula(mm$mt)
        cll$start <- mm$start
        res <- prepmidas_r(mm$y,mm$X,mm$mt,Zenv,cll,args,mm$start,Ofunction,weight_gradients,mm$itr$lagsTable)
        class(res) <- "midas_r"
        res
    })
   
    if(show_progress) {
        cat("\nModel selection progress:\n")
        pb <- txtProgressBar(min=0,max=length(mrm),initial=0,style=3)
    }
    candlist <- mapply(function(l,i){
       if(show_progress) setTxtProgressBar(pb, i)
        out <- try(midas_r.fit(l),silent=TRUE)
        out        
    },mrm,1:length(mrm),SIMPLIFY=FALSE)
    if(show_progress)close(pb)   
    success <- sapply(candlist,class)

    if("try-error" %in% success)
    for(i in which(success=="try-error")) {
        cat("\n====================================\n")
        cat("The following model did not converge:\n")
        cat(gsub("[ ]+"," ",capture.output(cat(deparse(formula(mrm[[i]]))))))
        cat("\nWith the following error message:\n")
        cat(candlist[[i]])
        
        cat("--------------------------------------\n")        
    }
    candl <- candlist[success!="try-error"]
    
    make_ic_table(candl,IC,test)
}

##' @method update midas_r_ic_table
##' @export
update.midas_r_ic_table <- function(object,...) {
    do.call("make_ic_table",object[-1])
}

make_ic_table <- function(candlist,IC,test,...) {
    makelist <- function(x) {
        if(length(x)==1)list(x)
        else as.list(x)
    }
    makefun <- function(l) {
        lapply(l,function(ll) {
            if(is.function(ll))ll
            else eval(as.name(ll))
        })
    }
    ICfun <- makefun(makelist(IC))
    tfun <- makefun(makelist(test))

    dtest <- function(mm) {
        ##The error might be given if there are NaN values in the hessian, give NA instead of the error.
        out <- try(deriv_tests(mm),silent=TRUE)
        if(class(out)=="try-error")c(NA,NA)
        else unlist(out[c("first","second")])
    }
    tab <- lapply(candlist,function(mm) {
        c(sapply(ICfun,function(ic)ic(mm)),
          sapply(ICfun,function(ic)ifelse(is.null(mm$unrestricted),NA,ic(mm$unrestricted))),
          sapply(tfun,function(tt){
              tst <- try(tt(mm))
              ifelse(class(tst)=="try-error",NA,
              tst$p.value)
          }),
          dtest(mm),
          mm$convergence
          )        
    })
       
    tab <- do.call("rbind",tab)

    colnames(tab) <- c(paste(IC,"restricted",sep="."),paste(IC,"unrestricted",sep="."),paste(test,"p.value",sep="."),"First","Second","Convergence")

    tab <- data.frame(model=sapply(candlist,function(mod) {
        gsub("[ ]+"," ",capture.output(cat(deparse(formula(mod)))))
    }),tab)
    tab$First <- as.logical(tab$First)
    tab$Second <- as.logical(tab$Second)

    res <- list(table=tab,candlist=candlist,IC=IC,test=test,weights=tab[,1],lags=tab[,2])
    class(res) <- "midas_r_ic_table"
    res
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
    
    vars[[term.no]][[1]] <- as.name("mls")
    
    formulas <- vector("list",length(table$lags))
    starts <- formulas
    for(i in 1:length(formulas)) {
        res <- vars
        lt <- res[[term.no]]
        wght <- table$weights[[i]]
        if(is.character(wght)) {
            if(wght=="") lt <- lt[1:4]
            else {
                if(wght=="*")  lt[[5]] <- wght
                else lt[[5]] <- as.name(wght)
            }
        }
        else {
            if(!is.function(wght))stop("Supply either function name or a function")
            nmwght <- names(table$weights)[i]
            lt[[5]] <- as.name(nmwght)
            ##A bit of nasty hack. We rely on the fact that environments are special R objects, i.e. they are not replicated when passed to function.
            if(nmwght %in% ls(envir=Zenv))warning("Name ", nmwght, " is reserved, overwriting with special weight function")
            assign(nmwght,wght,Zenv)
        }
        lt[[3]] <- table$lags[[i]]
        res[[term.no]] <- lt
        formulas[[i]] <- variables_to_formula(res)
        
        wst <- list(table$starts[[i]])
        names(wst) <- varname
        starts[[i]] <- c(start,wst)
    }
    
    list(formulas=formulas,
         lags=lapply(table$lags,function(l){out <- list(l);names(out)<-varname;out}),
         weights=table$weights,
         starts=starts)
}

prepare_model_frame <- function(data,Zenv,cl,mf,pf) {
##Get the response of the model to get the number of observations
##Get the model.frame object, not evaluated!
##Prepare data if necessary    
    if(is.null(data)) {
        ee <- NULL
    } else {
        ee <- data_to_env(data)
        parent.env(ee) <- pf
    }
    assign("ee",ee,Zenv)
    
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
    resf <- eval(mf$formula,Zenv)
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
##' @param to a vector of length two, containing minimum and maximum lags, high frequency if \code{m=1}, low frequency otherwise.
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
    if(is.null(names(start)))names(start) <- rep("",length(start))
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

    normalize_starts <- function(x) {
        inds <- which(names(x$weights) %in% c("*",""))
        for (i in inds ) {
            x$starts[[i]] <- rep(x$starts[[i]],length.out=length(x$lags[[i]]))
        }
        x
    }
    ##For "" and "*" weights replicate the starts according to lags
    out <- normalize_starts(list(weights=weights,lags=lags,starts=starts))
    class(out) <- "lws_table"
    out
}

##' @export
##' @method print lws_table
print.lws_table <- function(x,...) {
    if(is.null(names(x)))names(x) <- c("weights","lags","starts")
    p1 <- function(x)capture.output(cat(deparse(x)))
    p2 <- function(x)capture.output(cat(deparse(round(x,4))))
    print(data.frame(weights=names(x$weights),lags=sapply(x$lags,p1),starts=sapply(x$starts,p2)))
}

##' Create table of weights, lags and starting values for Ghysels weight schema, see \link{amweights}
##'
##' Given weight function creates lags starting from \code{kmin} to \code{kmax} and replicates starting values for each low frequency lag.
##' @title Create table of weights, lags and starting values for Ghysels weight schema
##' @param weight the names of weight functions
##' @param type the type of Ghysels schema, \code{"A"}, \code{"B"} or \code{"C"}
##' @param from the high frequency lags from which to start the fitting
##' @param to to a vector of length two, containing minimum and maximum lags, high frequency if \code{m=1}, low frequency otherwise.
##' @param m the frequency ratio
##' @param start the starting values for the weights of the one low frequency lag
##' @return a \code{lws_table} object, a list with elements \code{weights}, \code{lags} and \code{starts}
##' @examples
##' expand_amidas("nealmon","A",0,c(1,2),12,c(0,0,0))
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
expand_amidas <- function(weight,type=c("A","B","C"),from=0,to,m,start) {
    lags <- lapply(to[1]:to[2],function(x)from:(x*m-1))
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
    ff <- expression(amweights(p,d,m,weight,type))
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
##' nlmn <- expand_weights_lags("nealmon",0,c(4,8),1,start=list(nealmon=rep(0,3)))
##' nbt <- expand_weights_lags("nbeta",0,c(4,8),1,start=list(nbeta=rep(0,4)))
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
##' Create weight and lag selection table for the aggregates based MIDAS regression model
##'
##' This function estimates models sequentialy increasing the midas lag from \code{kmin} to \code{kmax} and varying the weights of the last term of the given formula 
##' @title Weight and lag selection table for aggregates based MIDAS regression model 
##' @param formula the formula for MIDAS regression, the lag selection is performed for the last MIDAS lag term in the formula
##' @param data a list containing data with mixed frequencies
##' @param weights the names of weights used in Ghysels schema
##' @param wstart the starting values for the weights of the firs low frequency lag
##' @param type the type of Ghysels schema see \link{amweights}, can be a vector of types
##' @param start the starting values for optimisation excluding the starting values for the last term
##' @param from a named list, or named vector with high frequency (NB!) lag numbers which are the beginnings of MIDAS lag structures. The names should correspond to the MIDAS lag terms in the formula for which to do the lag selection. Value NA indicates lag start at zero
##' @param to to a named list where each element is a vector with two elements. The first element is the low frequency lag number from which the lag selection starts, the second is the low frequency lag number at which the lag selection ends. NA indicates lowest (highest) lag numbers possible.
##' @param IC the names of information criteria which should be calculated
##' @param test the names of statistical tests to perform on restricted model, p-values are reported in the columns of model selection table
##' @param Ofunction see \link{midasr}
##' @param weight_gradients see \link{midas_r}
##' @param ... additional parameters to optimisation function, see \link{midas_r}
##' @return a \code{midas_r_ic_table} object which is the list with the following elements:
##'
##' \item{table}{the table where each row contains calculated information criteria for both restricted and unrestricted MIDAS regression model with given lag structure}
##' \item{candlist}{the list containing fitted models}
##' \item{IC}{the argument IC}
##' \item{test}{the argument test}
##' \item{weights}{the names of weight functions}
##' \item{lags}{the lags used in models}
##' @examples
##'
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##' 
##' tb <- amidas_table(y~trend+fmls(x,12,12,nealmon),
##'                    data=list(y=y,x=x,trend=trend),
##'                    weights=c("nealmon"),wstart=list(nealmon=c(0,0,0)),
##'                    start=list(trend=1),type=c("A"),
##'                    from=0,to=c(1,2))
##'
##' 
##' @details This function estimates models sequentially increasing the midas lag from \code{kmin} to \code{kmax} and varying the weights of the last term of the given formula 
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
amidas_table <- function(formula,data,weights,wstart,type,start=NULL,from,to,IC=c("AIC","BIC"),test=c("hAh_test"),Ofunction="optim",weight_gradients=NULL,...) {
    Zenv <- new.env(parent=environment(formula))
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(missing(data))data <- NULL

    prep <- prepare_model_frame(data,Zenv,cl,mf,parent.frame())

    lti <- last_term_info(formula,prep$Zenv)
    lags <- lti$lags
 
    m <- lti$frequency
    
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
        expand_amidas(weight=w,type=t,from=from,to=to,m=m,start=s)
    },as.list(weights),as.list(type),as.list(wstart),SIMPLIFY=FALSE)
    names(tb) <- NULL
    
    table <- list(do.call("+",tb))
    names(table) <- lti$varname
    midas_r_ic_table(formula,data,start=start,table=table,IC=IC,test=test,Ofunction=Ofunction,weight_gradients = NULL,...)
}
    

add_expressions <- function(l) {
    if(!is.list(l))stop("The summands must be in a list")
    if(length(l)==1) {
        return(l[[1]])
    }
    else {
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
        return(base[[1]])
    }
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

##' Creates tables for different forecast horizons and table for combined forecasts
##'
##' Divide data into in-sample and out-of-sample. Fit different forecasting horizons for in-sample data. Calculate accuracy measures for individual and average forecasts.
##' 
##' @title Create table for different forecast horizons
##' @param formula initial formula for the 
##' @param data list of data
##' @param from a named list of starts of lags from where to fit. Denotes the horizon
##' @param to a named list for lag selections
##' @param insample the low frequency indexes for in-sample data
##' @param outsample the low frequency indexes for out-of-sample data
##' @param weights names of weight function candidates
##' @param wstart starting values for weight functions
##' @param start other starting values
##' @param IC name of information criteria to choose model from
##' @param seltype argument to modsel, \code{"restricted"} for model selection based on information criteria of restricted MIDAS model, \code{"unrestricted"} for model selection based on unrestricted (U-MIDAS) model.
##' @param ftype which type of forecast to use. 
##' @param test argument to modsel
##' @param measures the names of goodness of fit measures
##' @param fweights names of weighting schemes
##' @param ... additional arguments for optimisation method, see \link{midas_r}
##' @return  a list containing forecasts, tables of accuracy measures and the list with selected models
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
##' @examples
##' ### Sets a seed for RNG ###
##' set.seed(1001)  
##' ## Number of low-frequency observations
##' n<-250
##' ## Linear trend and higher-frequency explanatory variables (e.g. quarterly and monthly)
##' trend<-c(1:n)
##' x<-rnorm(4*n)
##' z<-rnorm(12*n)
##' ## Exponential Almon polynomial constraint-consistent coefficients
##' fn.x <- nealmon(p=c(1,-0.5),d=8)
##' fn.z <- nealmon(p=c(2,0.5,-0.1),d=17)
##' ## Simulated low-frequency series (e.g. yearly)
##' y<-2+0.1*trend+mls(x,0:7,4)%*%fn.x+mls(z,0:16,12)%*%fn.z+rnorm(n)
##' ##Do not run
##' ## cbfc<-select_and_forecast(y~trend+mls(x,0,4)+mls(z,0,12),
##' ## from=list(x=c(4,8,12),z=c(12,24,36)),
##' ## to=list(x=rbind(c(14,19),c(18,23),c(22,27)),z=rbind(c(22,27),c(34,39),c(46,51))),
##' ## insample=1:200,outsample=201:250,
##' ## weights=list(x=c("nealmon","almonp"),z=c("nealmon","almonp")),
##' ## wstart=list(nealmon=rep(1,3),almonp=rep(1,3)),
##' ## IC="AIC",
##' ## seltype="restricted",
##' ## ftype="fixed",
##' ## measures=c("MSE","MAPE","MASE"),
##' ## fweights=c("EW","BICW","MSFE","DMSFE")
##' ## )
##' 
select_and_forecast<- function(formula,data,from,to,
                               insample,outsample,
                               weights,wstart,start=NULL,
                               IC="AIC",
                               seltype=c("restricted","unrestricted"),
                               test="hAh_test",
                               ftype=c("fixed","recursive","rolling"),
                               measures=c("MSE","MAPE","MASE"),
                               fweights=c("EW","BICW","MSFE","DMSFE"),
                               ...) {
    seltype <- match.arg(seltype)
    ftype <- match.arg(ftype)
    if(length(setdiff(fweights,c("EW","BICW","MSFE","DMSFE")))>0) {
        stop("Supported weight schemes are EW, BICW, MSFE, DMSFE")
    }
    
    Zenv <- new.env(parent=environment(formula))
    formula <- as.formula(formula)

    if(missing(data)||is.null(data)) dataenv <- Zenv
    else dataenv <- data_to_env(data)

    ##Change this, there is a cleaner way
    mt <- terms(formula)
    yname <- all.vars(mt[[2]])       
    m <- get_frequency_info(mt, Zenv)
    nms <- names(m)
    fullsample <- lapply(nms,function(nm)eval(as.name(nm),dataenv))
    names(fullsample) <- nms
    
    nmx <- names(fullsample)
    nmx <- nmx[nmx!=yname]
    
    datasplit <- split_data(fullsample,insample,outsample)
    indata <- datasplit$indata
    outdata <- datasplit$outdata
    
    wperm <- do.call("expand.grid",weights)
    nperm <- colnames(wperm)
    
    fhtab <- vector("list",length(from[[1]]))
    for(h in 1:length(from[[1]])) {
        res <- vector("list",nrow(wperm))
        for(i in 1:nrow(wperm)) {
            res[[i]] <- lapply(nperm,function(nm){
                wname <- as.character(wperm[i,nm])
                expand_weights_lags(wname,from[[nm]][h],to[[nm]][h,],1,start=wstart[wname])
            })
            names(res[[i]]) <- nperm
        }
        fhtab[[h]] <- res
    }

    modno <- sapply(fhtab,length)
    nmodno <- sum(modno)
    cat("\nModel selection progress:\n")    
    pb <- txtProgressBar(min=0,max=nmodno,initial=0,style=3)
   
    bestm <- mapply(function(fh,prog)
                    mapply(function(tb,i){
                        out <- modsel(midas_r_ic_table(formula,data=indata,start=start,table=tb,IC=IC,test=test,show_progress=FALSE),IC=IC,type=seltype,print=FALSE,...)
                        setTxtProgressBar(pb, i)
                        out
                    },fh,as.list(prog+1:length(fh)),SIMPLIFY=FALSE),
                    fhtab,as.list(c(0,cumsum(modno)[-length(modno)])),SIMPLIFY=FALSE)
    close(pb)   
    
    avgforc <- lapply(bestm,function(l)average_forecast(l,data=fullsample,insample=insample,outsample=outsample,type=ftype,fweights=fweights,measures=measures))
    tabfh <- do.call("rbind",lapply(avgforc,function(l)l$accuracy$individual))

    tboutc <- lapply(avgforc,function(l)l$accuracy$average)
    
    tabh <- do.call("rbind",mapply(function(tb,h){
        data.frame(Horizon=h,tb)
    },tboutc,1:length(tboutc),SIMPLIFY=FALSE))
    
    tabh[order(tabh$Horizon,tabh$Scheme),]
        
    list(accuracy=list(individual=tabfh,average=tabh),bestlist=bestm,forecasts=avgforc)
}

    
lf_range_to_hf <- function(range,m) {
    unlist(lapply(range,function(lf)(lf-1)*m+1:m))
}

MSE <- function(o,p) {
    mean((o-p)^2)
}

MAPE <- function(o,p) {
    mean(abs((o-p)/o)*100)
}

MASE <- function(o,p) {
    mean(abs(o-p)/mean(abs(diff(o))))
}
##' Splits mixed frequency data into in-sample and out-of-sample datasets given the indexes of the low frequency data
##'
##' It is assumed that data is a list containing mixed frequency data. Then given the indexes of the low frequency data the function splits the data into two subsets. 
##' @title Split mixed frequency data into in-sample and out-of-sample
##' @param data a list containing mixed frequency data 
##' @param insample the low frequency indexes for in-sample data
##' @param outsample the low frequency indexes for out-of-sample data
##' @return a list with elements \code{indata} and \code{outdata} containing respectively in-sample and out-of-sample data sets
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
##' @examples
##'
##' #Monthly data
##' x <- 1:24
##' #Quartely data
##' z <- 1:8
##' #Yearly data
##' y <- 1:2
##' split_data(list(y=y,x=x,z=z),insample=1,outsample=2)
split_data <- function(data,insample,outsample) {
    fullsample <- data_to_list(data)
    nr <- sapply(fullsample,length)
    m <- nr %/% min(nr)
    indata <- mapply(function(var,freq){
        var[lf_range_to_hf(insample,freq)]
    },fullsample,as.list(m),SIMPLIFY=FALSE)
    
    outdata <- mapply(function(var,freq){
        var[lf_range_to_hf(outsample,freq)]
    },fullsample,as.list(m),SIMPLIFY=FALSE)        

    list(indata=indata,outdata=outdata)
}
##' Average MIDAS model forecasts using specified weighting scheme. Produce in-sample and out-of-sample accuracy measures. 
##'
##' Given the data, split it to in-sample and out-of-sample data. Then given the list of models, reestimate each model with in-sample data and produce out-of-sample forecast. Given the forecasts average them with the specified weighting scheme. Then calculate the accuracy measures for individual and average forecasts.
##'
##' The forecasts can be produced in 3 ways. The \code{"fixed"} forecast uses model estimated with in-sample data. The \code{"rolling"} forecast reestimates model each time by increasing the in-sample by one low frequency observation and dropping the first low frequency observation. These reestimated models then are used to produce out-of-sample forecasts. The \code{"recursive"} forecast differs from \code{"rolling"} that it does not drop observations from the beginning of data.
##' 
##' @title Average forecasts of MIDAS models
##' @param modlist a list of \code{midas_r} objects
##' @param data a list with mixed frequency data
##' @param insample the low frequency indexes for in-sample data
##' @param outsample the low frequency indexes for out-of-sample data
##' @param type a string indicating which type of forecast to use. 
##' @param fweights names of weighting schemes
##' @param measures names of accuracy measures
##' @param show_progress logical, TRUE to show progress bar, FALSE for silent evaluation
##' @return a list containing forecasts and tables of accuracy measures
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @export
##' @examples
##' set.seed(1001)  
##' ## Number of low-frequency observations
##' n<-250
##' ## Linear trend and higher-frequency explanatory variables (e.g. quarterly and monthly)
##' trend<-c(1:n)
##' x<-rnorm(4*n)
##' z<-rnorm(12*n)
##' ## Exponential Almon polynomial constraint-consistent coefficients
##' fn.x <- nealmon(p=c(1,-0.5),d=8)
##' fn.z <- nealmon(p=c(2,0.5,-0.1),d=17)
##' ## Simulated low-frequency series (e.g. yearly)
##' y<-2+0.1*trend+mls(x,0:7,4)%*%fn.x+mls(z,0:16,12)%*%fn.z+rnorm(n)
##' mod1 <- midas_r(y ~ trend + mls(x, 4:14, 4, nealmon) + mls(z, 12:22, 12, nealmon),
##'                 start=list(x=c(10,1,-0.1),z=c(2,-0.1)))
##' mod2 <- midas_r(y ~ trend + mls(x, 4:20, 4, nealmon) + mls(z, 12:25, 12, nealmon),
##'                 start=list(x=c(10,1,-0.1),z=c(2,-0.1)))
##' 
##' ##Calculate average forecasts
##' avgf <- average_forecast(list(mod1,mod2),
##'                         data=list(y=y,x=x,z=z,trend=trend),
##'                         insample=1:200,outsample=201:250,
##'                         type="fixed",                            
##'                         measures=c("MSE","MAPE","MASE"),
##'                         fweights=c("EW","BICW","MSFE","DMSFE"))
average_forecast <- function(modlist,
                             data, insample, outsample,
                             type = c("fixed", "recursive", "rolling"),
                             fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                             measures = c("MSE", "MAPE", "MASE"),
                             show_progress=TRUE) {

    #if(length(modlist)==1)stop("Need more than 1 model to produce average forecasts")
    if(missing(data))stop("Data need to be supplied for forecasting")
    
    type <- match.arg(type)
 
    last_in<- length(insample)
    if(insample[last_in]>outsample[1])stop("The in-sample and out-of-sample indexes should not overlap") 
    if(outsample[1]-insample[last_in]!=1)stop("There should be no gaps between in-sample and out-of-sample indexes")
   
    datasplit <- split_data(data,insample,outsample)

    indata <- datasplit$indata
    outdata <- datasplit$outdata
    
    ##Check that the response variables of the modlist are the same
    ynames <- sapply(modlist,function(mod)all.vars(terms(mod)[[2]]))
    yname <- unique(unlist(ynames))
    if(length(yname)>1)stop("The response variable must be the same for all the models")

    outy <- outdata[[yname]]
    
    reeval <- function(candlist,redata) {
        lapply(candlist,function(mod) {
            ##Setup all the necessary info
            if(inherits(mod,"midas_r_np")) {
                do.call("midas_r_np",list(formula(mod),data=redata),envir=mod$Zenv)
            } else {             
                out <- update(mod, data = redata) #do.call("midas_r",list(formula(mod),data=redata,start=mod$start.list,Ofunction="optim",method="BFGS",control=list(maxit=0)),envir=mod$Zenv)
#            ##Run optimisation with the original model settings
#                out$argmap.opt <- mod$argmap.opt
#                out$start.opt <- coef(mod)
#                do.call("midas_r",list(out,start=coef(mod)),envir=out$Zenv)
                                        #                midas_r.fit(out)
                out
            }
        })
    }

    bestm <- reeval(modlist,indata)
    
    inf <- lapply(bestm,function(mod) cbind(mod$model[,1],fitted(mod)))
       
    if(type=="fixed") {
        outf <- lapply(bestm,function(mod)                  
                       cbind(outdata[[yname]],point_forecast.midas_r(mod,newdata=outdata,method="static"))
                       )       
    }
    else {
        if(type%in%c("rolling")) {
            fulls <- c(insample,outsample)            
        }
        outm <- matrix(NA,nrow=length(outsample),ncol=length(modlist))
        if(show_progress) {
            cat("\nDoing", type, "forecast :\n")    
            pb <- txtProgressBar(min=0,max=length(outsample),initial=0,style=3)
        }
        
        for(i in 1:length(outsample)) {
            newout <- outsample[i]
            if(i>1) {
                if(type=="recursive") newin <- c(insample,outsample[1:(i-1)])
                else {                    
                    newin <- fulls[1:last_in+i-1]
                }
            }
            else newin <- insample
            splitnew <- split_data(data,newin,newout)
            emod <- reeval(modlist,splitnew$indata)
            outm[i,] <- sapply(emod,point_forecast.midas_r,newdata=splitnew$outdata,method="static")
            if(show_progress) setTxtProgressBar(pb, i)
        }
        if(length(modlist)>1) {
            outf <- lapply(data.frame(outm),function(x)cbind(outdata[[yname]],x))}
        else {
            outf <- list(cbind(outdata[[yname]],outm[,1]))
        }
        if(show_progress)close(pb)
    }
    
    msrfun <- lapply(measures,function(msr)eval(as.name(msr)))

    calcmsr <- function(ll) {
        sapply(msrfun,function(msr) {
            sapply(ll,function(a)msr(a[,1],a[,2]))
        })
    }
    outstat <- calcmsr(outf)
    instat <- calcmsr(inf)
    if(is.null(dim(outstat)))outstat <- matrix(outstat,nrow=1)
    if(is.null(dim(instat)))instat <- matrix(instat,nrow=1)
    
    colnames(outstat) <- paste0(measures,".out-of-sample")
    colnames(instat) <- paste0(measures,".in-sample")
    modi <- sapply(bestm,function(mod)gsub("[ ]+"," ",capture.output(cat(deparse(formula(mod))))))

    tabfh <- data.frame(Model=modi,outstat,instat)
    rownames(tabfh) <- NULL
    EW <- function(hh) {
        n <- length(hh)
        rep(1/n,n)
    }
    BICW <- function(hh) {
        bicm <- sapply(hh,BIC)
        bicm <- bicm-min(bicm)
        ebic <- exp(-bicm)
        sebic <- sum(ebic)
        if(sebic==0)EW(hh)
        else ebic/sebic        
    }
    MSFEd <- function(hh,delta) {
        mi <- sapply(hh,function(xx) {
            sum((xx[,1]-xx[,2])^2*delta^(length(outy):1-1))
        })
        imi <- 1/mi
        imi/sum(imi)
    }
    MSFE <-  function(hh) MSFEd(hh,1)
    DMSFE <- function(hh) MSFEd(hh,0.9)


    w1 <- EW(bestm)
    w2 <- BICW(modlist)
    w3 <- MSFE(outf)
    w4 <- DMSFE(outf)

    fc <- sapply(outf,function(m)m[,2])
       
    allwlist <- list(w1,w2,w3,w4)
    names(allwlist) <- c("EW","BICW","MSFE","DMSFE")
    wlist <- allwlist[fweights]
    outc <- lapply(wlist,function(ww)cbind(outf[[1]][,1],apply(fc,1,function(r)sum(r*ww))))
    names(outc) <- fweights 
    
    tboutc <- calcmsr(outc)
    colnames(tboutc) <- measures
        
    tabh <- data.frame(Scheme=fweights,tboutc)
    rownames(tabh) <- NULL

    list(forecast=fc,
         avgforecast=sapply(outc,function(m)m[,2]),
         accuracy=list(
             individual=tabfh,
             average=tabh))
}

