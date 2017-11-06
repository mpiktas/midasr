##' @export
##' @method print midas_r
print.midas_qr <- function(x, digits=max(3,getOption("digits")-3), ...) {
    print.midas_r(x)
}

##' @export
##' @method coef midas_r
coef.midas_qr <- function(object, midas = FALSE, term_names = NULL, ...) {
    if(is.null(term_names)) coef.midas_r(object, midas, ...)
    else {
        if(length(object$tau) > 1) {
            mc <- object$midas_coefficients
            if(length(setdiff(term_names,names(object$term_info)))>0) stop("Some of the term names are not present in estimated MIDAS regression")
            if(midas) {
                res <- lapply(object$term_info[term_names], function(x) object$midas_coefficients[x$midas_coef_index,])
            } else {
                res <- lapply(object$term_info[term_names], function(x) object$coefficients[x$coef_index,])
            }
            names(res) <- NULL
            if(length(res) ==1) return(res[[1]])
            else return(unlist(res))
        } else coef.midas_r(object, midas, term_names, ...)
    }        
 
}

