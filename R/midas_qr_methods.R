##' @export
##' @method print midas_r
print.midas_qr <- function(x, digits=max(3,getOption("digits")-3), ...) {
    print.midas_r(x)
}
