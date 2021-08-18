##' Check whether non-linear least squares restricted MIDAS regression problem has converged
##'
##' Computes the gradient and hessian of the optimisation function of restricted MIDAS regression and
##' checks whether the conditions of local optimum are met. Numerical estimates are used.
##' @param x \code{\link{midas_r}} object
##' @param tol a tolerance, values below the tolerance are considered zero
##' @return a list with gradient, hessian of optimisation function and convergence message
##' @rdname deriv_tests
##' @seealso midas_r
##' @export
##' @author Vaidotas Zemlys
deriv_tests <- function(x, tol = 1e-6) UseMethod("deriv_tests")

#' @rdname deriv_tests
#' @method deriv_tests midas_r
#' @export
deriv_tests.midas_r <- function(x, tol = 1e-6) {
  gr <- x$gradient(coef(x))
  hess <- x$hessian(coef(x))
  egh <- eigen(hess)$values

  first <- sqrt(sum(gr^2)) < tol
  second <- !any(egh < 0) & !any(abs(egh) < tol)
  list(first = first, second = second, gradient = gr, eigenval = egh)
}
