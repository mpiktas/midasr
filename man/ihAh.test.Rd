\name{ihAh.test}
\alias{ihAh.test}
\title{Test restrictions on coefficients of MIDAS regression with I(1) regressors}
\usage{
  ihAh.test(x)
}
\arguments{
  \item{x}{a MIDAS regression model with restricted
  coefficients, estimated with \code{\link{imidas_r}}}
}
\value{
  a \code{htest} object
}
\description{
  Perform a test whether the restriction on MIDAS
  regression coefficients with I(1) regressors holds.
}
\details{
  The test is chosen depending on the model. For
  \code{"onestep"} MIDAS regression test \code{hAh.test} is
  used. For \code{"twosteps"} a modified version of the
  test \code{ihAh.nls.test} is used. For \code{"reduced"}
  special test \code{ihAh.Td.test} is calculated.
}
\examples{
theta.h0 <- function(p, dk) {
  i <- (1:dk-1)/100
  pol <- p[3]*i + p[4]*i^2
 (p[1] + p[2]*i)*exp(pol)
}
theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)

xx <- simplearma.sim(list(ar=1),1500*12,1,12)
y <- midas.sim(500,theta0,xx,1)
x <- window(xx,start=start(y))

imr <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,start=list(x=c(-0.1,10,-10,-10)))

imr.t0 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="reduced",start=list(x=c(-0.1,10,-10,-10)))

imr.t2 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="twosteps",start=list(x=c(-0.1,10,-10,-10)))

ihAh.test(imr)
#ihAh.test(imr.t0)
#ihAh.test(imr.t2)
}
\author{
  Benediktas Bilinskas, Virmantas Kvedaras, Vaidotas Zemlys
}
\seealso{
  hAh.test, hAhr.test
}

