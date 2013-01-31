\name{ihAh.Td.test}
\alias{ihAh.Td.test}
\title{Test restrictions on coefficients of MIDAS regression with I(1) regressors}
\usage{
  ihAh.Td.test(x)
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
  regression coefficients with I(1) regressors holds, in
  case where twostep estimation is used.
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

imr.t0 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="reduced",start=list(x=c(-0.1,10,-10,-10)))

ihAh.Td.test(imr.t0)
}
\author{
  Benediktas Bilinskas, Virmantas Kvedaras, Vaidotas Zemlys
}

