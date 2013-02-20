library(midasr)

theta.h0 <- function(p, dk) {
    i <- (1:dk-1)/100
    pol <- p[3]*i + p[4]*i^2
    (p[1] + p[2]*i)*exp(pol)
}

theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)

xx <- simplearma.sim(list(ar=1),1500*12,1,12)
y <- midas.sim(500,theta0,xx,1)

x <- window(xx,start=start(y))
dx <- c(NA,diff(x))

pp1 <- function(p,d)cumsum(theta.h0(p,d))

###Do the transformation by hand
mr <- midas_r(y~fmls(dx,4*12-1,12,pp1)+mls(x,4*12,12)-1,start=list(dx=c(-0.1,10,-10,-10)))

###Use the designated function
imr <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,start=list(x=c(-0.1,10,-10,-10)))

###Test the restriction. The usual test hAh can be used.
hAh.test(imr)
