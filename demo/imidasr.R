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


###Estimate with additional step, with unrestricted last coefficient

###By hand
mu <- midas_u(y~fmls(dx,4*12-1,12)+mls(x,4*12,12)-1)
cf <- "mls(x, 4 * 12, 12)"
u <- mu$model[,1]-mu$model[,cf]*coef(mu)[cf]
uu <- c(rep(NA,length(y)-length(u)),u)
mr.t2 <- midas_r(uu~fmls(dx,4*12-1,12,pp1)-1,start=list(dx=c(-0.1,10,-10,-10)))

##With the special function
imr.t2 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="twosteps",start=list(x=c(-0.1,10,-10,-10)))

###Estimate with additional step, without the last coefficient

###By hand
mu.t0 <- midas_u(y~fmls(dx,4*12-2,12)+mls(x,4*12-1,12)-1)
cf <- "mls(x, 4 * 12 - 1, 12)"
u0 <- mu.t0$model[,1]-mu.t0$model[,cf]*coef(mu.t0)[cf]
uu0 <- c(rep(NA,length(y)-length(u0)),u0)
mr.t0 <- midas_r(uu0~fmls(dx,4*12-2,12,pp1)-1,start=list(dx=c(-0.1,10,-10,-10)))

###With special function
imr.t0 <- imidas_r(y~fmls(x,4*12-1,12,theta.h0)-1,model="reduced",start=list(x=c(-0.1,10,-10,-10)))


###Add additional regressors for the first example

z <- ts(simplearma.sim(list(ar=0.5),length(y),1,frequency(y)),start=start(y))
yy <- y+z*0.5+1

###By hand
mr2 <- midas_r(yy~fmls(dx,4*12-1,12,pp)+mls(x,4*12,12)+z,start=list(dx=c(-0.1,10,-10,-10)))

###With special function
imr2 <- imidas_r(yy~fmls(x,4*12-1,12,theta.h0)+z,start=list(x=c(-0.1,10,-10,-10)))

