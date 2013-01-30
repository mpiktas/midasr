library(midasr)
theta.h0 <- function(p, dk) {
   i <- (1:dk-1)/100
   pol <- p[3]*i + p[4]*i^2
   (p[1] + p[2]*i)*exp(pol)
}

##Generate coefficients
theta0 <- theta.h0(c(-0.1,10,-10,-10),4*12)


##Generate the predictor variable
xx <- simplearma.sim(list(ar=0.6),3000*12,1,12)

aa <- foreach(n=c(50,100,200,500,1000,1500,2000)) %do% {
    y <- midas.auto.sim(n,theta0,c(0.5),xx,1,n.start=100)
    x <- window(xx,start=start(y))
    midas_r(y~mls(y,1,1)+fmls(x,4*12-1,12,theta.h0),start=list(x=c(-0.1,10,-10,-10)))
}

sapply(aa,function(x)c(nrow(x$model),coef(x)))

bb <- foreach(n=c(50,100,200,500,1000,1500,2000)) %do% {
    y <- midas.auto.sim(n,theta0,c(0.5,0.1),xx,1,n.start=100)
    x <- window(xx,start=start(y))
    midas_r(y~mls(y,1:2,1)+fmls(x,4*12-1,12,theta.h0),start=list(x=c(-0.1,10,-10,-10)))
}

sapply(bb,function(x)c(nrow(x$model),coef(x)))
