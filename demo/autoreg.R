library(midasr)
theta_h0 <- function(p, dk) {
   i <- (1:dk-1)/100
   pol <- p[3]*i + p[4]*i^2
   (p[1] + p[2]*i)*exp(pol)
}

##Generate coefficients
theta0 <- theta_h0(c(-0.1,10,-10,-10),4*12)


xx <- ts(arima.sim(model = list(ar = 0.6), 3000*12), frequency = 12)

aa <- lapply(c(50,100,200,500,1000,1500,2000), function(n) {
    y <- midas_auto_sim(n, 0.5, xx, theta0, n_start = 200)
    x <- window(xx,start=start(y))
    midas_r(y~mls(y,1,1)+fmls(x,4*12-1,12,theta_h0),start=list(x=c(-0.1,10,-10,-10)))  
})

sapply(aa,function(x)c(nrow(x$model),coef(x)))

bb <- lapply(c(50,100,200,500,1000,1500,2000), function(n) {
  y <- midas_auto_sim(n, c(0.5, 0.1), xx, theta0, n_start = 200)
  x <- window(xx, start=start(y))
  midas_r(y~mls(y,1:2,1)+fmls(x,4*12-1,12,theta_h0),start=list(x=c(-0.1,10,-10,-10)))  
})

sapply(bb,function(x)c(nrow(x$model),coef(x)))
