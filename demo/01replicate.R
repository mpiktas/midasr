###Replicates the table of the article
#library(midasr)
rm(list=ls())
set.seed(13)

library(devtools)

load_all("midasr")

theta.h0 <- function(p, dk) {
    i <- (1:dk-1)/100
    pol <- p[3]*i + p[4]*i^2
    (p[1] + p[2]*i)*exp(pol)
}

grad.h0<-function(p, dk){
    alpha <- p[1]
    beta <- p[2]
    lambda <- c(p[3],p[4])
    index <- c(1:dk)
    i <- (index-1)/100
    pol <- poly(i,2,raw=TRUE) %*%lambda
    a <- (alpha+beta*i)*exp(pol)
    cbind(a, a*i, a*i*(alpha+beta*i), a*i^2*(alpha+beta*i))
}


sd.x <- 1
sd.y <- 1

n <- 500
m <- 12

n.x.max <- (n+1000)*m

k0 <- 3

theta0 <- theta.h0(c(-0.1,10,-10,-10),(k0+1)*m)

x <- simplearma.sim(list(ar=0.6),n.x.max,sd.x,m)

y <- midas.sim(n,theta0,x,sd.y)

y <- y[-k0:-1]

yx <- model.matrix.midas(y,x,k0)

mu <- midas.u(y,x,k0)

mr <- midas.r(y,x,theta.h0,c(-0.1,10,-10,-10),dk=(k0+1)*m)

hAh.test(mu,mr,grad.h0,dk=(k0+1)*m)
