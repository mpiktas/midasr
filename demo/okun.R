rm(list=ls())

library(midasr)
data("USunempr")
data("USrealgdp")

y <- diff(log(USrealgdp))
x <- window(diff(USunempr),start=1949)

##Restriction used in Kvedaras, Zemlys (2012)
theta.h0 <- function(p, dk) {
   i <- (1:dk-1)/100
   (p[1]*i)*exp(p[2]*i + p[3]*i^2+p[4]*i^3)
}

##Almon lag restriction
theta.h1 <- function(p, dk) {
   i <- (1:dk-1)/100
    p[1]*exp(p[2]*i+p[3]*i^2+p[4]*i^3)/sum(exp(p[2]*i+p[3]*i^2+p[4]*i^3))

}

##Additional variables: intercept and trend
exog<-ts(cbind(rep(1,length(y)),1:length(y)),start=start(y),frequency=frequency(y))

##Starting zero values for additional variables
exos<-rep(0,ncol(exog))

###Estimate unrestricted and restricted models with different number of lags

alli <- foreach(k=c(0,1,2,3)) %do% {	
	mu <- midas.u(y,x,exo=exog,k)
	mr <- midas.r(y,x,resfun=theta.h0,exo=exog,start=list(resfun=c(10,-10,-10,-10),exo=exos),dk=(k+1)*12)
	mr1 <- midas.r(y,x,resfun=theta.h1,exo=exog,start=list(resfun=c(10,-10,-10,-10),exo=exos),dk=(k+1)*12)
    list(ur=mu,kz=mr,al=mr1)
}

###Get p-values

sapply(alli,with,c(hAh.test(kz)$p.value,hAh.test(al)$p.value))

###Get coefficients

##Unrestricted
lapply(alli,with,coef(ur))
##KZ restriction
lapply(alli,with,coef(kz))
##GH restriction
lapply(alli,with,coef(al))

###Get restriction parameters

#Kvedaras, Zemlys
sapply(alli,with,kz$parameters)
#Almon lag
sapply(alli,with,al$parameters)

##Plot the coefficients

dev.new()
par(mfrow=c(2,2))

lapply(alli,with,{
    plot(coef(ur))
    points(coef(kz),pch=16,col="red")
    points(coef(al),pch=16,col="blue")
})
