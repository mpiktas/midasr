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

##Low freqency variables:
ldt <- data.frame(y=y,trend=1:length(y))

##High frequency variables
hdt <- data.frame(x=x)

###Estimate unrestricted and restricted models with different number of lags

alli <- foreach(k=c(0,1,2,3)) %do% {	
	mu <- midas_u(y~mdslag(x,k)+trend,ldt,hdt)
	mr <- midas_r(y~mdslag(x,k,theta.h0)+trend,ldt,hdt,start=list(theta.h0=c(10,-10,-10,-10)),dk=(k+1)*12)
	mr1 <- midas_r(y~mdslag(x,k,theta.h1)+trend,ldt,hdt,start=list(theta.h1=c(10,-10,-10,-10)),dk=(k+1)*12)
    list(ur=mu,kz=mr,al=mr1,k=k)
}

###Get p-values

sapply(alli,with,c(hAh.test(kz)$p.value,hAh.test(al)$p.value))

###Get coefficients of MIDAS regression

lapply(alli,with,ur$midas.coefficients)
##KZ restriction
lapply(alli,with,kz$midas.coefficients)
##GH restriction
lapply(alli,with,al$midas.coefficients)

###Get restriction parameters

#Kvedaras, Zemlys
sapply(alli,with,coef)
#Almon lag
sapply(alli,with,coef)

##Plot the coefficients

dev.new()
par(mfrow=c(2,2))

lapply(alli,with,{
    plot(coef(ur),xlab="Lags",ylab="Coefficients",main=paste("k = ",k,sep=""))
    points(kz$midas.coefficients,pch=16,col="red")
    points(al$midas.coefficients,pch=16,col="blue")
})
