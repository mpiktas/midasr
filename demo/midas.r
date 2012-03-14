###################
### ANTRAS--MIX ###
###################
rm(list=ls())
cores <- 16
########################################################################
m <- c(12)# m is: L=(k+1)m
k0 <- c(3)# k is: L=(k+1)m
rho <- c(0.6) # egzogeninio x autokoreliacija
########################################################################
n<-500 ##stebejimu skaicius
type<-"ar" ##  x dgp
fix.k <-0 ## 1- fiksuotas k=k0; 0 - k estimated using info
info.type <- 0 ## n^alpha--0, AIC--1, BIC--2, HQ--3, AICC--4
if(info.type==0){laipsn <- 0.33} else {laipsn <- 0}##  kmax=n^laipsn arba Schwert arba:
n.is.kiek <- 10 ## nustatyti kmax=n/(m*n.is.kiek)-1, jei
n.is.kiek.sharp <- 1 ## jei 1 - sharp, jei ne - tik max kontrolei naudojamas n.is.kiek

### funkcijos parametrai ###
a.0<- -0.1
b.0<-10
l.0<-c(-10,-10)
sd.x <- 1
sd.y <- 1

### grid search NLS'e ###
grid <- 0 ### i kiek intervalu sudalinti, jei 0 - nedaromas grid s.
grid.0<- expand.grid(alpha=seq(-0.3, 0.1, len = grid),beta=seq(0.05, 20, len = grid), lambda.1=seq(-30, -10, len = grid),lambda.2=seq(-30, -0.1, len = grid))
kiek.0 <- nrow(grid.0)
#######################################################################
#######################################################################
#######################################################################
librar <- "/home/scratch/lib64/R/library"
library("multicore",lib=librar)
library("iterators",lib=librar)
library("foreach",lib=librar)
library("doMC",lib=librar)
library("MASS",lib=librar)
registerDoMC(cores)
#################
### Funkcijos ###
#################
theta.h0 <- function(index,alpha,beta,lambda.1,lambda.2) {
    i <- (index-1)/100
    pol <- lambda.1*i+lambda.2*i^2
    (alpha+beta*i)*exp(pol)
}
theta.h1 <- function(index,alpha,beta,lambda.1,lambda.2) {
    i <- (index-1)/100
    pol <- lambda.1*i+lambda.2*i^2
    epol <- exp(pol)
    alpha+beta*epol/sum(epol)
}
make.g.h0<-function(coefs,k.par){
    alpha<-coefs[1]
    beta<-coefs[2]
    lambda<-c(coefs[3],coefs[4])
    index <- c(1:k.par)
    i <- (index-1)/100
    pol <- poly(i,2,raw=TRUE) %*%lambda
    a<-(alpha+beta*i)*exp(pol)
    cbind(a,a*i,a*i*(alpha+beta*i),a*i^2*(alpha+beta*i))
}
make.g.h1<-function(coefs,k.par){
    alpha<-coefs[1]
    beta<-coefs[2]
    lambda<-c(coefs[3],coefs[4])
    index <- c(1:k.par)
    i <- (index-1)/100
    pol <- poly(i,2,raw=TRUE) %*%lambda
    epol<-exp(pol)
    b0 <- epol/sum(epol)
    b1 <- sum(epol*i)/sum(epol)
    b2 <- sum(epol*i^2)/sum(epol)
    cbind(index/index,b0,beta*b0*(i-b1),beta*b0*(i^2-b2))
}
dgp.x <- function(n.x,rho,sd,burn.in=300) {
    e <- rnorm(n.x+burn.in,sd=sd.x)
    if(type=="ar"){res <- filter(e,filter=rho,method="recursive",sides=1)}
    if(type=="ma"){res<-c(0,e[2:length(e)]+rho*e[1:(length(e)-1)])}
    res[-burn.in:-1]
}
dgp <- function(n,n.x,kmax,theta,rho,m){
                theta.d <- cbind(theta[1:n.x],1:n.x)[order(-(1:n.x)),1]
                x <- dgp.x(n.x,rho,sd.x)
                idx <- m*c((n.x/m-n+1):(n.x/m))
                X <- foreach(h.x=0:((kmax+1)*m-1), .combine='cbind')%do%{
                         x[idx-h.x]
                     }
                y <-foreach(h.y= 1:n,.combine='c')%do%{
                             t(x[1:(n.x-(n-h.y)*m)])%*%theta.d[((n-h.y)*m+1):n.x]

                       }+rnorm(n,sd=sd.y)
                yx <- cbind(y,X)
            }
gen.IC <- function(kmax,yx,n.n){
    res <- foreach(k.l=0:(kmax-1),.combine=rbind) %do% {##? be -1
        mod <- lm(yx[(k.l+1):n.n,1]~yx[(k.l+1):n.n,2:((k.l+1)*m+1)]-1)
        k.par <- length(coef(mod))
        N <- length(resid(mod))
        sigma2 <- sum(resid(mod)^2)/(N-k.par)
        rezu <- c(k.l,AIC(mod,k=2),AIC(mod,k=log(N)),AIC(mod,k=2*log(log(N))),AIC(mod,k=2*N/(N-k.par-1)))
    }
    if(is.matrix(res)==FALSE){res<-rbind(res,res)}
    colnames(res) <- c("Lag","AIC","BIC","HQ","AICC")
    res
}
###########################################################################
###########################################################################
###########################################################################

###################
### GENERAVIMAS ###
###################
lamb <- c(a.0,b.0,l.0[1],l.0[2]) ## apribojimo funkcijos parametrai
if(laipsn==0) kmax <- k0+1 else kmax <- trunc(n^laipsn)
## max eiles nustatymas, jei k parenkamas pagal IC##
if(fix.k==0){
    kmax <- trunc(n^laipsn)
    if(laipsn==0){kmax <- trunc(m*(n/100)^0.25)}
    if(kmax*m>n/n.is.kiek){kmax <- trunc(n/(m*n.is.kiek))} #virsut.riba
    if(n.is.kiek.sharp==1){kmax <- trunc(n/(m*n.is.kiek))}
}
n.x.max <- (n+1000)*m ##kiek generuoti x stebejimu
theta0 <- theta.h0(1:n.x.max,alpha=a.0,beta=b.0,lambda.1=l.0 [1],lambda.2=l.0[2])## generuota funkcija
isrink <- c((1:((k0+1)*m))/(1:((k0+1)*m)),1:(length(theta0)-(k0+1)*m)*0)
theta <- theta0*isrink
yx <- dgp(n,n.x.max,max(kmax,k0+1),theta,rho,m)
##################
### VERTINIMAS ###
##################
ifelse(fix.k==1, k.x <- k0,{## k=k0, arba k pagal IC
    ifelse(laipsn>0, k.x <- trunc(n^laipsn),{
        IC6<-gen.IC(kmax,yx,n)
        k.x <- IC6[apply(IC6[,2:5],2,which.min)[info.type],1]})
})
n.e <- (k.x+1)*m
###########
### OLS ###
###########
#jei lm
#mod <- lm(yx[(k.x+1):n,1]~yx[(k.x+1):n,2:(n.e+1)]-1)
#se2 <- sum(resid(mod)^2)/(n-n.e-1)
### Duomenu masyvai-be lm ###
X <- yx[(k.x+1):n,2:(n.e+1)]
y <- yx[(k.x+1):n,1]
OLS<-ginv(t(X)%*%X)%*%t(X)%*%y
se2 <- t(y-X%*%OLS)%*%(y-X%*%OLS)/(n-n.e)
XtX<-t(X)%*%X
P<-t(chol(XtX))
###########
### NLS ###
###########
g <- function(lamb) {y-X%*%theta.h0(1:n.e,lamb[1],lamb[2],lamb[3],lamb[4])}
fn0 <- function(lamb) {t(g(lamb))%*%g(lamb)}
opt.h0 <- optim(lamb,fn0,control=list(maxit=10000))
if(grid>0){
    val <- 9e99
    zin <- 0
#    foreach(zin=1:kiek.0)%dopar%{ ### padaryti su foreach
    while(zin<kiek.0){
        zin <- zin+1
        opt <- optim(grid.0[zin,],fn0,control=list(maxit=1000))
        if(val>opt$value){
            opt.h0 <- opt
            val <- opt$value
        }
#        opt.h0
    }
}
e0 <- opt.h0[[1]]
NLS0 <- theta.h0(1:n.e,e0[1],e0[2],e0[3],e0[4])
##############
### Testas ###
##############
h0 <- NA
h0.try <- try({
    D0<-make.g.h0(e0,n.e)
    h.0<-t(P)%*%(OLS-NLS0)
    Delta.0 <- D0%*%ginv(t(D0)%*%XtX%*%D0)%*%t(D0)
    A0<-diag(n.e)-t(P)%*%Delta.0%*%P
    h0<-t(h.0[1:n.e])%*%A0[1:n.e,1:n.e]%*%h.0[1:n.e]/se2
})
vec <- c(h0,k.x,kmax,e0,n.e,n.e-length(e0),n,se2)
names(vec) <- c("H0","k.x","kmax","lamb.1","lamb.2","alpha","beta","n.e","df","n","se2")
vec
