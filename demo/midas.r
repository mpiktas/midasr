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
