\docType{data}
\name{oos_prec}
\alias{oos_prec}
\title{Out-of-sample prediction precision data on simulation example}
\format{A \code{data.frame} object with four columns. The first column indicates the sample size, the second the type of constraint, the third the value of the precision measure and the fourth the type of precision measure.}
\description{
The code in the example generates the out-of-sample
prediction precision data for correctly and incorrectly
constrained MIDAS regression model compared to
unconstrained MIDAS regression model.
}
\examples{
## Do not run:
## set.seed(1001)

## gendata<-function(n) {
##     trend<-c(1:n)
##     z<-rnorm(12*n)
##     fn.z <- nealmon(p=c(2,0.5,-0.1),d=17)
##     y<-2+0.1*trend+mls(z,0:16,12)\%*\%fn.z+rnorm(n)
##     list(y=as.numeric(y),z=z,trend=trend)
## }

## nn <- c(50,100,200,300,500,750,1000)

## data_sets <- lapply(n,gendata)

## mse <- function(x) {
##     mean(residuals(x)^2)
## }

## bnorm <- function(x) {
##     sqrt(sum((midas_coef(x)-c(2,0.1,nealmon(p=c(2,0.5,-0.1),d=17)))^2))
## }

## rep1 <- function(n) {
##     dt <- gendata(round(1.25*n))
##     ni <- n
##     ind <- 1:ni
##     mind <- 1:(ni*12)
##     indt<-list(y=dt$y[ind],z=dt$z[mind],trend=dt$trend[ind])
##     outdt <- list(y=dt$y[-ind],z=dt$z[-mind],trend=dt$trend[-ind])
##     um <- midas_r(y~trend+mls(z,0:16,12),data=indt,start=NULL)
##     nm <- midas_r(y~trend+mls(z,0:16,12,nealmon),data=indt,start=list(z=c(1,-1,0)))
##     am <- midas_r(y~trend+mls(z,0:16,12,almonp),data=indt,start=list(z=c(1,0,0,0)))
##     modl <- list(um,nm,am)
##     names(modl) <- c("um","nm","am")
##     list(norms=sapply(modl,bnorm),
##          mse=sapply(modl,function(mod)mean((forecast(mod,newdata=outdt)-outdt$y)^2)))
## }

## repr <- function(n,R) {
##     cc <- lapply(1:R,function(i)rep1(n))
##     list(norms=t(sapply(cc,"[[","norms")),mse=t(sapply(cc,"[[","mse")))
## }

## res <- lapply(nn,repr,R=1000)

## norms <- data.frame(nn,t(sapply(lapply(res,"[[","norms"),function(l)apply(l,2,mean))))
## mses <- data.frame(nn,t(sapply(lapply(res,"[[","mse"),function(l)apply(l,2,mean))))


## msd <- melt(mses[-1,],id=1)
## colnames(msd)[2] <- "Constraint"
## nmd <- melt(norms[-1,],id=1)
## colnames(nmd)[2] <- "Constraint"

## msd$Type <- "Mean squared error"
## nmd$Type <- "Distance from true values"
## oos_prec <- rbind(msd,nmd)
## oos_prec$Type <- factor(oos_prec$Type,levels=c("Mean squared error","Distance from true values"))
}
\keyword{datasets}

