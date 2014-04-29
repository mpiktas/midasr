##' Package for estimating and testing MIDAS regression.
##' 
##' The main feature of this package is function \code{\link{hAh.test}}
##' which performs test whether coefficients of MIDAS regression have certain functional form
##' 
##' @name midasr-package
##' @aliases midasr
##' @docType package
##' @title Estimating and testing MIDAS regression
##' @author Virmantas Kvedaras \email{virmantas.kvedaras@@mif.vu.lt}, Vaidotas Zemlys (maintainer) \email{zemlys@@gmail.com}
##' @keywords package
NULL

##' US monthly unemployment rate
##'
##' The monthly unemployment rate for United States from 1948 to 2011.
##' 
##' @name USunempr
##' @docType data
##' @format A \code{\link{ts}} object.
##' @source \href{http://data.bls.gov/timeseries/LNS14000000}{U.S. Bureau of Labor Statistics}
##' @keywords datasets
NULL


##' US annual gross domestic product in billions of chained 2005 dollars
##'
##' The annual gross domestic product in billions of chained 2005 dollars for US from 1948 to 2011.
##' 
##' @name USrealgdp
##' @docType data
##' @format A \code{\link{ts}} object.
##' @source \href{http://www.bea.gov/national/xls/gdplev.xls}{U.S. Department of Commerce, Bureau of Economic Analysis}
##' @keywords datasets
NULL

##' Realized volatility of S&P500 index
##'
##' Realized volatility of S&P500(Live) index of the period 2000 01 03 - 2013 11 22
##'
##' @name rvsp500
##' @docType data
##' @format A \code{data.frame} object with two columns. First column contains date id, and the second the realized volatility for S&P500 index.
##' @source \href{http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip}{http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip}
##' @references Heber, Gerd and Lunde, Asger, and Shephard, Neil and Sheppard, Kevin \emph{Oxford-Man Institute's realized library}, Oxford-Man Institute, University of Oxford (2009)
##' @keywords datasets
##' @examples
##' ## Do not run:
##' ## Download the data from
##' ## http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip
##' ## It contains the file OxfordManRealizedVolatilityIndices.csv.
##' 
##' ## rvi <- read.csv("OxfordManRealizedVolatilityIndices.csv",check.names=FALSE,skip=2)
##' ## ii <- which(rvi$DateID=="20131112")
##' ## rvsp500 <- na.omit(rvi[1:ii,c("DataID","SPX2.rv")]
NULL

##' Out-of-sample prediction precision data on simulation example
##'
##' The code in the example generates the out-of-sample prediction precision data for correctly and incorrectly constrained MIDAS regression model compared to unconstrained MIDAS regression model. 
##' @name oos_prec
##' @docType data
##' @format A \code{data.frame} object with four columns. The first column indicates the sample size, the second the type of constraint, the third the value of the precision measure and the fourth the type of precision measure.
##' @keywords datasets
##' @examples
##' ## Do not run:
##' ## set.seed(1001)
##' 
##' ## gendata<-function(n) {
##' ##     trend<-c(1:n)
##' ##     z<-rnorm(12*n)
##' ##     fn.z <- nealmon(p=c(2,0.5,-0.1),d=17)
##' ##     y<-2+0.1*trend+mls(z,0:16,12)%*%fn.z+rnorm(n)
##' ##     list(y=as.numeric(y),z=z,trend=trend)
##' ## }
##' 
##' ## nn <- c(50,100,200,300,500,750,1000)
##' 
##' ## data_sets <- lapply(n,gendata)
##' 
##' ## mse <- function(x) {
##' ##     mean(residuals(x)^2)
##' ## }
##' 
##' ## bnorm <- function(x) {
##' ##     sqrt(sum((midas_coef(x)-c(2,0.1,nealmon(p=c(2,0.5,-0.1),d=17)))^2))
##' ## }
##' 
##' ## rep1 <- function(n) {
##' ##     dt <- gendata(round(1.25*n))
##' ##     ni <- n    
##' ##     ind <- 1:ni
##' ##     mind <- 1:(ni*12)
##' ##     indt<-list(y=dt$y[ind],z=dt$z[mind],trend=dt$trend[ind])
##' ##     outdt <- list(y=dt$y[-ind],z=dt$z[-mind],trend=dt$trend[-ind])
##' ##     um <- midas_r(y~trend+mls(z,0:16,12),data=indt,start=NULL)
##' ##     nm <- midas_r(y~trend+mls(z,0:16,12,nealmon),data=indt,start=list(z=c(1,-1,0)))
##' ##     am <- midas_r(y~trend+mls(z,0:16,12,almonp),data=indt,start=list(z=c(1,0,0,0)))
##' ##     modl <- list(um,nm,am)
##' ##     names(modl) <- c("um","nm","am")
##' ##     list(norms=sapply(modl,bnorm),mse=sapply(modl,function(mod)mean((forecast(mod,newdata=outdt)-outdt$y)^2)))
##' ## }
##' 
##' ## repr <- function(n,R) {
##' ##     cc <- lapply(1:R,function(i)rep1(n))
##' ##     list(norms=t(sapply(cc,"[[","norms")),mse=t(sapply(cc,"[[","mse")))
##' ## }
##' 
##' ## res <- lapply(nn,repr,R=1000)
##' 
##' ## norms <- data.frame(nn,t(sapply(lapply(res,"[[","norms"),function(l)apply(l,2,mean))))
##' ## mses <- data.frame(nn,t(sapply(lapply(res,"[[","mse"),function(l)apply(l,2,mean))))
##' 
##' 
##' ## msd <- melt(mses[-1,],id=1)
##' ## colnames(msd)[2] <- "Constraint"
##' ## nmd <- melt(norms[-1,],id=1)
##' ## colnames(nmd)[2] <- "Constraint"
##' 
##' ## msd$Type <- "Mean squared error"
##' ## nmd$Type <- "Distance from true values"
##' ## oos_prec <- rbind(msd,nmd)
##' ## oos_prec$Type <- factor(oos_prec$Type,levels=c("Mean squared error","Distance from true values"))
NULL

##' United States gross domestic product, quarterly, seasonaly adjusted annual rate.
##'
##' United States gross domestic product, quarterly, seasonaly adjusted annual rate. Retrieved from FRED, symbol "GDP" at 2014-04-25.
##' 
##' @name USqgdp
##' @docType data
##' @format A \code{\link{ts}} object.
##' @source \href{http://research.stlouisfed.org/fred2/series/GDP}{FRED, Federal Reserve Economic Data, from the Federal Reserve Bank of St. Louis}
##' @keywords datasets
##' @examples
##' ## Do not run:
##' ## library(quantmod)
##' ## USqgdp <- ts(getSymbols("GDP",src="FRED",auto.assign=FALSE),start=c(1947,1),frequency=4)
NULL

##' United States total employment non-farms payroll, monthly, seasonally adjusted.
##'
##' United States total employment non-farms payroll, monthly, seasonally adjusted. Retrieved from FRED, symbol "PAYEMS" at 2014-04-25.
##' 
##' @name USpayems
##' @docType data
##' @format A \code{\link{ts}} object.
##' @source \href{http://research.stlouisfed.org/fred2/series/PAYEMS}{FRED, Federal Reserve Economic Data, from the Federal Reserve Bank of St. Louis}
##' @keywords datasets
##' @examples
##'
##' ## Do not run:
##' ## library(quantmod)
##' ## USpayems <- ts(getSymbols("PAYEMS",src="FRED",auto.assign=FALSE),start=c(1939,1),frequency=12)
NULL
