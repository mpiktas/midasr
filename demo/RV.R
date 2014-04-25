##More details about the models can be found in the article
##"The statistical content and empirical testing of the MIDAS restrictions"
##by Virmantas Kvedaras and Vaidotas Zemlys

library(midasr)
data(rvsp500)

ii <- which(rvsp500$DateID=="20120522")
y <- as.numeric(rvsp500[1:ii,2])

allh <- lapply(c(5,10,20,40),function(h){
    rvh <- filter(c(rep(0,h),y),c(rep(1,h),rep(0,h+1)))
    rvh <- rvh[-h:-1]
    y <- y[1:length(rvh)]
    midas_r(midas_r(rvh~fmls(y,70-1,1,nealmon),start=list(y=rep(0,3))),Ofunction="nls")
})

####Compute the derivative test                
dtest <- lapply(allh,deriv_tests)

###The first derivative tests, gradient is zero
sapply(dtest,with,first)
                
###The second derivative tests, hessian is positive definite
sapply(dtest,with,second)

###View summaries
lapply(allh,summary)

##Precompute the meat matrix for robust testing. Takes some time to compute!!!
PHI <- lapply(allh,function(x)meatHAC(x$unrestricted,prewhite=TRUE,weights=weightsAndrews))

###Apply hAh test
lapply(allh,hAh.test)

##Apply robust hAh test with precomputed PHI
mapply(hAhr.test,allh,PHI,SIMPLIFY=FALSE)

##Parameter j is superfluous, j=0 means no logarithm transformation was
##applied, j=1 means that logarithm transformation was applied. The graph
##is made to be the same as in the aforementioned article.

graph <- function(x,phi,j,h) {
    cfur <- coef(x$unrestricted)
    cfur <- cfur[grep("fmls",names(cfur))]
    cfre <- weight_coef(x)
    k <- length(cfur)
    sdval <- sqrt(diag(sandwich(x$unrestricted,meat=phi)))
    sdval <- sdval[grep("fmls",names(sdval))]
    pv0hac <- hAhr.test(x,PHI=phi)$p.value
    plot(c(0:(k - 1)), c(cfur), col = "black", ylab = "Beta coefficients", xlab = "h")
    title(main = sprintf("k(H=%.0f,j=%.0f) = %.0f: p-val.(hAh_HAC) < %.2f", h, j, 
        k, max(pv0hac, 0.01)), cex.main = 1, font.main = 4, col.main = "black")
    
    points(c(0:(k - 1)), cfre[1:k], type = "l", col = "blue")
    points(c(0:(k - 1)), cfur[1:k] + 2 * sdval[1:k], type = "l", col = "red", lty = 2)
    points(c(0:(k - 1)), cfur[1:k] - 2 * sdval[1:k], type = "l", col = "red", lty = 2)
}

dev.new()
par(mfrow=c(2,2))

mapply(graph,allh,PHI,as.list(rep(0,4)),as.list(c(5,10,20,40)),SIMPLIFY=FALSE)

