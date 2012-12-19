##More details about the models can be found in the article
##"The statistical content and empirical testing of the MIDAS restrictions"
##by Virmantas Kvedaras and Vaidotas Zemlys

library(midasr)
data("USunempr")
data("USrealgdp")

y <- diff(log(USrealgdp))
x <- window(diff(USunempr),start=1949)
trend <- 1:length(y)

allk <- foreach(k=c(12,15,18,24)-1) %do% {	
    midas_r(midas_r(y~trend+fmls(x,k,12,nealmon),start=list(x=rep(0,3))),Rfunction="nls")
}
                                                
####Compute the derivative test                
dtest <- lapply(allk,deriv_tests)

###The first derivative tests, gradient is zero
sapply(dtest,with,first)
                
###The second derivative tests, hessian is positive definite
sapply(dtest,with,second)

###The minimal eigenvalue of hessian is borderline zero, yet positive.
sapply(dtest,with,min(eigenval))

###Apply hAh test
lapply(allk,hAh.test)

###Apply robust hAh test
lapply(allk,hAhr.test)

###View summaries
lapply(allk,summary)

##Plot the coefficients
dev.new()
par(mfrow=c(2,2))

lapply(allk,function(x){
    cfur <- coef(x$unrestricted)
    cfur <- cfur[grep("embedlf",names(cfur))]
    cfre <- restr_coef(x)
    k <- length(cfur)
    sdval <- sqrt(diag(vcovHAC(x$unrestricted)))
    sdval <- sdval[grep("embedlf",names(sdval))]

    plot(0:(k-1),cfur,col="black",ylab="Beta coefficients",xlab="h")
    title(main=sprintf("d = %.0f: p-val.(hAh_HAC) < %.2f", k, max(hAhr.test(x)$p.value, 0.01)), cex.main = 1, font.main = 4, col.main = "black")
    points(c(0:(k - 1)), cfre, type = "l", col = "blue")
    points(c(0:(k - 1)), cfur + 2 * sdval[1:k], type = "l", col = "red", lty = 2)
    points(c(0:(k - 1)), cfur - 2 * sdval[1:k], type = "l", col = "red", lty = 2)
})
