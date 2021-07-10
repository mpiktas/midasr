## More details about the models can be found in the article
## "The statistical content and empirical testing of the MIDAS restrictions"
## by Virmantas Kvedaras and Vaidotas Zemlys

library(midasr)
data("USunempr")
data("USrealgdp")

y <- diff(log(USrealgdp))
x <- window(diff(USunempr), start = 1949)
trend <- 1:length(y)

allk <- lapply(c(12, 15, 18, 24) - 1, function(k) {
  update(midas_r(y ~ trend + fmls(x, k, 12, nealmon), start = list(x = rep(0, 3))), Ofunction = "nls")
})

#### Compute the derivative test
dtest <- lapply(allk, deriv_tests)

### The first derivative tests, gradient is zero
sapply(dtest, with, first)

### The second derivative tests, hessian is positive definite
sapply(dtest, with, second)

### The minimal eigenvalue of hessian is borderline zero, yet positive.
sapply(dtest, with, min(eigenval))

### Apply hAh test
lapply(allk, hAh_test)

### Apply robust hAh test
lapply(allk, hAhr_test)

### View summaries
lapply(allk, summary)

## Plot the coefficients
dev.new()
par(mfrow = c(2, 2))


plot_info <- lapply(allk, function(x) {
  k <- length(coef(x, midas = TRUE, term = "x"))
  ttl <- sprintf("d = %.0f: p-val.(hAh_HAC) < %.2f", k, max(hAhr_test(x)$p.value, 0.01))
  plot_midas_coef(x, title = ttl)
})
