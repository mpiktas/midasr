## More details about the models can be found in the article
## "The statistical content and empirical testing of the MIDAS restrictions"
## by Virmantas Kvedaras and Vaidotas Zemlys

library(midasr)
data(rvsp500)

ii <- which(rvsp500$DateID == "20120522")
y <- log(as.numeric(rvsp500[1:ii, 2]))

nlmn <- function(p, d, m) {
  i <- 1:d / 100
  plc <- poly(i, degree = length(p) - 1, raw = TRUE) %*% p[-1]
  as.vector(p[1] * exp(plc) / sum(exp(plc)))
}

## Convergence is harder to achieve, so try to pick better starting values
prestart <- function(start, cfur, k) {
  fn0 <- function(p) {
    sum((cfur - nlmn(p, k))^2)
  }
  optim(start, fn0, method = "BFGS", control = list(maxit = 1000))$par
}

rvhk <- function(h, k) {
  rvh <- filter(c(rep(0, h), y), c(rep(1, h), rep(0, h + 1)))
  rvh <- rvh[-h:-1]
  y <- y[1:length(rvh)]
  mu <- midas_u(rvh ~ fmls(y, k, 1))
  cfur <- coef(mu)[grep("fmls", names(coef(mu)))]
  update(midas_r(rvh ~ fmls(y, k, 1, nlmn), start = list(y = prestart(c(0.2, -1, 1), cfur, k + 1))), Ofunction = "nls")
}

allh <- lapply(c(5, 10, 20, 40), rvhk, k = 69)

#### Compute the derivative test
dtest <- lapply(allh, deriv_tests, tol = 0.5)

### The first derivative tests, gradient is zero
sapply(dtest, with, first)

### The second derivative tests, hessian is positive definite
sapply(dtest, with, second)

### View summaries
lapply(allh, summary)

## Precompute the meat matrix for robust testing. Takes some time to compute!!!
PHI <- lapply(allh, function(x) meatHAC(x$unrestricted, prewhite = TRUE, weights = weightsAndrews))

### Apply hAh test
lapply(allh, hAh_test)

## Apply robust hAh test with precomputed PHI
mapply(hAhr_test, allh, PHI, SIMPLIFY = FALSE)

## Parameter j is superfluous, j=0 means no logarithm transformation was
## applied, j=1 means that logarithm transformation was applied. The graph
## is made to be the same as in the aforementioned article.

graph <- function(x, phi, j, h) {
  k <- length(coef(x, midas = TRUE, term_names = "y"))
  pv0hac <- hAhr_test(x, PHI = phi)$p.value
  ttl <- sprintf("k(H=%.0f,j=%.0f) = %.0f: p-val.(hAh_HAC) < %.2f", h, j, k, max(pv0hac, 0.01))
  plot_midas_coef(x, title = ttl, term_name = "y")
}

dev.new()
par(mfrow = c(2, 2))

plot_info <- mapply(graph, allh, PHI, as.list(rep(1, 4)), as.list(c(5, 10, 20, 40)), SIMPLIFY = FALSE)
