set.seed(131313)

#Generate the data sets for use in tests
#
n <- 250
nnbeta <- function(p, k) nbeta(c(1,p),k)
dgp_pl <- midas_pl_sim(n, m = 12, theta = nbeta(c(1.5, 2, 4), 24), gfun = function(x)0.25*x^3, ar.x = 0.9, ar.y = 0.5, n.start = 100)
dgp_si <- midas_si_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), gfun = function(x)0.03*x^3, ar.x = 0.9, ar.y = 0.5, n.start = 100)

a100 <- midas_sp(y~mlsd(x, 0:23, y, nbeta)+mlsd(y, 1:2, y) | z, 
                 bws = 1, degree = 1, data = dgp_pl,
                 start = list(x = c(1.5, 2, 4), y = c(0.5, 0)), 
                 method = "Nelder-Mead", control = list(maxit = 100))

b100 <- midas_sp(y~mlsd(y, 1:2, y) | mlsd(x, 0:23, y, nnbeta), 
         bws = 1, degree = 1, data = dgp_si,
         start = list(x = c(2, 4), y = c(0.5, 0)), 
         method = "Nelder-Mead", control = list(maxit = 100))
accuracy <- sqrt(.Machine$double.eps)

test_that("Plain and formula interface give the same results for PL", {
    skip_on_cran()  
    X <- mls(dgp_pl$x, 0:23, 12)
    
    mpl <- midas_pl_plain(dgp_pl$y, X, dgp_pl$z, p.ar = 2L, nbeta, degree = 1, 
                          start_bws = 1, start_x = c(1.5, 2, 4), start_ar = c(0.5, 0),
                          method = "Nelder-Mead", control = list(maxit = 100)) 
    
    mfr <- a100
    
    fmpl <- fitted(mpl)
    fmfr <- fitted(mfr)
    
    expect_true((sum(abs(mfr$coefficients - mpl$coefficients)) < accuracy) & (sum(abs(fmpl-fmfr)) < accuracy))
})

test_that("Plain and formula interface give the same results for SI", {
    skip_on_cran()
    X <- mls(dgp_si$x, 0:23, 12)
    
   mpl <- midas_si_plain(dgp_si$y, X, p.ar = 2L, nnbeta, degree = 1, start_bws = 1, start_x = c(2, 4), start_ar = c(0.5, 0),
                             method = "Nelder-Mead", control = list(maxit = 100)) 
    
   mfr <- midas_sp(y~mlsd(y, 1:2, y) | mlsd(x, 0:23, y, nnbeta), 
                           bws = 1, degree = 1, data = dgp_si,
                           start = list(x = c(2, 4), y = c(0.5, 0)), 
                           method = "Nelder-Mead", control = list(maxit = 100))
        
   cmap <- c(1,4:5, 2:3)
   
   fmpl <- fitted(mpl)
   fmfr <- fitted(mfr)
   
   expect_true((sum(abs(mfr$coefficients[cmap] - mpl$coefficients)) < accuracy)  & (sum(abs(fmpl-fmfr)) < accuracy))
})


test_that("Rearanging terms works", {
    skip_on_cran()
    mfr1 <- a100
    mfr2 <- midas_sp(y~mlsd(y, 1:2, y) + mlsd(x, 0:23, y, nbeta)| z, 
                   bws = 1, degree = 1, data = dgp_pl,
                   start = list(x = c(1.5, 2, 4), y = c(0.5, 0)), 
                   method = "Nelder-Mead", control = list(maxit = 100))
    
    cmap <- c(1,5:6, 2:4)
    
    
    expect_true(sum(abs(mfr1$coefficients[cmap]-mfr2$coefficients)) < accuracy)
})

test_that("Updating Ofunction works for sp", {
    skip_on_cran()
    a <- a100
    b <- update(a, Ofunction = "nls")
    c <- suppressWarnings(update(b, Ofunction = "optimx", method = c("BFGS", "spg"), itnmax = 10))
    
    expect_true(a$argmap_opt$Ofunction == "optim")
    expect_true(b$argmap_opt$Ofunction == "nls")
    expect_true(c$argmap_opt$Ofunction == "optimx")
    
    expect_true(inherits(b$opt, "nls"))
    expect_true(inherits(c$opt, "optimx"))
  
    expect_true(sum(abs(coef(a) - b$start_opt)) == 0)
    expect_true(sum(abs(coef(b) - c$start_opt)) == 0)
    
})

test_that("Updating Ofunction arguments  works", {
    skip_on_cran()
    a <- a100 
    b <- update(a, method = "CG", control = list(maxit = 5))
    
    expect_true(b$argmap_opt$method == "CG")
    
})

test_that("Updating data and starting values works",{
    skip_on_cran()
    spd <- dgp_si[c("y","x")]
    spd$y <- window(spd$y, start = 1, end = 200)
    spd$x <- window(spd$x, start = c(1,1), end = c(200,12))
    
    a <- b100
    
    b <- update(a, data = spd, start = list(x = c(2, 4), y = c(0.5, 0)), control = list(maxit = 1), method = "Nelder-Mead")
    
    m <- na.omit(cbind(spd$y, mlsd(spd$y,1:2, spd$y), mlsd(spd$x, 0:23, spd$y)))
    
    expect_lt(sum(abs(m - b$model)), accuracy)
})

test_that("Predicting works for PL", {
    skip_on_cran()
    r <- predict(a100, newdata = dgp_pl) - fitted(a100)
    expect_equivalent(sum(abs(r)), 0)
})



test_that("Predicting works for SI", {
    skip_on_cran()
    r <- predict(b100, newdata = dgp_si) - fitted(b100)
    expect_that(sum(abs(r)), equals(0))
})


test_that("Predicting works for pure SI model", {
    skip_on_cran()   
   bb <- midas_sp(y~ mlsd(x, 0:23, y, nnbeta), 
                     bws = 1, degree = 1, data = dgp_si,
                     start = list(x = c(2, 4)), 
                     method = "Nelder-Mead", control = list(maxit = 100))
    r <- predict(bb, newdata = dgp_si) - fitted(bb)
    expect_that(sum(abs(r)), equals(0))
}) 


test_that("g_np and g_np_mv works the same for numeric and for matrices in case of a vector", {
    skip_on_cran()
    oo <- midasr:::midas_sp_fit(a100)
    
    rn <- g_np(as.numeric(oo$y - oo$xi), as.numeric(oo$z), as.numeric(oo$z), coef(a100)[1], a100$degree)
    rm <- g_np_mv(oo$y - oo$xi, oo$z, oo$z, coef(a100)[1], a100$degree)
    expect_that(sum(abs(rn-rm)), equals(0))
})


test_that("cv_np works the same for numeric and for matrices in case of a vector", {
    oo <- midasr:::midas_sp_fit(a100)
    
    rn <- cv_np(as.numeric(oo$y - oo$xi), as.numeric(oo$z),coef(a100)[1], a100$degree)
    rm <- cv_np(oo$y - oo$xi, oo$z, coef(a100)[1], a100$degree)
    expect_that(sum(abs(rn - rm)), equals(0))
})
