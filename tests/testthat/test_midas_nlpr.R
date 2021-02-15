set.seed(131313)

#Generate the data sets for use in tests
#
nnbeta <- function(p, k) nbeta(c(1,p),k)

dgp_lstr <- midas_lstr_sim(250, m = 12, theta = nnbeta(c(2, 4), 24), intercept = 1, plstr = c(1.5, 1, 1, 1), ar.x = 0.9, ar.y = 0.5, n.start = 100)
dgp_mmm <-  midas_mmm_sim(250, m = 12, theta = nnbeta(c(2, 4), 24), intercept = 1, pmmm = c(1.5, 1), ar.x = 0.9, ar.y = 0.5, n.start = 100)

accuracy <- sqrt(.Machine$double.eps)

test_that("Plain and formula interface give the same results for LSTR", {
    z <- cbind(1, mls(dgp_lstr$y, 1:2, 1))
    colnames(z) <- c("Intercept", "y1", "y2")
    X <- mls(dgp_lstr$x, 0:23, 12)
    
    mpl <- midas_lstr_plain(dgp_lstr$y, X, z, nnbeta, start_lstr = c(1.5, 1, 1, 1), start_x = c(2, 4), start_z = c(1, 0.5, 0), method = "Nelder-Mead", control = list(maxit = 100)) 
    mfr <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_lstr, 
                              start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                           y = c(0.5, 0),
                                           `(Intercept)` = 1), Ofunction = "optim", method = "Nelder-Mead", control = list(maxit = 100))
    cmap <- c(6:9, 4:5, 1:3)
    
    expect_true(sum(abs(mfr$coefficients[cmap] - mpl$coefficients)) < 1e-10)
})

test_that("Plain and formula interface give the same results for MMM", {
    z <- cbind(1, mls(dgp_mmm$y, 1:2, 1))
    colnames(z) <- c("Intercept", "y1", "y2")
    X <- mls(dgp_mmm$x, 0:23, 12)
    
    mpl <- midas_mmm_plain(dgp_mmm$y, X, z, nnbeta, start_mmm = c(1.5, 1), start_x = c(2, 4), start_z = c(1, 0.5, 0), method = "Nelder-Mead") 
    mfr <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_mmm, 
                      start = list(x = list(r = c(2, 4), mmm = c(1.5, 1)),
                                   y = c(0.5, 0),
                                   `(Intercept)` = 1), Ofunction = "optimx", method = "Nelder-Mead")
    
    cmap <- c(6:7, 4:5, 1:3)
    
    expect_true(sum(abs(mfr$coefficients[cmap] - mpl$coefficients)) < 1e-10)
})

test_that("Rearanging terms works", {
    mfr1 <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_lstr, 
                      start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                   y = c(0.5, 0),
                                   `(Intercept)` = 1), Ofunction = "optim", method = "Nelder-Mead", control = list(maxit = 100))
   
    mfr2 <- midas_nlpr(y~ mlsd(x, 0:23, y, nnbeta) + mlsd(y, 1:2,  y), data = dgp_lstr, 
                       start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                    y = c(0.5, 0),
                                    `(Intercept)` = 1), Ofunction = "optim", method = "Nelder-Mead", control = list(maxit = 100))
    cmap <- c(1, 4:9, 2:3)
    
    expect_true(sum(abs(mfr1$coefficients[cmap] - mfr2$coefficients)) < 1e-10)
})

test_that("midas_nlpr works with the intercept term only", {
    mfr1 <- midas_nlpr(y~ mlsd(x, 0:23, y, nnbeta), data = dgp_lstr, 
                       start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                    `(Intercept)` = 1), Ofunction = "optimx", method = "Nelder-Mead", control = list(maxit = 10))
    
    
    expect_true(length(coef(mfr1)) == 7)
})


test_that("Updating Ofunction works for nplr", {
    a <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_lstr, 
                       start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                    y = c(0.5, 0),
                                    `(Intercept)` = 1), Ofunction = "optim", method = "Nelder-Mead", control = list(maxit = 5000))
    
    b <- update(a, Ofunction = "nls")
    c <- suppressWarnings(update(b, Ofunction = "optimx", method = c("Nelder-Mead", "BFGS", "spg")))
    
    expect_true(a$argmap_opt$Ofunction == "optim")
    expect_true(b$argmap_opt$Ofunction == "nls")
    expect_true(c$argmap_opt$Ofunction == "optimx")
    
    expect_true(inherits(b$opt, "nls"))
    expect_true(inherits(c$opt, "optimx"))
    
    expect_true(a$convergence == 0)
    expect_true(b$convergence == 0)
    expect_true(c$convergence == 0)
    
})

test_that("Updating Ofunction arguments  works", {
    a <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_lstr, 
                    start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                 y = c(0.5, 0),
                                 `(Intercept)` = 1), Ofunction = "optim", method = "Nelder-Mead", control = list(maxit = 5000))
    
    
    b <- update(a, method = "CG")
    
    expect_true(b$argmap_opt$method == "CG")
    
})

test_that("updating data and starting values works",{
  
    spd <- dgp_lstr[c("y","x")]
    spd$y <- window(spd$y, start = 1, end = 200)
    spd$x <- window(spd$x, start = c(1,1), end = c(200,12))
    
    a <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_lstr, 
                    start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                 y = c(0.5, 0),
                                 `(Intercept)` = 1), Ofunction = "optim", method = "Nelder-Mead", control = list(maxit = 5000))
    
    
    c <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = spd, 
                    start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                 y = c(0.5, 0),
                                 `(Intercept)` = 1), Ofunction = "optim", method = "Nelder-Mead", control = list(maxit = 5000))
    
    
    b <- update(a, data = spd, start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                                     y = c(0.5, 0),
                                                     `(Intercept)` = 1), Ofunction = "optim", method = "Nelder-Mead", control = list(maxit = 5000)) 
    
    
    expect_lt(sum(abs(coef(b) - coef(c))), accuracy)
})

test_that("LSTR standard errors work", {
    mfr <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_lstr, 
                      start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                   y = c(0.5, 0),
                                   `(Intercept)` = 1), Ofunction = "optimx", method = "Nelder-Mead", control = list(maxit = 5000))
    
    mfr1 <- update(mfr, Ofunction = "nls")
    
    s1 <- summary(mfr1)
    s2 <- summary(mfr1$opt)
    
    expect_true(sum(abs(s2$coefficients[, 2] - s1$coefficients[, 2])) < sqrt(accuracy))
  
})

test_that("MMM standard errors work", {
    
    mfr <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_mmm, 
                      start = list(x = list(r = c(2, 4), mmm = c(1.5, 1)),
                                   y = c(0.5, 0),
                                   `(Intercept)` = 1), Ofunction = "optimx", method = "Nelder-Mead", control = list(maxit = 5000))
    
    mfr1 <- update(mfr, Ofunction = "nls")
    
    s1 <- summary(mfr1)
    s2 <- summary(mfr1$opt)
    
    expect_true(sum(abs(s2$coefficients[, 2] - s1$coefficients[, 2])) < sqrt(accuracy))
    
})

test_that("Predicting works for LSTR", {
    
    mfr <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_lstr, 
                      start = list(x = list(r = c(2, 4), lstr = c(1.5, 1, 1, 1)),
                                   y = c(0.5, 0),
                                   `(Intercept)` = 1), Ofunction = "optimx", method = "Nelder-Mead", control = list(maxit = 5000))
    
    p1 <- predict(mfr, newdata = list(x = window(dgp_lstr$x, end = c(200,12)), y = window(dgp_lstr$y, end = 200))) 
    
    p2 <- predict(mfr)[1:198]
    
    expect_true(sum(abs(p1 - p2)) < accuracy)
})

test_that("Predicting works for MMM", {
    mfr <- midas_nlpr(y~mlsd(y, 1:2,  y) + mlsd(x, 0:23, y, nnbeta), data = dgp_mmm, 
                      start = list(x = list(r = c(2, 4), mmm = c(1.5, 1)),
                                   y = c(0.5, 0),
                                   `(Intercept)` = 1), Ofunction = "optimx", method = "Nelder-Mead", control = list(maxit = 5000))
    
    
    p1 <- predict(mfr, newdata = list(x = window(dgp_mmm$x, end = c(200,12)), y = window(dgp_mmm$y, end = 200))) 
    
    p2 <- predict(mfr)[1:198]
    
    expect_true(sum(abs(p1 - p2)) < accuracy)
})

