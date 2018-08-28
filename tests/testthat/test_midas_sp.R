set.seed(131313)

#Generate the data sets for use in tests
#
n <- 250
nnbeta <- function(p, k) nbeta(c(1,p),k)
dgp_pl <- midas_pl_sim(n, m = 12, theta = nbeta(c(1.5, 2, 4), 24), gfun = function(x)0.25*x^3, ar.x = 0.9, ar.y = 0.5, n.start = 100)
dgp_si <- midas_si_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), gfun = function(x)0.03*x^3, ar.x = 0.9, ar.y = 0.5, n.start = 100)



a100 <- midas_sp(y~mlsd(x, 0:23, y, nbeta)+mlsd(y, 1:2, y) | z, 
                 bws = 0, degree = 1, data = dgp_pl,
                 start = list(x = c(1.5, 2, 4), y = c(0.5, 0)), 
                 method = "Nelder-Mead", control = list(maxit = 100))

b100 <- midas_sp(y~mlsd(y, 1:2, y) | mlsd(x, 0:23, y, nnbeta), 
         bws = 0, degree = 1, data = dgp_si,
         start = list(x = c(2, 4), y = c(0.5, 0)), 
         method = "Nelder-Mead", control = list(maxit = 100))

test_that("Plain and formula interface give the same results for PL", {
    X <- mls(dgp_pl$x, 0:23, 12)
    
    mpl <- midas_pl_plain(dgp_pl$y, X, dgp_pl$z, p.ar = 2L, nbeta, degree = 1, 
                          start_bws = 0, start_x = c(1.5, 2, 4), start_ar = c(0.5, 0),
                          method = "Nelder-Mead", itnmax = 100) 
    
    mfr <- a100

    expect_true(sum(abs(mfr$coefficients - mpl$coefficients)) < 1e-10)
})

test_that("Plain and formula interface give the same results for SI", {
    X <- mls(dgp_si$x, 0:23, 12)
    
   mpl <- midas_si_plain(dgp_si$y, X, p.ar = 2L, nnbeta, degree = 1, start_bws = 0, start_x = c(2, 4), start_ar = c(0.5, 0),
                             method = "Nelder-Mead", itnmax = 100) 
    
   mfr <- midas_sp(y~mlsd(y, 1:2, y) | mlsd(x, 0:23, y, nnbeta), 
                           bws = 0, degree = 1, data = dgp_si,
                           start = list(x = c(2, 4), y = c(0.5, 0)), 
                           method = "Nelder-Mead", control = list(maxit = 100))
        
   cmap <- c(1,4:5, 2:3)
    expect_true(sum(abs(mfr$coefficients[cmap]-mpl$coefficients)) < 1e-10)
})



test_that("Rearanging terms works", {
    
    mfr1 <- a100
    mfr2 <- midas_sp(y~mlsd(y, 1:2, y) + mlsd(x, 0:23, y, nbeta)| z, 
                   bws = 0, degree = 1, data = dgp_pl,
                   start = list(x = c(1.5, 2, 4), y = c(0.5, 0)), 
                   method = "Nelder-Mead", control = list(maxit = 100))
    
    cmap <- c(1,5:6, 2:4)
    
    
    expect_true(sum(abs(mfr1$coefficients[cmap]-mfr2$coefficients)) < 1e-10)
})

test_that("Updating Ofunction works for sp", {
    
    a <- a100
    b <- update(a, Ofunction = "nls")
    c <- suppressWarnings(update(b, Ofunction = "optimx", method = c("BFGS", "spg"), itnmax = 10))
    
    expect_that(a$argmap_opt$Ofunction == "optim", is_true())
    expect_that(b$argmap_opt$Ofunction == "nls", is_true())
    expect_that(c$argmap_opt$Ofunction == "optimx", is_true())
    
    expect_that(inherits(b$opt, "nls"), is_true())
    expect_that(inherits(c$opt, "optimx"), is_true())
    
    expect_that(sum(abs(coef(a) - b$start_opt)) == 0, is_true())
    expect_that(sum(abs(coef(b) - c$start_opt)) == 0, is_true())
    
})

test_that("Updating Ofunction arguments  works", {
    a <- a100 
    b <- update(a, method = "CG", control = list(maxit = 5))
    
    expect_that(b$argmap_opt$method == "CG", is_true())
    
})

test_that("Updating data and starting values works",{
  
    spd <- dgp_si[c("y","x")]
    spd$y <- window(spd$y, start = 1, end = 200)
    spd$x <- window(spd$x, start = c(1,1), end = c(200,12))
    
    a <- b100
    
    b <- update(a, data = spd, start = list(x = c(2, 4), y = c(0.5, 0)), control = list(maxit = 1), method = "Nelder-Mead")
    
    m <- na.omit(cbind(spd$y, mlsd(spd$y,1:2, spd$y), mlsd(spd$x, 0:23, spd$y)))
    
    expect_that(sum(abs(m - b$model)), equals(0))
})
