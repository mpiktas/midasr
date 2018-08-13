set.seed(131313)

nnbeta <- function(p, k) nbeta(c(1,p),k)

do_n_si <- function(n) {
    dgp <- midas_si_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), gfun = function(x)0.03*x^3, ar.x = 0.9, ar.y = 0.5, n.start = 100)
    
    X <- mls(dgp$x, 0:23, 12)
    y <- as.numeric(dgp$y)
    si_mod <- midas_si_plain(y, X, p.ar = 2L, nnbeta, degree = 1, start_bws = 0, start_x = c(2, 4), start_ar = c(0.5, 0),
                             method = "Nelder-Mead", itnmax = 5000) 
    list(dgp = dgp, si = si_mod)
}

sim_si <- lapply(c(250,500), do_n_si)

dgp <- sim_si[[1]]$dgp

X <- mls(dgp$x, 0:23, 12)
y <- as.numeric(dgp$y)
source("~/Downloads/R-kodai-duom/np-SI-midas-1.R")

b.gam <- c(2,4) 
m <- 12
k <- 24
f.r <- function(p){nbeta(p = c(1,p,0), d = k)}
sid.250 <- f.si(y=y,xx=X,b.gam,b.ar=c(0.5,0),f.r,h.g=1,pln.lpsn=1,mx.iter=5000)
system.time(si_mod <- midas_si_plain(y, X, p.ar = 2L, nnbeta, degree = 1, start_bws = 0, start_x = c(2, 4), start_ar = c(0.5, 0),
                                     method = "Nelder-Mead", itnmax = 5000)) 
