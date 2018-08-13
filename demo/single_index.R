set.seed(131313)

nnbeta <- function(p, k) nbeta(c(1,p),k)

do_n_si <- function(n) {
    dgp <- midas_si_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), gfun = function(x)0.03*x^3, ar.x = 0.9, ar.y = 0.5, n.start = 100)
    
    X <- mls(dgp$x, 0:23, 12)
    
    si_mod <- midas_si_plain(dgp$y, X, p.ar = 2L, nnbeta, degree = 1, start_bws = 0, start_x = c(2, 4), start_ar = c(0.5, 0),
                             method = "Nelder-Mead", control = list(maxit = 5000)) 
    list(dgp = dgp, si = si_mod)
}

sim_si <- lapply(c(250,500), do_n_si)