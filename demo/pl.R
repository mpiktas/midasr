set.seed(131313)


do_n_pl <- function(n) {
    dgp <- midas_pl_sim(n, m = 12, theta = nbeta(c(1.5, 2, 4), 24), gfun = function(x)0.25*x^3, ar.x = 0.9, ar.y = 0.5, n.start = 100)
    
    X <- mls(dgp$x, 0:23, 12)
    
    pl_mod <- midas_pl_plain(dgp$y, X, dgp$z, p.ar = 2L, nbeta, degree = 1, start_bws = 0, start_x = c(1.5, 2, 4), start_ar = c(0.5, 0),
                             method = "Nelder-Mead", itnmax = 5000) 
    list(dgp = dgp, pl = pl_mod)
}

sim_pl <- lapply(c(250,500), do_n_pl)


