set.seed(131313)

nnbeta <- function(p, k) nbeta(c(1,p),k)

do_n_lstr <- function(n) {
    dgp <- midas_lstr_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), intercept = 1, plstr = c(1.5, 1, 1, 1), ar.x = 0.9, ar.y = 0.5, n.start = 100)

    z <- cbind(1, mls(dgp$y, 1:2, 1))
    colnames(z) <- c("Intercept", "y1", "y2")
    X <- mls(dgp$x, 0:23, 12)

    lstr_mod_plain <- midas_lstr_plain(dgp$y, X, z, nnbeta, start_lstr = c(1.5, 1, 1, 1), start_x = c(2, 4), start_z=c(1, 0.5, 0), 
                                 method = "Nelder-Mead", control= list( maxit= 5000))
    lstr_mod_nlpr <- midas_nlpr(y~mlsd(x, 0:23, y, nnbeta) + mlsd(y,1:2,y), data = dgp,
                                start = list(x = list(lstr = c(1.5, 1, 1, 1),
                                                      r = c(2, 4)), 
                                             y = c(0.5, 0),
                                             `(Intercept)` =  1))
    list(dgp = dgp, plain = lstr_mod_plain, nlpr = lstr_mod_nlpr)
}

sim_lstr <- lapply(c(250,500, 1000, 2000, 5000, 10000, 20000), do_n_lstr)

do_n_mmm <- function(n) {
    dgp <- midas_mmm_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), intercept = 1, pmmm = c(1.5, 1), ar.x = 0.9, ar.y = 0.5, n.start = 100)
    
    z <- cbind(1, mls(dgp$y, 1:2, 1))
    colnames(z) <- c("Intercept", "y1", "y2")
    X <- mls(dgp$x, 0:23, 12)
    
    mmm_plain <- midas_mmm_plain(dgp$y, X, z, nnbeta, start_mmm = c(1.5, 1), start_x = c(2, 4), start_z=c(1, 0.5, 0),
                           method = "Nelder-Mead", control = list(maxit = 5000)) 
    
    mmm_nlpr <- midas_nlpr(y~mlsd(x, 0:23, y, nnbeta) + mlsd(y,1:2,y), data = dgp,
               start = list(x = list(mmm = c(1.5, 1),
                                     r = c(2, 4)), 
                            y = c(0.5, 0),
                            `(Intercept)` =  1), 
               method = "Nelder-Mead", control = list(maxit = 5000))
    
    list(dgp = dgp, plain = mmm_plain, nlpr = mmm_nlpr)
}

sim_mmm <- lapply(c(250,500, 1000, 2000, 5000, 10000, 20000), do_n_mmm)






