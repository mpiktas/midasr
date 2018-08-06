set.seed(131313)

nnbeta <- function(p, k) nbeta(c(1,p),k)

do_n <- function(n) {
    dgp <- midas_lstr_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), intercept = 1, lstr = c(1.5, 1, log(1), 1), ar.x = 0.9, ar.y = 0.5, n.start = 100)

    z <- cbind(1, mls(dgp$y, 1:2, 1))
    colnames(z) <- c("Intercept", "y1", "y2")
    X <- mls(dgp$x, 0:23, 12)

    lstr <- midas_lstr_plain(dgp$y, X, z, nnbeta, start_lstr = c(1.5, 1, log(1), 1), start_x = c(2, 4), start_z=c(1, 0.5, 0)) 
    list(dgp = dgp, lstr = lstr)
}

sim_lstr <- lapply(c(250,500), do_n)