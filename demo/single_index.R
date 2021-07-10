set.seed(131313)

nnbeta <- function(p, k) nbeta(c(1, p), k)

do_n_si <- function(n) {
  dgp <- midas_si_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), gfun = function(x) 0.03 * x^3, ar.x = 0.9, ar.y = 0.5, n.start = 100)

  X <- mls(dgp$x, 0:23, 12)
  y <- as.numeric(dgp$y)
  si_mod <- midas_si_plain(y, X,
    p.ar = 2L, nnbeta, degree = 1, start_bws = 0, start_x = c(2, 4), start_ar = c(0.5, 0),
    method = "Nelder-Mead", itnmax = 5000
  )
  list(dgp = dgp, si = si_mod)
}

t1 <- system.time(sim_si <- lapply(c(250, 500), do_n_si))

t2 <- system.time(msi <- midas_sp(y ~ mlsd(y, 1:2, y) | mlsd(x, 0:23, y, nnbeta),
  bws = 0, degree = 1, data = sim_si[[1]]$dgp,
  start = list(x = c(2, 4), y = c(0.5, 0)), method = "Nelder-Mead", control = list(maxit = 5000)
))
