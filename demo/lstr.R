set.seed(131313)
n <- 250
burn_in <- 100

nnbeta <- function(p, k) nbeta(c(1,p),k)

dgp_lstr_250 <- midas_lstr_sim(n, m = 12, theta = nnbeta(c(2, 4), 24), linear = c(1, 1.5), lstr = c(log(1), 1, 1), ar.x = 0.9, ar.y = 0.5, n.start = 100)

z <- cbind(1, mls(dgp_lstr_250$y, 1:2, 1))
colnames(z) <- c("Intercept", "y1", "y2")
X <- mls(dgp_lstr_250$x, 0:23, 12)

mm <- na.omit(cbind(dgp_lstr_250$y,X))

lstr_250 <- midas_lstr_simple(dgp_lstr_250$y, X, z, nnbeta, start_lstr = c(1.5, log(1), 1, 1), start_x = c(2, 4), start_z=c(1, 0.5, 0), method = "Nelder-Mead", max.iter = 5000)
