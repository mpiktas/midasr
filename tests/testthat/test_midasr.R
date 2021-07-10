context("Testing midas_r")
set.seed(1001)

n <- 250
trend <- c(1:n)
x <- rnorm(4 * n)
z <- rnorm(12 * n)
fn_x <- nealmon(p = c(1, -0.5), d = 8)
fn_z <- nealmon(p = c(2, 0.5, -0.1), d = 17)
y <- 2 + 0.1 * trend + mls(x, 0:7, 4) %*% fn_x + mls(z, 0:16, 12) %*% fn_z + rnorm(n)

accuracy <- sqrt(.Machine$double.eps)

data_ts <- list(
  y = ts(as.numeric(y), frequency = 1),
  x = ts(x, frequency = 4),
  z = ts(z, frequency = 12),
  trend = trend
)


test_that("midas_r with start=NULL is the same as lm", {
  eq_u1 <- lm(y ~ trend + mls(x, k = 0:7, m = 4) + mls(z, k = 0:16, m = 12))

  eq_u2 <- midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), start = NULL)

  eq_u3 <- midas_u(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12))

  expect_lt(sum(abs(coef(eq_u1) - coef(eq_u2))), accuracy)
  expect_lt(sum(abs(coef(eq_u3) - coef(eq_u2))), accuracy)
})

test_that("midas_r with mlsd is the same as mls for start=NULL", {
  eq_u1 <- midas_u(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), data = data_ts)

  eq_u2 <- midas_r(y ~ trend + mlsd(x, 0:7, y) + mlsd(z, 0:16, y), start = NULL, data = data_ts)

  eq_u3 <- midas_u(y ~ trend + mlsd(x, 0:7, y) + mlsd(z, 0:16, y), data = data_ts)

  expect_lt(sum(abs(coef(eq_u1) - coef(eq_u2))), accuracy)
  expect_lt(sum(abs(coef(eq_u3) - coef(eq_u2))), accuracy)
})


test_that("midas_r without weights gives the same summary as midas_u", {
  a <- summary(lm(y ~ trend + mls(x, k = 0:7, m = 4) + mls(z, k = 0:16, m = 12)))

  b <- summary(midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), start = NULL), vcov = NULL)

  expect_lt(sum(abs(coef(a) - coef(b))), accuracy)
})

test_that("midas_r without start throws an error", {
  expect_error(midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12)))
})

test_that("midas_u picks up data from main R environment", {
  eq_u1 <- midas_u(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12))

  eq_u2 <- midas_u(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12),
    data = list(y = y, trend = trend, x = x, z = z)
  )

  expect_lt(sum(abs(coef(eq_u1) - coef(eq_u2))), accuracy)
})

test_that("midas_r picks up data from main R environment", {
  eq_u1 <- midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), start = NULL)

  eq_u2 <- midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12),
    data = list(y = y, trend = trend, x = x, z = z), start = NULL
  )

  expect_lt(sum(abs(coef(eq_u1) - coef(eq_u2))), accuracy)
})

test_that("midas_r and midas_u picks data from the parent environment with mlsd", {
  eq_u1 <- midas_u(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), data = data_ts)
  xx <- data_ts$x
  yy <- data_ts$y
  zz <- data_ts$z

  eq_u2 <- midas_r(yy ~ trend + mlsd(xx, 0:7, yy) + mlsd(zz, 0:16, yy), start = NULL)

  eq_u3 <- midas_u(yy ~ trend + mlsd(xx, 0:7, yy) + mlsd(zz, 0:16, yy))

  expect_lt(sum(abs(coef(eq_u1) - coef(eq_u2))), accuracy)
  expect_lt(sum(abs(coef(eq_u3) - coef(eq_u2))), accuracy)
})


test_that("NLS problem solution is close to the DGP", {
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
  expect_lt(sum(abs(coef(a) - c(2, 0.1, c(1, -0.5), c(2, 0.5, -0.1)))), 1)
  expect_lt(sum(abs(coef(a, midas = TRUE) - c(2, 0.1, fn_x, fn_z))), 1)
})

test_that("midas_r gives the same result for mlsd", {
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))

  b <- midas_r(y ~ trend + mlsd(x, 0:7, y, nealmon) + mlsd(z, 0:16, y, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)), data = data_ts)

  expect_lt(sum(abs(coef(a) - coef(b))), 1e-08)
  expect_lt(sum(abs(coef(a, midas = TRUE) - coef(b, midas = TRUE))), accuracy)
})



test_that("Deriv tests give positive results", {
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
  dt <- deriv_tests(a)
  expect_false(dt$first)
  expect_true(dt$second)
  expect_lt(sum(abs(dt$gradient)) / nrow(a$model), (0.002))
})

test_that("Updating Ofunction works", {
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
  b <- update(a, Ofunction = "nls")
  c <- update(b, Ofunction = "optimx", method = c("Nelder-Mead", "BFGS", "spg"))

  expect_true(a$argmap_opt$Ofunction == "optim")
  expect_true(b$argmap_opt$Ofunction == "nls")
  expect_true(c$argmap_opt$Ofunction == "optimx")

  expect_true(inherits(b$opt, "nls"))
  expect_true(inherits(c$opt, "optimx"))

  expect_true(a$convergence == 0)
  expect_true(b$convergence == 0)
  expect_true(c$convergence == 0)
})

test_that("Updating Ofunction works for mlsd", {
  a <- midas_r(y ~ trend + mlsd(x, 0:7, y, nealmon) + mlsd(z, 0:16, y, nealmon),
    start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)), data = data_ts
  )

  b <- update(a, Ofunction = "nls")
  c <- update(b, Ofunction = "optimx", method = c("Nelder-Mead", "BFGS", "spg"))

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
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
  b <- update(a, method = "CG")

  expect_true(b$argmap_opt$method == "CG")
})


test_that("update works with start=NULL", {
  eq_u1 <- midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), start = NULL)

  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
  eq_u2 <- update(a, start = NULL)

  expect_lt(sum(abs(coef(eq_u1) - coef(eq_u2))), accuracy)
})

test_that("updating gradient works", {
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
  b <- update(a, weight_gradients = list(nealmon = nealmon_gradient))

  expect_lt(sum(abs(b$term_info$x$gradient(c(1, 0.1, 0.1)) - nealmon_gradient(c(1, 0.1, 0.1), 8))), accuracy)
})

test_that("updating data and starting values works", {
  dt <- list(y = y, x = x, z = z, trend = trend)
  spd <- split_data(dt, insample = 1:200, outsample = 201:250)

  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon),
    start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
    data = dt
  )

  c <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon),
    start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
    data = spd$indata
  )
  b <- update(a, data = spd$indata, start = list(
    x = c(1, -0.5), z = c(2, 0.5, -0.1),
    `(Intercept)` = c$start_opt["(Intercept)"],
    trend = c$start_opt["trend"]
  ))


  expect_lt(sum(abs(coef(b) - coef(c))), accuracy)
})


test_that("Gradient passing works", {
  eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) +
    mls(z, 0:16, 12, nealmon),
  start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
  weight_gradients = list(nealmon = nealmon_gradient)
  )

  dt <- deriv_tests(eq_r2)
  expect_lt(sum(abs(dt$gradient)) / nrow(eq_r2$model), 1e-3)
})

test_that("Gradient passing works for default gradients", {
  eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) +
    mls(z, 0:16, 12, nealmon),
  start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
  weight_gradients = list()
  )

  expect_lt(sum(abs(eq_r2$term_info$x$gradient(c(1, 0.1, 0.1)) - nealmon_gradient(c(1, 0.1, 0.1), 8))), accuracy)
})

test_that("Gradient passing works with mlsd", {
  eq_r2 <- midas_r(y ~ trend + mlsd(x, 0:7, y, nealmon) +
    mlsd(z, 0:16, y, nealmon),
  start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
  data = data_ts,
  weight_gradients = list(nealmon = nealmon_gradient)
  )

  dt <- deriv_tests(eq_r2)
  expect_lt(sum(abs(dt$gradient)) / nrow(eq_r2$model), 1e-3)
})

test_that("Gradient passing works for default gradients with mlsd", {
  eq_r2 <- midas_r(y ~ trend + mlsd(x, 0:7, y, nealmon) +
    mlsd(z, 0:16, y, nealmon),
  start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
  data = data_ts,
  weight_gradients = list()
  )

  expect_lt(sum(abs(eq_r2$term_info$x$gradient(c(1, 0.1, 0.1)) - nealmon_gradient(c(1, 0.1, 0.1), 8))), accuracy)
})


test_that("Gradient passing works for nls", {
  eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) +
    mls(z, 0:16, 12, nealmon),
  start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
  weight_gradients = list(nealmon = nealmon_gradient), Ofunction = "nls"
  )

  dt <- deriv_tests(eq_r2)
  expect_lt(sum(abs(dt$gradient)) / nrow(eq_r2$model), 1e-3)
})


test_that("Term info gathering works", {
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
  expect_named(coef(a), c("(Intercept)", "trend", "x1", "x2", "z1", "z2", "z3"))
  expect_true(length(coef(a)) == 7)
  expect_true(length(coef(a, midas = TRUE)) == 27)
  expect_named(a$term_info, c("(Intercept)", "trend", "x", "z"))

  lgs <- lapply(a$term_info, "[[", "lag_structure")
  expect_true(lgs[["(Intercept)"]] == 0)
  expect_true(lgs[["trend"]] == 0)
  expect_identical(lgs[["x"]], 0:7)
  expect_identical(lgs[["z"]], 0:16)

  expect_identical(
    a$term_info[["x"]]$weight(coef(a, term_names = "x")),
    nealmon(coef(a, term_names = "x"), 8)
  )

  expect_identical(
    a$term_info[["z"]]$weight(coef(a, term_names = "z")),
    nealmon(coef(a, term_names = "z"), 17)
  )
})


test_that("AR* model works", {
  a <- midas_r(y ~ trend + mls(y, c(1, 4), 1, "*") + mls(x, 0:7, 4, nealmon)
    + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))

  cfx <- coef(a, term_names = "x")
  cfb <- a$term_info$x$weight(cfx)
  mcfx <- coef(a, midas = TRUE, term_names = "x")
  cfy <- coef(a, midas = TRUE, term_names = "y")
  expect_lt(sum(abs(mcfx[1:4] - cfb[1:4])), accuracy)
  expect_lt(sum(abs(c(cfb[5:8], rep(0, 4)) - cfb * cfy[1] - mcfx[5:12])), accuracy)
  expect_lt(sum(abs(-cfb * cfy[2] - mcfx[13:20])), accuracy)
})

test_that("AR* model works with gradient", {
  a <- midas_r(y ~ trend + mls(y, c(1, 4), 1, "*") + mls(x, 0:7, 4, nealmon)
    + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)), weight_gradients = list())
  b <- midas_r(y ~ trend + mls(y, c(1, 4), 1, "*") + mls(x, 0:7, 4, nealmon)
    + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))

  expect_lt(sum(abs(b$gradD(coef(a)) - a$gradD(coef(a)))), accuracy)
})


test_that("Midas_r_plain works", {
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)), Ofunction = "optimx")
  fn_a <- function(p, d) {
    c(nealmon(p[1:2], d = 8), nealmon(p[3:5], d = 17))
  }
  s <- midas_r_plain(y, cbind(mls(x, 0:7, 4), mls(z, 0:16, 12)), cbind(1, trend), fn_a,
    startx = a$start_opt[3:7], startz = a$start_opt[1:2]
  )

  expect_lt(
    sum(abs(coef(s)[c(6:7, 1:5)] - coef(a))),
    accuracy
  )
})

test_that("Midas_r_plain gradient works", {
  a <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)), Ofunction = "optimx", weight_gradients = list())
  fn_a <- function(p, d) {
    c(nealmon(p[1:2], d = 8), nealmon(p[3:5], d = 17))
  }
  gr_fn_a <- function(p, d) {
    gr1 <- nealmon_gradient(p[1:2], d = 8)
    gr2 <- nealmon_gradient(p[3:5], d = 17)
    cbind(
      rbind(gr1, matrix(0, nrow = 17, ncol = 2)),
      rbind(matrix(0, nrow = 8, ncol = 3), gr2)
    )
  }

  s <- midas_r_plain(y, cbind(mls(x, 0:7, 4), mls(z, 0:16, 12)), cbind(1, trend), fn_a,
    startx = a$start_opt[3:7], startz = a$start_opt[1:2],
    grw = gr_fn_a
  )

  expect_lt(
    sum(abs(coef(s)[c(6:7, 1:5)] - coef(a))),
    (1e-11)
  )
  expect_lt(
    sum(abs(s$gradient(coef(s))[c(6:7, 1:5)] - a$gradient(coef(s)[c(6:7, 1:5)]))),
    (1e-10)
  )
  expect_lt(
    sum(abs(s$gradD(coef(s))[c(26:27, 1:25), c(6:7, 1:5)] - a$gradD(coef(s)[c(6:7, 1:5)]))),
    accuracy
  )
})
