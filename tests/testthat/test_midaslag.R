context("Testing mls")
test_that("Embedding to low frequency works as expected",{
    x <- 1:24
    x1 <- matrix(2*(1:12),ncol=1)
    x2 <- cbind(2*(1:12),2*(1:12)-1)
    x3 <- rbind(c(NA,NA,NA),cbind(2*(2:12),2*(2:12)-1,2*(2:12)-2))
    x4 <- rbind(c(NA,NA),embed(x,2))
    expect_equivalent(fmls(x,0,2),x1)
    expect_equivalent(fmls(x,1,2),x2)
    expect_equivalent(fmls(x,2,2),x3)
    expect_equivalent(fmls(x,1,1),x4)
})

test_that("mls and fmls give the same results",{
    set.seed(123)
    e <- rnorm(100)
    a <- mls(e, 0:10, 10)
    b <- fmls(e, 10, 10)
    expect_that(sum(abs(a - b), na.rm = TRUE), equals(0))
})

test_that("fmls and dmls give the same results",{
    set.seed(123)
    e <- rnorm(100)
    a <- fmls(c(NA, diff(e)), 10, 10)
    b <- dmls(e, 10, 10)
    expect_that(sum(abs(a - b), na.rm = TRUE), equals(0))
})
