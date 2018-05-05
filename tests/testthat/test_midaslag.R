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

test_that("mlsd works the same as mls", {
    x <- c(1:144)
    y <- c(1:12)
    datex <- x
    datey <- (y-1)*12+1
    m1 <- mlsd(x, 0:5, datex, datey)
    m2 <- mls(x, 0:5, 12)
    
    expect_true(sum(abs(m1-m2)) < 1e-10)
    
})

test_that("mlsd works for the ts objects", {
    x <- ts(c(1:144), freq = 12)
    y <- ts(c(1:48), freq = 4)
    
    m1 <- mlsd(x, 0:7, x, y)
    m2 <- mls(x, 0:7, 3)
    
    expect_true(sum(abs(m1-m2), na.rm = TRUE) < 1e-10)
})

test_that("mlsd works for the xts objects", {
    library(xts)
    library(lubridate)
    data(sample_matrix)
    x <- as.xts(sample_matrix, descr='my new xts object')
    y <- xts(1:6, order.by = unique(floor_date(index(sample.xts), unit = "month")))
    
    m1 <- mlsd(as.numeric(x[,1]), 0:42, x, y)
    
})
