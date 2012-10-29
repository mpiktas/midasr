test_that("Embedding to low frequency works as expected",{
    x <- 1:24
    x1 <- matrix(2*(1:12),ncol=1)
    x2 <- cbind(2*(1:12),2*(1:12)-1)
    x3 <- rbind(c(NA,NA,NA),cbind(2*(2:12),2*(2:12)-1,2*(2:12)-2))
    x4 <- rbind(c(NA,NA),embed(x,2))
    expect_equivalent(embedlf.default(x,1,2),x1)
    expect_equivalent(embedlf.default(x,2,2),x2)
    expect_equivalent(embedlf.default(x,3,2),x3)
    expect_equivalent(embedlf.default(x,2,1),x4)
})
