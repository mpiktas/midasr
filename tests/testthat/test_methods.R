context("Testing midas_r methods")
set.seed(1001)
n<-250
trend<-c(1:n)
x<-rnorm(4*n)
z<-rnorm(12*n)
fn_x <- nealmon(p=c(1,-0.5),d=8)
fn_z <- nealmon(p=c(2,0.5,-0.1),d=17)
y<-2+0.1*trend+mls(x,0:7,4)%*%fn_x+mls(z,0:16,12)%*%fn_z+rnorm(n)


spd <- split_data(list(y = y, x = x, trend = trend, z = z), 1:200, 201:250 )

test_that("midas_r preserves ts attribute",{
    spd$indata$y <- ts(spd$indata$y)
    a <- midas_r(y~trend+mls(x,0:7,4,nealmon)+mls(z,0:16,12,nealmon),start=list(x=c(1,-0.5),z=c(2,0.5,-0.1)), data = spd$indata)
    fa <- forecast(a, newdata = spd$outdata)
    expect_true(inherits(fa$mean,"ts"))
    expect_true(identical(time(fa$mean),time(ts(201:250,start=201))))
    
    spd$indata$y <- ts(spd$indata$y, start = 1950, frequency = 4)
    a <- midas_r(y~trend+mls(x,0:7,4,nealmon)+mls(z,0:16,12,nealmon),start=list(x=c(1,-0.5),z=c(2,0.5,-0.1)), data = spd$indata)
    fa <- forecast(a, newdata = spd$outdata)
    expect_true(inherits(fa$mean,"ts"))
    expect_true(identical(time(fa$mean),time(ts(201:250,start=c(2000,1), frequency = 4))))
    
    spd$indata$y <- ts(spd$indata$y, start = 1990, frequency = 12)
    a <- midas_r(y~trend+mls(x,0:7,4,nealmon)+mls(z,0:16,12,nealmon),start=list(x=c(1,-0.5),z=c(2,0.5,-0.1)), data = spd$indata)
    fa <- forecast(a, newdata = spd$outdata)
    expect_true(inherits(fa$mean,"ts"))
    expect_that(sum(abs(time(fa$mean)-time(ts(201:250,start=c(2006, 9), frequency = 12)))), is_less_than(1e-10))
})