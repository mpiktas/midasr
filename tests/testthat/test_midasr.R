set.seed(1001)
n<-250
trend<-c(1:n)
x<-rnorm(4*n)
z<-rnorm(12*n)
fn_x <- nealmon(p=c(1,-0.5),d=8)
fn_z <- nealmon(p=c(2,0.5,-0.1),d=17)
y<-2+0.1*trend+mls(x,0:7,4)%*%fn_x+mls(z,0:16,12)%*%fn_z+rnorm(n)


test_that("midas_r with start=NULL is the same as lm",{
    eq_u1<-lm(y~trend+mls(x,k=0:7,m=4)+mls(z,k=0:16,m=12))
    
    eq_u2<-midas_r(y~trend+mls(x,0:7,4)+mls(z,0:16,12),start=NULL)
    
    eq_u3<-midas_u(y~trend+mls(x,0:7,4)+mls(z,0:16,12))
    
    expect_that(sum(abs(coef(eq_u1) - coef(eq_u2))), equals(0))
    expect_that(sum(abs(coef(eq_u3) - coef(eq_u2))), equals(0))
})

test_that("midas_r without weights gives the same summary as",{
    a <-summary(lm(y~trend+mls(x,k=0:7,m=4)+mls(z,k=0:16,m=12)))
    
    b <-summary(midas_r(y~trend+mls(x,0:7,4)+mls(z,0:16,12),start=NULL), vcov = NULL)
    
    expect_that(sum(abs(coef(a) - coef(b))), equals(0))
    
})

test_that("midas_r without start throws an error",{
    expect_that(midas_r(y~trend+mls(x,0:7,4)+mls(z,0:16,12)), throws_error())
    
})


test_that("midas_r picks up data from main R environment",{
    eq_u1<-midas_r(y~trend+mls(x,0:7,4)+mls(z,0:16,12),start=NULL)
    
    eq_u2<-midas_r(y~trend+mls(x,0:7,4)+mls(z,0:16,12), 
                   data = list(y = y, trend = trend, x = x, z = z), start = NULL)
    
    expect_that(sum(abs(coef(eq_u1) - coef(eq_u2))), equals(0))
    
})

test_that("NLS problem solution is close to the DGP", {
    a <- midas_r(y~trend+mls(x,0:7,4,nealmon)+mls(z,0:16,12,nealmon),start=list(x=c(1,-0.5),z=c(2,0.5,-0.1)))
    expect_that(sum(abs(coef(eq.r) - c(2,0.1,c(1,-0.5),c(2,0.5,-0.1)))), is_less_than(1))
    expect_that(sum(abs(coef(a,midas=TRUE)-c(2,0.1,fn_x,fn_z))), is_less_than(1))
})

test_that("Deriv tests give positive results", {
    a <- midas_r(y~trend+mls(x,0:7,4,nealmon)+mls(z,0:16,12,nealmon),start=list(x=c(1,-0.5),z=c(2,0.5,-0.1)))
    dt<-deriv_tests(a)
    expect_that(dt$first, is_false())
    expect_that(dt$second, is_true())
    expect_that(sum(abs(dt$gradient))/nrow(a$model), is_less_than(0.002))
})

test_that("Updating Ofunction works", {
    a <- midas_r(y~trend+mls(x,0:7,4,nealmon)+mls(z,0:16,12,nealmon),start=list(x=c(1,-0.5),z=c(2,0.5,-0.1)))
    b <- update(a, Ofunction = "nls")
    c <- update(b, Ofunction = "optimx", method = c("Nelder-Mead", "BFGS", "spg"))
    
    expect_that(a$argmap_opt$Ofunction == "optim", is_true())
    expect_that(b$argmap_opt$Ofunction == "nls", is_true())
    expect_that(c$argmap_opt$Ofunction == "optimx", is_true())
    
    expect_that(inherits(b$opt, "nls"), is_true())
    expect_that(inherits(c$opt, "optimx"), is_true())
    
    expect_that(a$convergence == 0, is_true())
    expect_that(b$convergence == 0, is_true())
    expect_that(c$convergence == 0, is_true())
    
})

test_that("Gradient passing works", {
    eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + 
                         mls(z, 0:16, 12, nealmon),
                     start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)),
                     weight_gradients=list(nealmon = nealmon_gradient))
    
    dt <- deriv_tests(eq_r2)
    expect_that(sum(abs(dt$gradient))/nrow(a$model), is_less_than(1e-3))
    })

test_that("Term info gathering works", {
    a <- midas_r(y~trend+mls(x,0:7,4,nealmon)+mls(z,0:16,12,nealmon),start=list(x=c(1,-0.5),z=c(2,0.5,-0.1)))
    expect_that(names(coef(a)), is_equivalent_to(c("(Intercept)", "trend", "x1", "x2", "z1", "z2", "z3")))         
    expect_that(length(coef(a)), equals(7))
    expect_that(length(coef(a, midas = TRUE)), equals(27))
    expect_that(names(a$term_info), is_equivalent_to(c("(Intercept)", "trend", "x", "z")))
    
    lgs <- lapply(a$term_info,"[[","lag_structure")
    expect_that(lgs[["(Intercept)"]], is_equivalent_to(0))    
    expect_that(lgs[["trend"]], is_equivalent_to(0))
    expect_that(lgs[["x"]], is_equivalent_to(0:7))    
    expect_that(lgs[["z"]], is_equivalent_to(0:16))
    
    expect_that(a$term_info[["x"]]$weight(coef(a, term_names = "x")), 
                is_equivalent_to(nealmon(coef(a, term_names = "x"), 8)))
    
    expect_that(a$term_info[["z"]]$weight(coef(a, term_names = "z")), 
                is_equivalent_to(nealmon(coef(a, term_names = "z"), 17)))
    
})
