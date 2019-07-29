
library(RobStatTM)
set.seed(123)
x <- 1:100
y <- 2:101 + rnorm(100, sd=.00001)

y <- 2:101
x <- 1:100
a0 <- lmrobM(y~x)
b0 <- lmrobM(x~y)


library(RobStatTM)
y <- 2:101
x <- 1:100
a0 <- lmrobM(y~x)
b0 <- lmrobM(x~y)
a <- lmrobdetMM(y~x)
b <- lmrobdetMM(x~y)
h <- lmrobdetDCML(y~x)
w <- lmrobdetDCML(x~y)


library(RobStatTM)
y <- 101:200
x <- 1:100
a0 <- lmrobM(y~x)
b0 <- lmrobM(x~y)
a <- lmrobdetMM(y~x)
b <- lmrobdetMM(x~y)
h <- lmrobdetDCML(y~x)
w <- lmrobdetDCML(x~y)



summary(abs(y-2:101))

y <- 2:101
x <- 1:100
robust::lmRob(y~x)
robust::lmRob(x~y)

# x <- cbind(rep(1,100), x)
x <- as.matrix(x)
p <- 2; n <- 100
control <- RobStatTM::lmrobdet.control()
pyinit(x=x, y=y, intercept=TRUE,
       cc=1.54764,
       psc_keep=control$psc_keep*(1-(p/n)), resid_keep_method=control$resid_keep_method,
       resid_keep_thresh = control$resid_keep_thresh, resid_keep_prop=control$resid_keep_prop,
       maxit = control$py_maxit, eps=control$py_eps,
       mscale_maxit = control$mscale_maxit, mscale_tol = control$mscale_tol,
       mscale_rho_fun="bisquare")
