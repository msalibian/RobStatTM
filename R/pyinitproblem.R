
library(pyinit)
set.seed(123)
y <- rnorm(100)

X <- matrix(rep(1, 100), 100, 1)

pyinit(X=X, y=y, intercept=FALSE, deltaesc=.5, cc.scale=1.54764,
       prosac=.5, clean.method='threshold', C.res=2, py.nit=10, en.tol=1e-5)$objF

pyinit(X=X, y=y, intercept=FALSE, deltaesc=.5, cc.scale=1.54764,
       prosac=1/100-1e-6, clean.method='threshold', C.res=2, py.nit=10, en.tol=1e-5)


X <- matrix(rnorm(100, mean=1, sd=.00001), 100, 1)

