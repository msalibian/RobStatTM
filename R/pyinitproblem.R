
library(pyinit)
set.seed(123)
y <- rnorm(100)
X <- matrix(rep(1, 100), 100, 1)
pyinit(X=X, y=y, intercept=FALSE, deltaesc=.5, cc.scale=1.54764,
       prosac=1/100^(1.5), clean.method='threshold', C.res=2, py.nit=10, en.tol=1e-5)

