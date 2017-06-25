
library(pyinit)
set.seed(123)
X <- matrix(rnorm(100),100,1)

y <- rnorm(100)

# one covariate, no intercept
# (2p+3 = 5) coefficients are scalars
pyinit(X=X, y=y, intercept=FALSE, deltaesc=0.5, cc.scale=1.547, 
       prosac=.5,clean.method="threshold", C.res = 2, 
       prop=.2, py.nit = 20, en.tol=1e-5, 
       mscale.rho.fun ="bisquare")


# one covariate, add an intercept
# coefficients are still scalars, and there's only ONE candidate?
pyinit(X=X, y=y, intercept=TRUE, deltaesc=0.5, cc.scale=1.547, 
       prosac=.5,clean.method="threshold", C.res = 2, 
       prop=.2, py.nit = 20, en.tol=1e-5, 
       mscale.rho.fun ="bisquare")

# add intercept by hand, and set "intercept = FALSE"
# then it looks OK, but 2p+3 = 8? 
X <- cbind(rep(1,100), X)
pyinit(X=X, y=y, intercept=FALSE, deltaesc=0.5, cc.scale=1.547, 
       prosac=.5,clean.method="threshold", C.res = 2, 
       prop=.2, py.nit = 20, en.tol=1e-5, 
       mscale.rho.fun ="bisquare")



