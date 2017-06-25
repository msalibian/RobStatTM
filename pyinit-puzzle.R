
X <- matrix(rnorm(100),100,1)

X <- cbind(rep(1,100), X)

y <- rnorm(100)

out <- pyinit(X=X, y=y, intercept=FALSE, deltaesc=0.5, cc.scale=1.547, 
              prosac=.5,clean.method="threshold", C.res = 2, 
              prop=.2, py.nit = 20, en.tol=1e-5,
              mscale.rho.fun ="bisquare")

out


