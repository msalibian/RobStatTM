
# problem is prosac is too high
# take prosac < (p / n) / 4 ?
rm(list=ls())
library(pyinit)
data(coleman, package='robustbase')
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

x <- model.matrix(Y ~ .  , data=co2)
y <- co2$Y

n <- nrow(x)
p <- ncol(x)
dee <- .5*(1-(p/n))

control <- list(seed = NULL, tuning.chi = 1.5477, bb = 0.5, # 50% Breakdown point
                            tuning.psi = 3.4434, # 85% efficiency
                            max.it = 500, refine.tol = 1e-7, rel.tol = 1e-7,
                            solve.tol = 1e-7, trace.lev = 0, mts = 1000,
                            compute.rd = FALSE, psi = 'bisquare',
                            split.type = "f",
                            cov = FALSE, initial='S', method='MM', subsampling='simple',
                            candidates = 'PY', 
                            # pyinit control 
                            prosac = 0.1, clean.method = 'threshold', 
                            C.res = 2, prop = .2, py.nit = 20, en.tol = 1e-5, 
                            mscale.maxit = 200, mscale.tol = 1e-08, 
                            mscale.rho.fun = 'bisquare')
                            

a <- pyinit(X=x, y=y, intercept=FALSE, deltaesc=dee, 
            cc.scale=control$tuning.chi, 
            prosac=control$prosac, clean.method=control$clean.method, 
            C.res = control$C.res, prop=control$prop, 
            py.nit = control$py.nit, en.tol=control$en.tol, 
            mscale.maxit = control$mscale.maxit, mscale.tol = control$mscale.tol,
            mscale.rho.fun=control$mscale.rho.fun)


