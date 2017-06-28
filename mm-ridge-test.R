
# Check MM-ridge regression

library(pense)

# compare S-ridge as computed with pense
# and that computed with ezequiel's package

# devtools::install_github('esmucler/mmlasso')
library(mmlasso)

n <- 50
p <- 20
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(2, 5), rep(0, p-5))) + rnorm(n, sd=.5)
a <- sridge(x=x, y=y, cualcv.S=5, numlam.S=30, niter.S=50, normin=0, 
            denormout=0, alone=1, ncores=1)

b <- pense(X=x, y=y, alpha=0, lambda=0, nlambda=1, standardize=TRUE, 
           control=pense.control(mscale.delta = 0.3))

d <- pense(X=x, y=y, alpha=0, nlambda=50, standardize=TRUE, 
           control=pense.control(mscale.delta = 0.3))

