
## This example bombs in my linux machine
##

library(pense)
library(mmlasso)
n <- 80
p <- 50
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 5), rep(0, p-5))) + rnorm(n, sd=.5)
a <- sridge(x=x, y=y, cualcv.S=5, numlam.S=30, niter.S=50, normin=0,
            denormout=0, alone=1, ncores=4)
b0 <- pense(X=x, y=y, alpha=0, standardize=TRUE, lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta))
# bad results, warnings
as.vector(b0$coeff)
# the next ones do not work
b1 <- pense(X=x, y=y, alpha=0, standardize=TRUE, # lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta))
b1 <- pense(X=x, y=y, alpha=0)


# Check MM-ridge regression


rho <- function(u, cc=1.5477) {
  w <- as.numeric( abs(u) <= cc )
  v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
  v <- v*6/cc^2
  return(v)
}


rhoint <- function(e)
  return(integrate(function(a, cc) rho(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value)


find.tuning.chi <- function(delta, low=.5, upp=10) {
  return( uniroot( function(e) (rhoint(e)-delta), lower=low, upper=upp)$root )
}

# devtools::install_github('gcohenfr/PENSE-package', ref='develop', auth_token='d49f60185d9de80d37dacef7f93b14df907e6af7')
library(pense)

# compare S-ridge as computed with pense
# and that computed with ezequiel's package

# devtools::install_github('esmucler/mmlasso')
library(mmlasso)

# library(glmnet)

n <- 50 # 500
p <- 20 # 10
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(2, 5), rep(0, p-5))) + rnorm(n, sd=.5)
# mmlasso
a <- sridge(x=x, y=y, cualcv.S=5, numlam.S=30, niter.S=50, normin=0,
            denormout=0, alone=1, ncores=4)
# Optimal lambda
a$lamda
# right-hand side of the M-scale used?
a$delta

# pense ridge?
b0 <- pense(X=x, y=y, alpha=0, standardize=TRUE, lambda=1e-9, initial='cold',
           options=pense_options(delta=a$delta))

cbind(a$coef, as.vector(b0$coef[,1]))
c(a$scale, b0$scale)

## MM-ridge?

g <- mstep(b0, complete_grid=TRUE)
cbind(a$coef, as.vector(b0$coef[,1]), g$coefficients[,1])

c(a$scale, b0$scale, g$scale)

# d <- pense::elnet(X=x, y=y, alpha=0, lambda=1e-9, addLeading1s=TRUE)

# as.vector(d$coef[,1])

re <- as.vector(y - x %*% (b0$coef[,1])[-1] - b0$coef[1,1])
mean(rho(re/b0$scale))

cc <- find.tuning.chi(a$delta)
mean(rho(rnorm(1e6), cc=cc))
re2 <- as.vector(y - x %*% a$coef[-1] - a$coef[1])
sum(rho(re2/a$scale, cc=cc))/(n - a$edf)
mean(rho(re2/a$scale, cc=cc))
a$delta

# b1 <- pense(X=x, y=y, alpha=0, lambda=seq(0, 5, length=20), standardize=TRUE,
#            control=pense.control(mscale.delta = 0.49))
#
#
# d <- cv.glmnet(x=x, y=y, lambda=seq(0, 5, length=20), family='gaussian',
#                intercept=TRUE, alpha=0) # , nlambda=50)

