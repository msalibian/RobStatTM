
# Auxiliary functions 
# Tukey's rho
rho <- function(u, cc=1.5477) {
  w <- as.numeric( abs(u) <= cc )
  v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
  v <- v*6/cc^2
  return(v)
}

# Integral of rho with respect to a N(0,1)
rhoint <- function(e)
  return(integrate(function(a, cc) rho(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value)

# find constant for a desired BDP (delta)
find.tuning.chi <- function(delta, low=.5, upp=10) {
  return( uniroot( function(e) (rhoint(e)-delta), lower=low, upper=upp)$root )
}

# devtools::install_github('esmucler/mmlasso')
library(mmlasso)

# data
n <- 50 # 500
p <- 20 # 10
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(2, 5), rep(0, p-5))) + rnorm(n, sd=.5)

# compute S-ridge
a <- sridge(x=x, y=y, cualcv.S=5, numlam.S=30, niter.S=50, normin=0,
            denormout=0, alone=1, ncores=4)

# sanity check
# find tuning constant for adjusted delta
cc <- find.tuning.chi(a$delta)
# get residuals
re2 <- as.vector(y - x %*% a$coef[-1] - a$coef[1])
# check the M-scale equation
mean(rho(re2/a$scale, cc=cc))
a$delta
sum(rho(re2/a$scale, cc=cc))/(n - a$edf)

