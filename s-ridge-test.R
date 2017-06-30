
# Auxiliary functions 
# Tukey's rho
rho <- function(u, cc=1.5477) {
  w <- as.numeric( abs(u) <= cc )
  v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
  v <- v*6/cc^2
  return(v)
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
# get residuals
re2 <- as.vector(y - x %*% a$coef[-1] - a$coef[1])
# check the M-scale equation
mean(rho(re2/a$scale))
a$delta
sum(rho(re2/a$scale, cc=cc))/(n - a$edf)

