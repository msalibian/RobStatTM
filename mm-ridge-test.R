
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
# > a$coef
# [1]  0.56393725  4.93835073  4.52520916  4.81057373  4.02522828  3.38954222 -0.24710486 -0.76467020
# [9]  0.18671219  0.82159307  0.87416926  1.01648880  1.45493286 -0.05151432  1.48146142 -0.66388817
# [17] -0.60210763  0.27785843  0.30947435  0.26140655 -0.84834424 -0.82912644  0.71738652 -0.37778783
# [25]  0.17111429 -0.25918776 -0.73510224  0.33117573 -0.65798889 -0.29323547 -0.01326912  0.96685920
# [33] -0.45921634 -0.49199308  1.14162844  0.29633986  0.97414776 -0.67190607  1.74029930 -0.71029622
# [41]  0.80130134  0.33768376 -0.29952025  0.58196743 -0.29347014 -0.11529802  2.05014463 -1.43114321
# [49]  0.91086385  1.00114588 -0.04155250
a <- list(delta=0.2584052, lamda=11.42803)
b0 <- pense(X=x, y=y, alpha=0, standardize=TRUE, lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta))
# warnings
# 12: In pen_s_reg(std_data$xs, std_data$yc, alpha = alpha,  ... :
#                    PENSE did not converge for lambda = 8.838164

# bad results
as.vector(b0$coeff)
# [1] -0.3129575054  0.7120673657  0.5462823997  0.6968572031  0.2965445688  0.4881687731
# [7] -0.0419745712 -0.0188042829 -0.0244133721  0.1717148662  0.2422456816  0.1121325511
# [13]  0.2879602240 -0.1262394940  0.0678288025 -0.0325247768 -0.2796469626  0.0788652877
# [19] -0.1782404369  0.2057230278  0.0570493001  0.1354804073  0.0135514081  0.0567303180
# [25]  0.0128586139  0.0009493423 -0.2189636073  0.0701574901 -0.0091319053 -0.0733079254
# [31]  0.0651547632  0.1729963321  0.0890290482 -0.0183685251  0.1539796881  0.0056640699
# [37]  0.0354984634  0.0049737764  0.3181100881 -0.0539147355 -0.1239917759 -0.0978902874
# [43]  0.0665372294 -0.1422687666  0.0257184038  0.0262372269  0.1479607236 -0.1557955018
# [49] -0.1678854116  0.0895567412  0.0521972164

# the next ones do not work
b1 <- pense(X=x, y=y, alpha=0, standardize=TRUE, # lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta))
# Error in seq.default(log(lmin), log(lmax), length.out = nlambda) : 
#   'from' must be a finite number
b1 <- pense(X=x, y=y, alpha=0)
# Error in seq.default(log(lmin), log(lmax), length.out = nlambda) : 
#   'from' must be a finite number


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

