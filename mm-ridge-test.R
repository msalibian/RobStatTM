
## This example bombs in my linux machine
##

detach(package:glmnet)
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
b0 <- pense(X=x, y=y, alpha=1e-3, standardize=TRUE, lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta))
# bad results
as.vector(b0$coeff)
# [1] -3.132728e-01  7.119411e-01  5.459432e-01  6.965050e-01  2.959825e-01  4.878504e-01 -4.091106e-02
# [8] -1.776365e-02 -2.355424e-02  1.709769e-01  2.414012e-01  1.115445e-01  2.873882e-01 -1.255835e-01
# [15]  6.687172e-02 -3.154495e-02 -2.787938e-01  7.805673e-02 -1.774030e-01  2.051673e-01  5.620912e-02
# [22]  1.349359e-01  1.277071e-02  5.575485e-02  1.228093e-02  9.559815e-05 -2.183080e-01  6.925380e-02
# [29] -8.225870e-03 -7.242959e-02  6.428259e-02  1.720992e-01  8.834833e-02 -1.762954e-02  1.530641e-01
# [36]  4.732597e-03  3.425310e-02  4.057803e-03  3.176702e-01 -5.325698e-02 -1.228324e-01 -9.702156e-02
# [43]  6.566646e-02 -1.413014e-01  2.502185e-02  2.539448e-02  1.472589e-01 -1.549121e-01 -1.672725e-01
# [50]  8.878876e-02  5.150094e-02
b0$scale
# scale
# 14.02308

b1 <- pense(X=x, y=y, alpha=1e-9, standardize=TRUE, # lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta),
            en_options = en_options_dal(),
            init_options = initest_options(psc_method = "exact"),
            lambda_min_ratio = 1e-14)
b1$lambda_opt
# [1] 0.05914019
as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
# [1]  0.33962275  5.99644725  5.03823998  4.40996312  2.43687458  3.78214838  0.01263743  0.18398226  0.51121438
# [10]  0.25186233  0.59657544  0.96749781  0.44741765  0.83437221 -0.20781775  0.81470387 -1.82038559  0.58098590
# [19]  0.25428982  0.31830200  2.18900085  0.16076587  0.69837584 -1.43335399 -0.53735748  0.59151405 -1.35005488
# [28] -0.78569416 -0.58204451 -1.53155351  0.22011564  2.33020436 -0.61899168 -0.09000444 -0.44586815  0.14299942
# [37]  0.93085374  1.03178662  0.65323495 -0.95956967 -0.78936411 -0.98212292  0.61499226  1.03615548  0.24066852
# [46]  0.12322585  0.32033421 -0.85020645 -0.26709221  0.62086677  1.36270552
b1$scale
# [1]  0.9173721  1.0554860  1.1910704  1.4652329  1.6827527  1.9659139  2.2887926  2.6909651  3.1760354  3.7562661
# [11]  4.3138571  5.0493622  5.8975691  6.6458776  7.4845584  8.3709777  9.2396228 10.0941083 10.9384758 11.7385632
# [21] 12.4721064 13.1257456 13.6925434 14.1719315 14.5685482 14.8905783 15.1479875 15.3511366 15.5098350 15.6328131
# [31] 15.7275099 15.8000735 15.8554654 15.8976222 15.9296339 15.9539023 15.9722738 15.9861587 15.9966362 16.0044659
# [41] 16.0102592 16.0145134 16.0174642 16.0194617 16.0207600 16.0216353 16.0220831 16.0220855 16.0220855 16.0220855

## how do the "classical" ones work?
##

detach(package:pense)
library(glmnet)
d <- glmnet(x=x, y=y, alpha=0, family='gaussian', standardize=TRUE, intercept=TRUE)
d2 <- cv.glmnet(x=x, y=y, alpha=0, nfolds=5, family='gaussian', standardize=TRUE, intercept=TRUE)
as.vector(coef(d2, s='lambda.min') )
# [1] -0.44758993  6.84071517  6.43792473  6.11671099  5.80955408  5.63654031  0.27685178  0.18691536  0.22341984
# [10]  0.21037653  0.31719134  0.12330850  0.35881440 -0.18028603  0.26846078 -0.06692925 -0.21768027  0.10069179
# [19] -0.30636298  0.38365965 -0.16083863  0.07292750 -0.06625486 -0.37753467 -0.25670810  0.43897691 -0.76350800
# [28] -0.43616188 -0.20168467 -0.32243167 -0.14685351  0.28685218 -0.11230243  0.24661119 -0.06925378  0.12252345
# [37]  0.50052367  0.24918111  0.43636937 -0.31912413 -0.32392971 -0.34499128  0.16741021  0.16953670 -0.24950992
# [46]  0.04629432  0.29358932 -0.59842793 -0.44200748  0.15908943  0.02100474


# # Error in seq.default(log(lmin), log(lmax), length.out = nlambda) :
# #   'from' must be a finite number
# b1 <- pense(X=x, y=y, alpha=0)
# # Error in seq.default(log(lmin), log(lmax), length.out = nlambda) :
# #   'from' must be a finite number
#

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

