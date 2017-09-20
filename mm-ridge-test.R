
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
b0 <- pense(X=x, y=y, alpha=1e-9, standardize=TRUE, lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta))
# bad results
# > as.vector(b0$coeff)
# [1] -0.3129452967  0.7121311101  0.5461833697  0.6968806820  0.2965058615  0.4881604731 -0.0419973497
# [8] -0.0187386321 -0.0243840548  0.1716804345  0.2422068259  0.1122685594  0.2880183202 -0.1262386725
# [15]  0.0677409486 -0.0326192023 -0.2795226942  0.0789656278 -0.1781787701  0.2058248845  0.0570528406
# [22]  0.1354692827  0.0136279562  0.0567609559  0.0129764541  0.0008582584 -0.2190313445  0.0702381527
# [29] -0.0092658433 -0.0732404819  0.0650584955  0.1730236542  0.0890287023 -0.0185502604  0.1539539257
# [36]  0.0057937459  0.0353332836  0.0049614192  0.3182054685 -0.0540393799 -0.1239302509 -0.0979559030
# [43]  0.0664725209 -0.1421960800  0.0257849657  0.0261492142  0.1480458235 -0.1559357246 -0.1680001431
# [50]  0.0896082852  0.0522610487
b0$scale
# scale
# 14.01923

b1 <- pense(X=x, y=y, alpha=1e-9, standardize=TRUE, # lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta))
# bad results?
as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
# [1] -2.510333e-01  4.929228e-04  3.895520e-04  4.941845e-04  1.967635e-04  3.536791e-04 -4.005779e-05
# [8] -4.223980e-06 -2.525764e-05  1.292631e-04  1.832549e-04  7.477593e-05  2.164555e-04 -9.971194e-05
# [15]  5.124504e-05 -5.541110e-06 -2.076651e-04  5.012369e-05 -1.432063e-04  1.583820e-04  6.275588e-05
# [22]  9.596635e-05  9.139492e-06  3.571185e-05  9.430201e-06  3.666528e-06 -1.571654e-04  5.546936e-05
# [29] -4.307249e-06 -7.372866e-05  4.755529e-05  1.158672e-04  7.330656e-05 -1.019428e-05  1.175726e-04
# [36]  8.749737e-06  8.490345e-06  7.032498e-06  2.331767e-04 -4.832811e-05 -7.352491e-05 -6.609487e-05
# [43]  5.763831e-05 -1.192467e-04  1.964519e-05  1.299670e-05  1.131480e-04 -1.085018e-04 -1.373534e-04
# [50]  6.775017e-05  4.249457e-05
b1$scale
# [1] 16.02063 16.02098 16.02125 16.02146 16.02161 16.02173 16.02182 16.02188 16.02193 16.02197 16.02200 16.02202
# [13] 16.02204 16.02205 16.02206 16.02206 16.02207 16.02207 16.02208 16.02208 16.02208 16.02208 16.02208 16.02208
# [25] 16.02208 16.02208 16.02208 16.02208 16.02208 16.02208 16.02209 16.02209 16.02209 16.02209 16.02209 16.02209
# [37] 16.02209 16.02209 16.02209 16.02209 16.02209 16.02209 16.02209 16.02209 16.02209 16.02209 16.02209 16.02209
# [49] 16.02209 16.02209

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

