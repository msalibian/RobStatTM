
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

# try to reproduce srdige results:
b0 <- pense(X=x, y=y, alpha=0, standardize=TRUE, lambda=a$lamda/nrow(x), initial='cold',
            options=pense_options(delta=a$delta),
            en_options = en_options_dal(),
            init_options = initest_options(psc_method = "exact"))
summary(b0$lambda)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1429  0.1429  0.1429  0.1429  0.1429  0.1429
as.vector(b0$coeff)
# [1] -0.56188256  5.35716464  4.45576453  3.83223186  2.57619773  3.27338500
# [7]  0.59675519  0.62648567  0.30245205  0.52771038  0.20992737  0.20518861
# [13]  0.18054032  0.10602851  0.33677310  0.35973216 -1.07022291 -0.65389808
# [19] -0.15588146  1.11045653  0.76979123  0.11396488  0.47741425  0.69250752
# [25] -0.14571474  0.74995238 -2.53510696 -1.29878442 -0.63056945 -0.41223131
# [31]  0.50595632  1.98362399  0.01130509 -0.66930877  0.77941751 -0.22201319
# [37]  0.39108893  0.34179955  2.48496102 -0.51311224 -0.82466886 -0.72035204
# [43] -0.53518098 -0.27894737 -0.32652570  0.64891310  0.02405424 -1.43743238
# [49] -0.31805878  0.88437940  1.01986135
b0$scale
# scale
# 2.869959


# a better S-ridge?
b1 <- pense(X=x, y=y, alpha=0, standardize=TRUE, # lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta),
            en_options = en_options_dal(),
            init_options = initest_options(psc_method = "exact")) #,
            #lambda_min_ratio = 1e-14)
summary(b1$lambda)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0019    0.0610    1.9337  155.8592   60.9775 1914.6384 
b1$lambda_opt
# [1] 0.01377938
as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
# [1] -0.257865402  6.212883217  5.593609847  5.875579238  5.234586110  5.124939559 -0.275942458
# [8] -0.847478717  0.032437231  0.852340246  0.395747959 -0.040367337 -0.030968139  0.438158608
# [15]  0.010852623  0.838641428 -0.942082755  0.029412031  0.002874029  1.143015109 -0.873840206
# [22] -0.339214340  0.511856609 -0.330053776 -0.240110088  1.021959134 -1.104889252 -0.505542362
# [29] -1.540149232 -1.252338673  0.465236057 -0.099322511  0.178829052  0.128735506  0.430600112
# [36]  0.230572202  0.813284759  1.094679058 -0.051851676 -0.556864101 -0.728261017 -1.309801010
# [43]  1.041464717 -0.295965949 -0.430933697  0.620052550  0.212928294 -1.560282646 -0.158097917
# [50]  0.357026514  0.331336697
b1$scale[b1$lambda == b1$lambda_opt]
# [1] 0.6874806

system.time(
  b1 <- pense(X=x, y=y, alpha=0, standardize=TRUE, options=pense_options(delta=a$delta))
)
summary(b1$lambda)
b1$lambda_opt
as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
b1$scale[b1$lambda == b1$lambda_opt]






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

