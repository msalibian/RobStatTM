
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
b0 <- pense(X=x, y=y, alpha=1e-9, standardize=TRUE, lambda=a$lamda/nrow(x), initial='cold',
            options=pense_options(delta=a$delta),
            en_options = en_options_dal(),
            init_options = initest_options(psc_method = "exact"))
summary(b0$lambda)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1429  0.1429  0.1429  0.1429  0.1429  0.1429
as.vector(b0$coeff)
# [1]  0.099931035  5.127215364  4.103291233  3.787796392  2.825800485  3.168253233  0.104835737  0.668650140
# [9]  0.129847273  0.797997409 -0.087012260  0.371518988  0.240088638  0.157939412  0.560327771  0.004954685
# [17] -1.138523415 -0.180187505 -0.645891444  0.745576307  0.808967413  0.325132700  0.690614906  0.454194649
# [25] -0.291426977  0.454984838 -2.403653159 -0.927144384 -0.867659328 -0.239478900  0.557485650  1.363572729
# [33]  0.017722912 -1.097852964  1.160609899  0.073932434  0.723816860  0.246106887  2.683474710  0.043096216
# [41] -0.998344691 -0.617306568 -0.321294108 -0.553822742 -0.544468600  0.618603074  0.038421376 -1.735686118
# [49]  0.084093221  0.922850352  0.631343409
b0$scale
# scale
# 3.109854

b1 <- pense(X=x, y=y, alpha=1e-9, standardize=TRUE, # lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta),
            en_options = en_options_dal(),
            init_options = initest_options(psc_method = "exact"),
            lambda_min_ratio = 1e-14)
summary(b1$lambda)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.000e+00 1.000e+00 2.019e+03 7.944e+08 6.277e+06 1.915e+10
b1$lambda_opt
# [1] 0.0007136994
as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
# [1] -0.080953012  7.155850862  7.235214622  6.946541552  6.825808803  6.849400621  0.076929068 -0.004258557
# [9] -0.030491592 -0.120645557  0.245846521  0.001669768 -0.200533405 -0.061332792 -0.053580889 -0.047906441
# [17]  0.190418658  0.154841152 -0.225974463  0.228098995  0.017502687 -0.184653276 -0.004838187 -0.053616424
# [25] -0.093718216  0.079654356  0.107040429 -0.167083521  0.023726712  0.023808286  0.049977395 -0.198843205
# [33] -0.158741526  0.131401514 -0.112381100 -0.305308955  0.130415008  0.004798759 -0.028150028 -0.044722125
# [41]  0.093824018 -0.069953968 -0.092755322 -0.145048917  0.006396454 -0.067620122  0.080991754 -0.168580189
# [49]  0.080324727  0.042955790 -0.025792814
b1$scale[b1$lambda == b1$lambda_opt]
# [1] 0.06009107





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

