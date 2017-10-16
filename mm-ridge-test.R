
# Install pense
# devtools::install_github("dakep/pense-rpkg", ref = "develop", force=TRUE)

# Install pyinit
# devtools::install_github("dakep/pyinit", ref = "master") # , force=TRUE)


## This example bombs in my linux machine
##

# detach(package:glmnet)
library(pense)
library(mmlasso)





n <- 80
p <- 50
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 5), rep(0, p-5))) + rnorm(n, sd=.5)
system.time(
  a.s <- sridge(x=x, y=y, cualcv.S=5, numlam.S=30, niter.S=50, normin=0,
            denormout=0, alone=1, ncores=4)
)

# try to reproduce srdige results, pick optimal lambda

system.time(
  b0 <- pense(X=x, y=y, alpha=0, standardize=TRUE, initial='cold', # lambda=a$lamda/nrow(x), 
              options=pense_options(delta=a.s$delta), ncores=4, 
              init_options = initest_options( # psc_method = "exact", 
                maxit = 1, # nr. of refinement steps in ENPY
                maxit_pense_refinement = 1 # nr. of PENSE refinements on all initial estimates before we pick the best 5 candidates
              ))
)

# more expensive, better S-ridge
system.time(
  b0.1 <- pense(X=x, y=y, alpha=0, standardize=TRUE, # initial='cold', # lambda=a$lamda/nrow(x), 
              options=pense_options(delta=a.s$delta), ncores=4,
              # init_options = initest_options( # psc_method = "exact", 
              #   maxit = 1, # nr. of refinement steps in ENPY
              #   maxit_pense_refinement = 1 # nr. of PENSE refinements on all initial estimates before we pick the best 5 candidates
              # )
              )
)


cbind(as.vector(a.s$coef), coef(b0), coef(b0.1))


## MM-LASSO

a <- mmlasso(x=x, y=y, ncores=4)
# b1 <- pense(X=x, y=y, alpha=0, nlambda=30, ncores=2) 
b2 <- pensem(b0.1, alpha=1, nlambda=30, ncores=4)

cbind(as.vector(a$coef.MMLasso), coef(b2))




# Another example
n <- 200
p <- 100
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 25), rep(0, p-25))) + rnorm(n, sd=1.5)

a <- mmlasso(x=x, y=y, ncores=4)
b1 <- pense(X=x, y=y, alpha=0, nlambda=30, ncores=4)
b2 <- pensem(b1, alpha=1, nlambda=30, ncores=4)
cbind(as.vector(a$coef.MMLasso), coef(b2))








# a better S-ridge?
system.time(
b1 <- pense(X=x, y=y, alpha=0, standardize=TRUE, # lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta),
            en_options = en_options_dal(),
            init_options = initest_options(psc_method = "exact")) #,
            #lambda_min_ratio = 1e-14)
)

system.time(
b1 <- pense(X=x, y=y, alpha=0, standardize=TRUE, # lambda=a$lamda, initial='cold',
            options=pense_options(delta=a$delta),
            init_options = initest_options(maxit=1, maxit_pense_refinement=1)) #,
#lambda_min_ratio = 1e-14)
)
a$coef
as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])


summary(b1$lambda)
b1$lambda_opt
as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
b1$scale[b1$lambda == b1$lambda_opt]

# Win
# user  system elapsed 
# 132.93    0.10  134.05 
# There were 50 or more warnings (use warnings() to see the first 50)
# > summary(b1$lambda)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0019    0.0610    1.9337  155.8592   60.9775 1914.6384 
# > b1$lambda_opt
# [1] 0.01826748
# > as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
# [1] -0.268118982  6.124693335  5.479921016  5.802224408  5.084074750  4.969363031 -0.274478511 -0.863227108
# [9]  0.029051857  0.890374412  0.418708030 -0.013776165 -0.028066569  0.484207896  0.008920756  0.886381106
# [17] -1.002247976  0.031423460 -0.009417850  1.204134302 -0.907760013 -0.337956634  0.548441991 -0.366367489
# [25] -0.238992297  1.091238943 -1.191603686 -0.532714817 -1.624271053 -1.314018251  0.489604173 -0.066118251
# [33]  0.194560874  0.120391313  0.480647068  0.278434277  0.854657690  1.169466959 -0.035295640 -0.590824710
# [41] -0.783109110 -1.388988857  1.095519014 -0.304125507 -0.454783769  0.649150350  0.201214173 -1.655888217
# [49] -0.214979744  0.390543444  0.365550743
# > b1$scale[b1$lambda == b1$lambda_opt]
# [1] 0.7728777

# Linux
# user  system elapsed
# user  system elapsed
# 180.776 795.112 245.299
# There were 50 or more warnings (use warnings() to see the first 50)
# > summary(b1$lambda)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0019    0.0610    1.9340  155.9000   60.9800 1915.0000
# > b1$lambda_opt
# [1] 0.01377938
# > as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
# [1] -0.257865438  6.212883217  5.593609814  5.875579211  5.234586081
# [6]  5.124939520 -0.275942455 -0.847478702  0.032437240  0.852340264
# [11]  0.395747977 -0.040367321 -0.030968148  0.438158618  0.010852610
# [16]  0.838641431 -0.942082754  0.029412018  0.002874041  1.143015171
# [21] -0.873840211 -0.339214341  0.511856621 -0.330053786 -0.240110077
# [26]  1.021959158 -1.104889266 -0.505542349 -1.540149259 -1.252338675
# [31]  0.465236032 -0.099322545  0.178829068  0.128735498  0.430600107
# [36]  0.230572195  0.813284790  1.094679086 -0.051851691 -0.556864143
# [41] -0.728261041 -1.309801015  1.041464717 -0.295965958 -0.430933725
# [46]  0.620052543  0.212928293 -1.560282666 -0.158097928  0.357026528
# [51]  0.331336714
# > b1$scale[b1$lambda == b1$lambda_opt]
# [1] 0.6874806



system.time(
  b1 <- pense(X=x, y=y, alpha=0, standardize=TRUE, options=pense_options(delta=a$delta))
)
summary(b1$lambda)
b1$lambda_opt
as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
b1$scale[b1$lambda == b1$lambda_opt]

# Windows
# user  system elapsed 
# 563.56    0.24  576.45 
# > summary(b1$lambda)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0019    0.0610    1.9337  155.8592   60.9775 1914.6384 
# > b1$lambda_opt
# [1] 0.005914019
# > as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
# [1]  0.21831799  6.95768346  5.71701138  5.96815432  4.26217773  5.15288984  1.11894491 -0.02227110  0.35026196
# [10] -0.06387154  0.58899811  0.52206940  0.37268119 -0.06401331  0.87112159 -0.21254100 -2.46299331  0.84307669
# [19] -0.26098040  1.07134889 -0.50477414 -0.08289497  0.38244137 -0.98500711 -0.26607553  0.64197412 -2.27387589
# [28] -0.74507596 -1.07266345  0.99181954  0.68863472  0.08978709 -0.99144962 -0.93299054  0.08375649 -0.71554213
# [37]  0.41921825  0.03594240  0.74546749 -0.56040232 -0.54954768 -0.58049474  0.43687705  1.42860082 -0.77560776
# [46]  0.08406832 -0.02274564 -1.14561040 -0.27014853 -0.35032377  0.32538113
# > b1$scale[b1$lambda == b1$lambda_opt]
# [1] 0.444755
# > sessionInfo()
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252
# [4] LC_NUMERIC=C                    LC_TIME=English_Canada.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] mmlasso_1.3.4 pense_1.0.0   Matrix_1.2-10
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.12      lattice_0.20-35   codetools_0.2-15  foreach_1.4.3     MASS_7.3-47       plyr_1.8.4       
# [7] grid_3.4.1        gtable_0.2.0      scales_0.5.0      ggplot2_2.2.1     rlang_0.1.2       lazyeval_0.2.0   
# [13] robustHD_0.5.1    perry_0.2.0       doParallel_1.0.10 robustbase_0.92-7 iterators_1.0.8   tools_3.4.1      
# [19] munsell_0.4.3     DEoptimR_1.0-8    parallel_3.4.1    compiler_3.4.1    colorspace_1.3-2  tibble_1.3.4     
 
# Linux
# user  system elapsed
# 858.464 670.060 753.594
# There were 50 or more warnings (use warnings() to see the first 50)
# > summary(b1$lambda)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0019    0.0610    1.9340  155.9000   60.9800 1915.0000
# > b1$lambda_opt
# [1] 0.01826748
# > as.vector(b1$coeff[, b1$lambda == b1$lambda_opt])
# [1] -0.268150156  6.124670867  5.479865045  5.802177838  5.083979300
# [6]  4.969298030 -0.274453261 -0.863237644  0.029037443  0.890409546
# [11]  0.418715672 -0.013737519 -0.028063133  0.484241176  0.008875158
# [16]  0.886437160 -1.002319599  0.031407947 -0.009390470  1.204267971
# [21] -0.907841266 -0.337972615  0.548503829 -0.366403587 -0.238968505
# [26]  1.091303477 -1.191666383 -0.532717236 -1.624375871 -1.314034171
# [31]  0.489622257 -0.066221722  0.194580523  0.120395834  0.480654388
# [36]  0.278398531  0.854704571  1.169484691 -0.035310982 -0.590871622
# [41] -0.783136429 -1.389035890  1.095573700 -0.304112473 -0.454849920
# [46]  0.649152126  0.201210402 -1.655877217 -0.215001902  0.390544140
# [51]  0.365551967
# > b1$scale[b1$lambda == b1$lambda_opt]
# [1] 0.7728996
# > sessionInfo()
# R version 3.3.1 (2016-06-21)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.3 LTS
# 
# locale:
#   [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8
# [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8
# [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#   [1] mmlasso_1.3.4 pense_1.0.0   Matrix_1.2-6
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.11      MASS_7.3-45       munsell_0.4.3     doParallel_1.0.10
# [5] colorspace_1.2-6  lattice_0.20-33   foreach_1.4.3     plyr_1.8.3
# [9] tools_3.3.1       parallel_3.3.1    grid_3.3.1        gtable_0.2.0
# [13] perry_0.2.0       iterators_1.0.8   robustHD_0.5.1    lazyeval_0.2.0
# [17] assertthat_0.1    tibble_1.2        ggplot2_2.2.1     codetools_0.2-14
# [21] robustbase_0.92-7 compiler_3.3.1    DEoptimR_1.0-8    scales_0.4.1




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

