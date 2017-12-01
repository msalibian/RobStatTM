

# ## Only run once
# # Install pense package
library(devtools)
install_github("dakep/pense-rpkg", ref = "develop", force=TRUE)
#
# # Install pyinit package
# library(devtools)
# install_github("dakep/pyinit", ref = "master", force=TRUE)


library(pense)

n <- 80
p <- 50
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 5), rep(0, p-5))) + rnorm(n, sd=.5)


# Compute the S-ridge (EN with alpha = 0)
b0 <- pense(X=x, y=y, alpha=0, ncores=4, nlambda=30)

# Compute MM-LASSO starting from the S-ridge (and the assoc. scale)
b2 <- pensem(b0, alpha=1, nlambda=30, ncores=4)

# Compare with Ezequiel's MM-LASSO
library(mmlasso)
a <- mmlasso(x=x, y=y, ncores=4)

round(cbind(as.vector(a$coef.MMLasso), coef(b2)), 3)

# Again with different values of n & p

n <- 75
p <- 100
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 25), rep(0, p-25))) + rnorm(n, sd=.5)

a <- mmlasso(x=x, y=y, ncores=2)
b1 <- pense(X=x, y=y, alpha=0, nlambda=30, ncores=2) #init_options = initest_options(maxit=1, maxit_pense_refinement=1))
b2 <- pensem(b1, alpha=1, nlambda=30, ncores=2)
round(cbind(as.vector(a$coef.MMLasso), coef(b2)), 3)

# install.packages('cellWise')
# data(glass, package='cellWise')
# data(glassVessels, package='locout')

load('glassVessels.RData')
x <- glass$data
# table(glass$group)
# matplot(t(x), col='gray', type='l', lty=1)
y <- glass$chemical[,'PbO']
x <- x[, 15:500]
library(glmnet)
set.seed(123)
tmp <- cv.glmnet(y=y, x=x, nfolds=10, alpha=1)
plot(tmp)
as.vector( coef(tmp, s='lambda.1se') )

library(mmlasso)
set.seed(1)
a <- mmlasso(x=x, y=y, ncores=4)
a$coef.MMLasso[ a$coef.MMLasso != 0 ]
a$coef.MMLasso.ad[ a$coef.MMLasso.ad !=0 ]
which(a$coef.MMLasso!=0)
# [1]   1  13 130 322 323 357 359 388
which(a$coef.MMLasso.ad!=0)
# [1]   1  13 322 359

set.seed(123)
a2 <- mmlasso(x=x, y=y, ncores=4)
a2$coef.MMLasso[ a2$coef.MMLasso != 0 ]
a2$coef.MMLasso.ad[ a2$coef.MMLasso.ad !=0 ]
which(a2$coef.MMLasso!=0)
which(a2$coef.MMLasso.ad!=0)

library(pense)
# Compute the S-ridge (EN with alpha = 0)
b0 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
# Compute MM-LASSO starting from the S-ridge (and the assoc. scale)
b2 <- pensem(b0, alpha=1, nlambda=30, ncores=4)

