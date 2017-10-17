

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

cbind(as.vector(a$coef.MMLasso), coef(b2))

# Again with different values of n & p

n <- 50
p <- 100
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 25), rep(0, p-25))) + rnorm(n, sd=.5)

a <- mmlasso(x=x, y=y, ncores=2)
b1 <- pense(X=x, y=y, alpha=0, nlambda=30, ncores=2) #init_options = initest_options(maxit=1, maxit_pense_refinement=1))
b2 <- pensem(b1, alpha=1, nlambda=30, ncores=2)
cbind(as.vector(a$coef.MMLasso), coef(b2))




