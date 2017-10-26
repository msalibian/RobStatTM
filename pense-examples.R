
library(devtools)
install_github("dakep/pense-rpkg", ref = "develop", force=TRUE)




library(pense)
n <- 80
p <- 50
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 5), rep(0, p-5))) + rnorm(n, sd=.5)
# Compute the MM-ridge (EN with alpha = 0)
set.seed(123)
b0 <- pensem(x=x, y=y, alpha=0, ncores=4, nlambda=30)
coef(b0)

set.seed(123)
h0 <- pense(x=x, y=y, alpha=0, ncores=4, nlambda=30)
h1 <- pensem(h0, alpha=0, ncores=4, nlambda=30)

cbind(coef(b0), coef(h1))

coef(h0)


library(pense)
n <- 80
p <- 50
set.seed(123)
x <- matrix(rnorm(n*p), n, p)
y <- as.vector( x %*% c(rep(7, 5), rep(0, p-5))) + rnorm(n, sd=.5)

# Compute the MM-ridge (EN with alpha = 0)
b0 <- pense(X=x, y=y, alpha=0, ncores=4, nlambda=30)
b1 <- pensem(b0, alpha=0, ncores=4, nlambda=30)


