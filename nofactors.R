
library(RobustStatistics)
set.seed(123)
n <- 200
x1 <- rnorm(n)
x2 <- runif(n)
x3 <- rexp(n, rate=.2)
x4 <- rnorm(n)
x5 <- rnorm(n, sd=2)
x6 <- rbinom(n, size=3, prob=.3)
x7 <- rt(n, df=2)
x8 <- rpois(n, lambda=5)
f1 <- as.factor( LETTERS[ rbinom(n, size=2, prob=.3) + 1] )
f2 <- factor( c('Low', 'Med', 'High')[ rbinom(n, size=2, prob=.3) + 1], levels=c('Low', 'Med', 'High'))
y <- rnorm(n, sd=.5) + 2*x1  # - 3*as.numeric(f2)
d <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, x6=x6, x7=x7, x8=x8, f1=f1, f2=f2, y = y)

mf <- model.frame(y ~ x1 + x2 + f1 + f2, data=d)
#mf <- model.frame(y ~ f1 + f2, data=d)
split <- splitFrame(mf=mf, type='f')
Z <- split$x1
X <- split$x2
n <- nrow(X)
q <- ncol(Z)
p <- ncol(X)
if(p==0) print('No continuous variables')


