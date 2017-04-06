
set.seed(123)
n <- 50
beta <- 0
x <- rnorm(n)
u <- rnorm(n)
y <- x * beta + u
x <- c(x, rep(10, 3))
y <- c(y, rep(20, 3))

source('R/DCML.R')

# trimmed MSE function
# returns the average of the (1-alpha)100%
# smallest elements in "x", each of them squared
tm <- function(x, alpha) {
  n <- length(x)
  n0 <- floor(alpha * n)
  n <- n - n0
  return( mean( (sort(x^2))[1:n] ) )
}

bb <- seq(-2, 3, length=100)
nbb <- length(bb)
lms.ob <- lts.ob <- s.ob <- ls.ob <- rep(NA, nbb)
for(i in 1:nbb) {
  re <- y - bb[i]*x
  s.ob[i] <- mscale(re, tol=1e-7)
  ls.ob[i] <- mean(re^2)
  lms.ob[i] <- median(re^2)
  lts.ob[i] <- tm(re, alpha=.5)
}

s.ob <- s.ob/max(s.ob)
ls.ob <- ls.ob/max(ls.ob)
lts.ob <- lts.ob/max(lts.ob)
lms.ob <- lms.ob/max(lms.ob)

plot(bb, ls.ob, type='l', lwd=3, col='gray30')
lines(bb, s.ob, lwd=3, lty=1, col='gray50')
lines(bb, lts.ob, lwd=3, lty=1, col='gray70')
lines(bb, lms.ob, lwd=3, lty=1, col='black')

