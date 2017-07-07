
# Generate data
# X1 is very important
library(RobustStatistics)
n <- 50
set.seed(123)
x1 <- rnorm(n)
x2 <- rnorm(n, sd=1.5)
y <- rnorm(n, sd=.5) + 3*x1 + 5*x2

# adjusted r^2 seems too low
tmp2 <- lmrobdet(y~x1 + x2)
tmp2$r.squared
INVTR2(tmp2$r.squared, tmp2$control$tuning.psi)
summary(lm(y~x1+x2))$r.squared



# Compute it by hand:

# objective function at the MM-estimator
(s2 <- sum(rho(resid(tmp2)/tmp2$scale, cc=tmp2$control$tuning.psi)))

# minimum of G(alpha) = sum(rho( (y-alpha) / sigma ) )
#
# where sigma is the scale obtained with the full model
#
# we use IRWLS starting from median(y)

a <- as.vector(refine.sm(x=matrix(rep(1,n), n, 1), y=y, initial.beta=median(y),
  initial.scale=tmp2$scale, k=200, conv=1, cc=tmp2$control$tuning.psi, step='M')$beta)

# objective function at this minimum
(s02 <- sum(rho((y-a)/tmp2$scale, cc=tmp2$control$tuning.psi)))

# check that the minimum is better than at the median, say:
(sum(rho((y-median(y))/tmp2$scale, cc=tmp2$control$tuning.psi)))

# check that it is indeed a minimum

oo <- tt <- seq(-1, 1, length=500)
for(j in 1:500) oo[j] <- sum(rho((y-tt[j])/tmp2$scale, cc=tmp2$control$tuning.psi))
plot(tt, oo, type='l', ylim=c(38, 40))
abline(v=a)
abline(h=s02)

# now the R^2
(s02 - s2)/s02

# With LS is very high, as you'd expect
summary(lm(y~x1))$r.squared

# Using lmrob in robustbase also gives a much higher value
tmp2 <- lmrob(y~x1, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
summary(tmp2)$r.squared



