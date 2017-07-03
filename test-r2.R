

n <- 50
set.seed(123)
x1 <- rnorm(n)
y <- rnorm(n, sd=.1) + x1

# refine.sm(x=matrix(rep(1,n), n, 1), y=y, initial.beta=median(y), initial.scale=mad(y), k=20,
#                       conv=1, cc=lmrobdet.control()$tuning.psi, step='M')

tmp2 <- lmrobdet(y~x1)
tmp0 <- lmrobdet(y~1)

# refine.sm(x=matrix(rep(1,n), n, 1), y=y, initial.beta=median(y), initial.scale=mad(y), k=20,
#           conv=1, cc=lmrobdet.control()$tuning.psi, step='M')

tmp2$r.squared

(s2 <- sum(rho(resid(tmp2)/tmp2$scale, cc=tmp2$control$tuning.psi)))
(s02 <- sum(rho(resid(tmp0)/tmp2$scale, cc=tmp2$control$tuning.psi)))
# (s02 <- sum(rho((y-median(y))/tmp2$scale, cc=tmp2$control$tuning.psi)))
(s02 - s2)/s02

# with robustbase
tmp2 <- lmrob(y~x1, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
tmp0 <- lmrob(y~1, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))

summary(tmp2)$r.squared

(s2 <- sum(rho(resid(tmp2)/tmp2$scale, cc=tmp2$control$tuning.psi)))
(s02 <- sum(rho(resid(tmp0)/tmp2$scale, cc=tmp2$control$tuning.psi)))
(s02 - s2)/s02


