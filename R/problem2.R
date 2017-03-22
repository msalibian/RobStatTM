
# How the data were created
# 
# n <- 500
# p <- 5
# lf <- 4
# set.seed(123456)
# x1 <- rnorm(n)
# x2 <- rexp(n)
# x3 <- runif(n)
# x4 <- runif(n) # rbinom(n, size=4, prob=.7)
# f1 <- as.factor(LETTERS[rbinom(n, size=lf, prob=.5) + 1])
# f2 <- as.factor(rbinom(n, size=lf, prob=.5) + 1)
# set.seed(77)
# y2 <- rnorm(n, sd=.5) + 7*x2 # 2*(as.numeric(f2) - 1)
# dat2 <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, f1=f1, f2=f2, y=y2)
# # y <- rnorm(n, sd=.5) + 0 # 2*(as.numeric(f2) - 1)
# # dat <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, f1=f1, f2=f2, y=y)
# save(list=c('dat', 'dat2'), file='problem500.RData')





# It does not seem to be a problem of my version of the code
rm(list=ls())
source('DCML_FINAL_NEW.R')
load('problem500.RData')
# otra transformacion, ningun problema si dejo el intercept
# # horrible si no
dat2$y <- dat$y + 8 * dat$x1 - 10 * dat$x2

# Y ~ N( 0, 0.5^2 )
mf <- model.frame(y ~ .-1, data=dat)
int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1)
a <- splitFrame(mf, type='f')
Z <- a$x1 
if(int.present) Z <- Z[, -1]
X <- a$x2
y <- dat$y
options(warn=-1)
smpy.v <- SM_PY(y=y, X=X, Z=Z, intercept=int.present)
smpy.v$scale

# Y ~ 7 * X_2 + N( 0, 0.5^2 )
mf <- model.frame(y ~ .-1, data=dat2)
int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1)
a <- splitFrame(mf, type='f')
Z <- a$x1
if(int.present) Z <- Z[, -1]
X <- a$x2
y <- dat2$y
options(warn=-1)
smpy.v2 <- SM_PY(y=y, X=X, Z=Z, intercept=int.present)
smpy.v2$scale

cbind(round(smpy.v$coef, 4),  round(smpy.v2$coef, 4))


# > round(cbind(coef(m2), coef(m22)), 4)
# [,1]    [,2]
# x1   0.0256  0.0183
# x2  -0.0108  6.9997
# x3   0.0417  0.0532
# x4   0.0717  0.0760
# f1A -0.1395 -0.1366
# f1B  0.0074 -0.0240
# f1C  0.0569  0.0363
# f1D -0.0290 -0.0637
# f1E -0.0727 -0.1503
# f22 -0.0382 -0.0429
# f23 -0.0008  0.0088
# f24  0.0234  0.0138
# f25 -0.0978 -0.0965
# > 

  mscale(smpy.v2$res, 1e-7, .487)
mscale(smpy.v$res, 1e-7, .487)



# With Matias' code (problem is less severe, but it's there)
rm(list=ls())
library(robustbase)
library(pyinit)
library(quantreg)
source('lmrob2.R')
source('DCML.R')
source('refineSM.R')

load('problem500.RData')
dat2$y <- dat$y - 7  * dat$x2 # + 8 * dat$x1
# otra transformacion, ningun problema si dejo el intercept
# # horrible si no

lm(formula=the.f, data=dat)
lm(formula=the.f, data=dat2)

the.f <- formula(y ~ . -1)
# PY candidates + SM
# Y ~ N( 0, 0.5^2 )
m2 <- lmrob2(formula=the.f, control=lmrob2.control(candidates='PY', initial='SM', refine.PY=500), data=dat) #)$coef # MMPY
m2$scale
# [1] 0.496774

# PY candidates + SM
# Y ~ 7 * X_2 + N( 0, 0.5^2 )
m22 <- lmrob2(formula=the.f, control=lmrob2.control(candidates='PY', initial='SM', refine.PY=500), data=dat2) #)$coef # MMPY
m22$scale
# [1] 0.5033762

round(c(m2$scale, m22$scale), 4)
round(cbind(coef(m2), coef(m22)), 4)




# with sub-sampling candidates + SM it works well
# Y ~ N( 0, 0.5^2 )
set.seed(123)
m1 <- lmrob2(formula=the.f, control=lmrob2.control(candidates='SS', initial='SM', refine.PY=500), data=dat) #)$coef # MMPY
m1$scale
# [1] 0.4870675

# with sub-sampling candidates + SM it works well
# Y ~ 7 * X_2 + N( 0, 0.5^2 )
set.seed(123)
m12 <- lmrob2(formula=the.f, control=lmrob2.control(candidates='SS', initial='SM', refine.PY=500), data=dat2) #)$coef # MMPY
m12$scale
# [1] 0.4870614

# with PY candidates + "regular S"  it works well
set.seed(123)
m0 <- lmrob2(formula=the.f, control=lmrob2.control(candidates='PY', initial='S', prosac=.02, refine.PY=500), data=dat) #)$coef # MMPY
m0$scale
# [1] 0.4845817

set.seed(123)
m02 <- lmrob2(formula=the.f, control=lmrob2.control(candidates='PY', initial='S', prosac=.02, refine.PY=500), data=dat2) #)$coef # MMPY
m02$scale
# [1] 0.4845817
# 

round(c(m0$scale, m02$scale, m1$scale, m12$scale, m2$scale, m22$scale), 4)
round(cbind(coef(m0), coef(m02), coef(m1), coef(m12), coef(m2), coef(m22)), 4)
