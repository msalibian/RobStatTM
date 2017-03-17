
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

# Y ~ N( 0, 0.5^2 )
summary(lm(y~., data=dat))$sigma

# Y ~ 7 * X_2 + N( 0, 0.5^2 )
summary(lm(y~., data=dat2))$sigma


# Y ~ N( 0, 0.5^2 )
mf <- model.frame(y ~ . , data=dat)
a <- splitFrame(mf, type='f')
Z <- a$x1 
X <- a$x2
y <- dat$y
options(warn=-1)
smpy.v <- SM_PY(y=y, X=X, Z=Z, intercept=FALSE)
smpy.v$scale

# Y ~ 7 * X_2 + N( 0, 0.5^2 )
mf <- model.frame(y ~ . , data=dat2)
a <- splitFrame(mf, type='f')
Z <- a$x1
X <- a$x2
y <- dat2$y
options(warn=-1)
smpy.v2 <- SM_PY(y=y, X=X, Z=Z, intercept=FALSE)
smpy.v2$scale


rm(list=ls())
library(robustbase)
library(pyinit)
library(quantreg)
source('lmrob2.R')
source('DCML.R')
source('refineSM.R')

load('problem500.Rdata')

summary(lm(y~., data=dat))$sigma
summary(lm(y~., data=dat2))$sigma

m2 <- lmrob2(y ~ ., control=lmrob2.control(candidates='PY', initial='SM', refine.PY=500), data=dat) #)$coef # MMPY
m2$scale

m22 <- lmrob2(y ~ ., control=lmrob2.control(candidates='PY', initial='SM', refine.PY=500), data=dat2) #)$coef # MMPY
m22$scale

# with sub-sampling candidates it works well
m1 <- lmrob2(y ~ ., control=lmrob2.control(candidates='SS', initial='SM', refine.PY=500), data=dat) #)$coef # MMPY
m1$scale

m12 <- lmrob2(y ~ ., control=lmrob2.control(candidates='SS', initial='SM', refine.PY=500), data=dat2) #)$coef # MMPY
m12$scale

# with PY candidates + "regular S"  it works well
m0 <- lmrob2(y ~ ., control=lmrob2.control(candidates='PY', initial='S', prosac=.02, refine.PY=500), data=dat) #)$coef # MMPY
m0$scale

m02 <- lmrob2(y ~ ., control=lmrob2.control(candidates='PY', initial='S', prosac=.02, refine.PY=500), data=dat2) #)$coef # MMPY
m02$scale

