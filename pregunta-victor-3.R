rm(list=ls())
library(robustbase)
library(pyinit)
library(quantreg)
source('R/lmrob2.R')
source('R/DCML.R')
data(coleman)
# Create simple data set with factors
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
n <- nrow(co2)
# Matias' version of SM+PY
set.seed(123)
(m2 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='PY', initial='SM'), data=co2))$coef
dee <- .5*(1-(length(m2$coef)/n))
cc <-  find.tuning.chi(.5) # find.tuning.chi(dee)
mscale(u=m2$resid, tol=1e-7, delta=dee, tuning.chi = cc)
m2$scale
# m2$resid
set.seed(123)
(m1 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='SS', initial='SM'), data=co2))$coef
m1$scale
mscale(u=m1$resid, tol=1e-7, delta=dee, tuning.chi = cc)
# m1$resid 
set.seed(123)
(m0 <- lmrob(Y ~ .-1, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
m0$scale
mscale(u=m0$resid, tol=1e-7, delta=dee, tuning.chi = cc)
# m0resid

# Check with Victor's new version of SM_PY

# Clean all
rm(list=ls())
data(coleman)
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

source('R/DCML_FINAL.R')
mf <- model.frame(Y ~ . -1, data=co2)
a <- splitFrame(mf, type='f')
Z <- a$x1 # x1 = factors, x2 = continuous, if there'
#s an intercept it's in x1
X <- a$x2
y <- co2$Y
(smpy.victor <- SM_PY(y=y, X=X, Z=Z, intercept=FALSE))
n <- nrow(Z)
dee <- .5*(1-(length(smpy.victor$coef)/n))
cc <-  find.tuning.chi(dee)
# smpy.victor
smpy.victor$scale
mscale(u=smpy.victor$resid, ep=1e-7, delta=.3)*find.tuning.chi(.5)/cc
mad(smpy.victor$resid)


