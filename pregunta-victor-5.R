# Clean all
rm(list=ls())
data(coleman)
co2 <- coleman
nadd <- 20
co2 <- rbind(as.matrix(co2), matrix(rnorm(ncol(co2)*nadd), nadd, ncol(co2)))
co2 <- as.data.frame(co2)
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=3, prob=.3)+1])

dim(co2)

source('DCML_FINAL.R')
mf <- model.frame(Y ~ . -1, data=co2)
a <- splitFrame(mf, type='f')
Z <- a$x1 # x1 = factors, x2 = continuous
X <- a$x2
y <- co2$Y
options(warn=-1)
smpy.victor <- SM_PY(y=y, X=X, Z=Z, intercept=FALSE)
smpy.victor$scale


n <- nrow(Z)
dee <- .5*(1-(length(smpy.victor$coef)/n))
# inconsistent scale
mscale(u=smpy.victor$resid, ep=1e-7, delta=dee)
# fixed
mscale(u=smpy.victor$resid, ep=1e-7, delta=dee)*1.547/2.560842
mad(smpy.victor$resid)


# compare PY-SM-initial with subsampling-SM-initial
rm(list=ls())
library(robustbase)
library(pyinit)
library(quantreg)
source('lmrob2.R')
source('DCML.R')
data(coleman)
co2 <- coleman
nadd <- 20
co2 <- rbind(as.matrix(co2), matrix(rnorm(ncol(co2)*nadd), nadd, ncol(co2)))
co2 <- as.data.frame(co2)
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=3, prob=.3)+1])

dim(co2)

n <- nrow(co2)
# Matias' version of SM+PY
set.seed(123)
# PY initial
m2 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='PY', initial='SM'), data=co2)
m2$scale
# subsampling initial
(m1 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='SS', initial='SM'), data=co2))$init$coef
m1$scale
