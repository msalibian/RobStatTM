rm(list=ls())
library(robustbase)
library(pyinit)
library(quantreg)
source('lmrob2.R')
source('DCML.R')
source('refineSM.R')
# data(coleman)
# # Create simple data set with factors
# co2 <- coleman
# set.seed(123)
# co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
data(coleman)
co2 <- coleman
nadd <- 40
co2 <- rbind(as.matrix(co2), matrix(rnorm(ncol(co2)*nadd), nadd, ncol(co2)))
co2 <- as.data.frame(co2)
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=3, prob=.3)+1])

# Matias' version of SM+PY
set.seed(123)
(m2 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='PY', initial='SM'), data=co2))$coef
set.seed(123)
(m1 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='SS', initial='SM'), data=co2))$coef
set.seed(123)
(m0 <- lmrob(Y ~ .-1, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
c(m0$scale, m1$scale, m2$scale)


# Matias' version of SM+PY
set.seed(123)
(m2 <- lmrob2(Y ~ ., control=lmrob2.control(candidates='PY', initial='SM'), data=co2))$coef
set.seed(123)
(m1 <- lmrob2(Y ~ ., control=lmrob2.control(candidates='SS', initial='SM'), data=co2))$coef
set.seed(123)
(m0 <- lmrob(Y ~ ., control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
c(m0$scale, m1$scale, m2$scale)




# Check with Victor's new version of SM_PY

# Clean all
rm(list=ls())
# data(coleman)
# co2 <- coleman
# set.seed(123)
# co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
data(coleman)
co2 <- coleman
nadd <- 40
co2 <- rbind(as.matrix(co2), matrix(rnorm(ncol(co2)*nadd), nadd, ncol(co2)))
co2 <- as.data.frame(co2)
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=3, prob=.3)+1])



source('DCML_FINAL.R')
mf <- model.frame(Y ~ . -1, data=co2)
a <- splitFrame(mf, type='f')
Z <- a$x1 # x1 = factors, x2 = continuous, if there'
#s an intercept it's in x1
X <- a$x2
y <- co2$Y
options(warn=-1)
(smpy.victor <- SM_PY(y=y, X=X, Z=Z, intercept=FALSE))$coef
set.seed(123)
(m0 <- lmrob(Y ~ .-1, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
c(m0$scale, smpy.victor$scale)

# n <- nrow(Z)
# dee <- .5*(1-(length(smpy.victor$coef)/n))
# cc <-  find.tuning.chi(.5) # find.tuning.chi(dee)
# # smpy.victor
# smpy.victor$scale
# mscale(u=smpy.victor$resid, ep=1e-7, delta=dee)*1.547/find.tuning.chi(dee) 
# mad(smpy.victor$resid)
# 
# s <- mscale(u=smpy.victor$resid, ep=1e-7, delta=dee)*1.547/find.tuning.chi(dee) 
# mean(rho(smpy.victor$resid/s, cc=find.tuning.chi(dee)))
# dee
# 
