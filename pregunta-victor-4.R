# Clean all
rm(list=ls())
data(coleman)
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

source('DCML_FINAL.R')
mf <- model.frame(Y ~ . -1, data=co2)
a <- splitFrame(mf, type='f')
Z <- a$x1 # x1 = factors, x2 = continuous
X <- a$x2
y <- co2$Y
options(warn=-1)
smpy.victor <- SM_PY(y=y, X=X, Z=Z, intercept=FALSE)
# Initial residuals
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -23.300  -4.154   0.000  11.150   7.778  95.980 
# Initial coeff
# X1salaryP  X1fatherWc   X1sstatus X1teacherSc X1motherLev                                     
# 16.586219    1.206977   -3.443043   18.903357    8.772294 -593.844926 -592.109992 -564.086065 

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
# Create simple data set with factors
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
n <- nrow(co2)
# Matias' version of SM+PY
set.seed(123)
# PY initial
m2 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='PY', initial='SM'), data=co2)
# subsampling initial
(m1 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='SS', initial='SM'), data=co2))$init$coef
