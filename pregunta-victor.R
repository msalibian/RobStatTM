
library(robustbase)
library(pyinit)
library(quantreg)

# delete everything
rm(list=ls())

# create a new data set
# with an artificial factor with 3 levels, called "educ"
# the new data set is called "co2" 
data(coleman)
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

# Matias' version of SM+PY
source('R/DCML.R')
# the argument is the model matrix plus the response variable
mf <- model.frame(Y ~ .  , data=co2)
y <- co2$Y
smpy.matias <- SMPY(mf=mf, y=y)
print(smpy.matias$coef)
# (Intercept)     salaryP    fatherWc     sstatus   teacherSc   motherLev       educB       educC 
# 30.52985182 -1.74514053  0.09131376  0.66582527  1.20618825 -4.25523553 -0.50495552 -0.30661485 

# compare with current implementation of the MS estimator in robustbase::lmrob()
rm(list=ls())
# create a new data set
# with an artificial factor with 3 levels, called "educ"
# the new data set is called "co2" 
data(coleman)
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
my.ctrl <- lmrob.control(tuning.chi = 1.5477, bb = 0.5, tuning.psi = 3.4434, init='M-S')
sm.current <- lmrob(Y ~ ., data=co2, control=my.ctrl)
print(sm.current$coef)
# (Intercept)     salaryP    fatherWc     sstatus   teacherSc   motherLev       educB       educC 
# 30.57076737 -1.72840071  0.09139785  0.66547912  1.20350537 -4.25695985 -0.52139165 -0.32012493 

# now try to run Victor's original code
rm(list=ls())
# create a new data set
# with an artificial factor with 3 levels, called "educ"
# the new data set is called "co2" 
data(coleman)
set.seed(123)
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
source('R/DCML_FINAL.R')

# we build the design matrices (X for cont, Z for factors)
mf <- model.frame(Y ~ .  , data=co2)
a <- splitFrame(mf, type='f') 
Z <- a$x1 # x1 = factors, x2 = continuous, if there's an intercept it's in x1
X <- a$x2
y <- co2$Y
# I think we need to remove the intercept from Z?
Z <- Z[, -1]
(smpy.victor <- SM_PY(y=y, X=X, Z=Z, intercept=TRUE))
# (Intercept)   XXsalaryP  XXfatherWc   XXsstatus XXteacherSc XXmotherLev     XXeducB     XXeducC 
# 30.52589694 -1.74676055  0.09130605  0.66585957  1.20644662 -4.25507104 -0.50330768 -0.30525631 


