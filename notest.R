# test script

# R CMD INSTALL --preclean --clean robustbroli 


# library(robustbroli)
library(robustbase)
library(pyinit)
data(coleman)
set.seed(123)
# source('R/lmrob2.R')

## Default for a very long time:
m1 <- lmrob2(Y ~ ., data=coleman)
X <- model.matrix(Y ~ ., data=coleman)
m2 <- old.MMPY(X=X, y=coleman$Y, intercept=FALSE)



m2 <- lmrob(Y ~ ., data=coleman)
S.init <- list(coef=coef(m2$init), scale=m2$scale)

m2 <- lmrob(Y ~ ., data=coleman, init=S.init, control=lmrob.control(trace.lev=1))


X <- model.matrix(Y ~ ., data=coleman)
MMPY(X=X, y=coleman$Y, intercept=FALSE)
