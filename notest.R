# test script

# R CMD INSTALL --preclean --clean robustbroli 


# library(robustbroli)
library(robustbase)
library(pyinit)
library(quantreg)
source('R/lmrob2.R')
source('R/DCML.R')
data(coleman)
set.seed(123)
# source('R/lmrob2.R')

## Default for a very long time:
m1 <- lmrob2(Y ~ ., data=coleman)
X <- model.matrix(Y ~ ., data=coleman)
m2 <- old.MMPY(X=X, y=coleman$Y, intercept=FALSE)
# names(m2$coef) <- names(m2$coefficients) <- names(coef(m1))
all.equal(coef(m1), coef(m2), check.attributes=FALSE)
all.equal(m1$cov[,], m2$cov[,], check.attributes=FALSE)


m2 <- lmrob(Y ~ ., data=coleman)
S.init <- list(coef=coef(m2$init), scale=m2$scale)

m2 <- lmrob(Y ~ ., data=coleman, init=S.init, control=lmrob.control(trace.lev=1))


X <- model.matrix(Y ~ ., data=coleman)
MMPY(X=X, y=coleman$Y, intercept=FALSE)

co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
X <- model.matrix(Y ~ . , data=co2)

mf <- model.frame(Y ~ . , data=co2)
a <- splitFrame(mf, type='f') # type = c("f","fi", "fii"))
head(a$x1) # x1 = factors, x2 = continuous
head(a$x2)

mf <- model.frame(Y ~ .*educ , data=co2)
a <- splitFrame(mf, type='fi') # type = c("f","fi", "fii"))
head(a$x1) # x1 = factors + interactions
head(a$x2) # x2 = continuous

mf <- model.frame(Y ~ . -1, data=co2)
a <- splitFrame(mf, type='f') # type = c("f","fi", "fii"))
Z <- a$x1
X <- a$x2
old.SMPY(X=X, y=coleman$Y, Z=Z, intercept=FALSE)


