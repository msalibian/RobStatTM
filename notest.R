# test script

# R CMD INSTALL --preclean --clean robustbroli 

rm(list=ls())

# library(robustbroli)
library(robustbase)
library(pyinit)
library(quantreg)
source('R/lmrob2.R')
source('R/DCML.R')
data(coleman)
set.seed(123)
# source('R/lmrob2.R')

# ## Default for a very long time:
# m1 <- lmrob2(Y ~ ., data=coleman)
# X <- model.matrix(Y ~ ., data=coleman)
# m2 <- old.MMPY(X=X, y=coleman$Y, intercept=FALSE)
# # names(m2$coef) <- names(m2$coefficients) <- names(coef(m1))
# all.equal(coef(m1), coef(m2), check.attributes=FALSE)
# all.equal(m1$cov[,], m2$cov[,], check.attributes=FALSE)
# # dee <- .5*(1-(ncol(X)/nrow(X)))
# # m3 <- lmrob(Y ~ ., data=coleman, 
# #             control=lmrob.control(tuning.chi = 1.5477, bb = .5, tuning.psi = 3.4434))
# # all.equal(coef(m1), coef(m3), check.attributes=FALSE)
# 
# ## Default for a very long time:
# m1 <- lmrob2(Y ~ .-1, data=coleman)
# X <- model.matrix(Y ~ .-1, data=coleman)
# m2 <- old.MMPY(X=X, y=coleman$Y, intercept=FALSE)
# # names(m2$coef) <- names(m2$coefficients) <- names(coef(m1))
# all.equal(coef(m1), coef(m2), check.attributes=FALSE)
# all.equal(m1$cov[,], m2$cov[,], check.attributes=FALSE)


## Default for a very long time:
m2 <- lmrob2(Y ~ ., data=coleman)
m1 <- lmrob2(Y ~ ., control=lmrob2.control(candidates="SS"), data=coleman)
m0 <- lmrob(Y ~ ., data=coleman, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
# names(m2$coef) <- names(m2$coefficients) <- names(coef(m1))
all.equal(coef(m1), coef(m0), check.attributes=FALSE)
all.equal(m1$cov[,], m0$cov[,], check.attributes=FALSE)

## Default for a very long time:
m1 <- lmrob2(Y ~ . - 1, control=lmrob2.control(candidates="SS"), data=coleman)
m0 <- lmrob(Y ~ . -1, data=coleman, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))


# With factors!

co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

# Matias' version of SM+PY
(m1 <- lmrob2(Y ~ ., control=lmrob2.control(initial='SM'), data=co2))$coef
(m0 <- lmrob(Y ~ ., control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
c(m1$scale, m0$scale)

(m1 <- lmrob2(Y ~ . , control=lmrob2.control(initial='S', candidates='SS', prosac=.1), data=co2))$coef
(m0 <- lmrob(Y ~ ., control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), data=co2))$coef
c(m1$scale, m0$scale)




# co2 <- coleman
# set.seed(123)
# co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
# X <- model.matrix(Y ~ . , data=co2)
# 
# mf <- model.frame(Y ~ . , data=co2)
# (int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1))
# a <- splitFrame(mf, type='f') # type = c("f","fi", "fii"))
# head(a$x1) # x1 = factors, x2 = continuous, if there's an intercept it's in x1!
# head(a$x2)
# 
# mf <- model.frame(Y ~ .*educ , data=co2)
# a <- splitFrame(mf, type='fi') # type = c("f","fi", "fii"))
# head(a$x1) # x1 = factors + interactions
# head(a$x2) # x2 = continuous
# 
# mf <- model.frame(Y ~ . -1, data=co2)
# # is there an intercept?
# (int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1))
# a <- splitFrame(mf, type='f') # type = c("f","fi", "fii"))
# Z <- a$x1
# X <- a$x2
# old.SMPY(X=X, y=coleman$Y, Z=Z, intercept=FALSE)
# svd(Z)$d
# 
# # set.seed(123456)
# # x1 <- rnorm(50)
# # x2 <- rnorm(50)
# # x3 <- runif(50)
# # dum <- as.factor(LETTERS[rbinom(50, size=7, prob=.3)+1])
# # y <- rnorm(50)
# # co2 <- data.frame(x1=x1, x2=x2, x3=x3, dum=dum, Y=y)
# 
# # SMPY "by hand"

## Create data set
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])
mf <- model.frame(Y ~ .  , data=co2)
y <- co2$Y
outlmrob3 <- SMPY(mf=mf, y=y)


control <- lmrob.control(tuning.chi = 1.5477, bb = 0.5, tuning.psi = 3.4434)
(int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1))
a <- splitFrame(mf, type='f') # type = c("f","fi", "fii"))
Zorig <- Z <- a$x1 # x1 = factors, x2 = continuous, if there's an intercept it's in x1!
X <- a$x2
if(int.present) Z <- Zorig[, -1]
n <- nrow(X)
q <- ncol(Z)
p <- ncol(X)
gamma <- matrix(0, q, p)
#Eliminate Z from X by L1 regression , obtaining a matrix X1
options(warn=-1)
for( i in 1:p) gamma[,i] <- coef(rq(X[,i]~Z-1))
options(warn=0)
X1 <- X - Z %*% gamma
pp <- p + if(int.present) 1 else 0
dee <- .5*(1-(pp/n))
initial <- pyinit(intercept=int.present, X=X1, y=y, 
              deltaesc=dee, cc.scale=control$tuning.chi, 
              prosac=.5, clean.method="threshold", C.res = 2, prop=.2, 
              py.nit = 20, en.tol=1e-5)
betapy <- initial$initCoef[,1]
sspy <- initial$objF[1]
# if(int.present) {
#   y1 <- y - betapy[1]
#   betapy <- betapy[-1]
# }
uu <- list(coef=betapy, scale=sspy)
if(int.present) {
  out0 <- lmrob(y~X1, control=control,init=uu) } else {
  out0 <- lmrob(y~X1 - 1, control=control,init=uu)
  }
# out0 <- lmrob(y~X1 - 1, control=control,init=uu)
beta <- out0$coeff
# after eliminating the influence of X1 we make an L1 regression using as covariables Z
y1 <- resid(out0) # y-X1%*%beta
options(warn=-1)
fi <- coef(rq(y1~Z-1))
options(warn=0)
# retransform the coefficients to the original matrices X and Z
# oo <- NULL
# if(int.present) oo <- fi[1]
# tt <- gamma[hh1:hh3,]
# if(p==1) tt <- matrix(tt,q,p)
# beta00 <- c(oo,beta[hh1:hh2],fi[hh1:hh3]-  tt%*%beta[hh1:hh2])
if(int.present) {
beta00 <- c(beta, fi - gamma %*% beta[-1])
} else { beta00 <- c(beta, fi - gamma %*% beta) }
res <- as.vector( y1-Z%*%fi )
# XX <- cbind(X,Z)
nc <- ncol(X) + ncol(Z) + if(int.present) 1 else 0
dee <- .5*(1-(nc/n))
ss <- mscale(u=res, tol=1e-5, delta=dee, tuning.chi=control$tuning.chi)
uu <- list(coef=beta00, scale=ss)
# Compute the MMestimator using lmrob
# control <- lmrob.control(tuning.chi = 1.5477, bb = dee, tuning.psi = 3.4434)
XX <- cbind(X,Z)
if(int.present) {
outlmrob <- lmrob(y~XX,control=control,init=uu)
} else {
  outlmrob <- lmrob(y~XX-1,control=control,init=uu)
}
control$method <- 'M' 
control$cov <- ".vcov.w"
# # lmrob() does the above when is.list(init)==TRUE, in particular:
#
XX2 <- model.matrix(attr(mf, 'terms'), mf)
outlmrob2 <- lmrob.fit(XX2, y, control, init=uu, mf=mf)

# doesn't work when int.present is FALSE
# gives different results when int.present is TRUE
# old.SMPY=function(y,X,Z, intercept=TRUE)
# old.SMPY(y=y, X=X, Z=Z, intercept=int.present)


all.equal(coef(outlmrob), coef(outlmrob2), check.attributes=FALSE)
all.equal(coef(outlmrob), coef(outlmrob3), check.attributes=FALSE)
all.equal(coef(outlmrob2), coef(outlmrob3), check.attributes=FALSE)
all.equal(outlmrob$cov[,], outlmrob2$cov[,], check.attributes=FALSE)
all.equal(outlmrob3$cov[,], outlmrob2$cov[,], check.attributes=FALSE)





# pretty much the same as
control <- lmrob.control(tuning.chi = 1.5477, tuning.psi = 3.4434, init='M-S')
o2 <- lmrob(Y ~ . , control=control, data=co2)
o2$coef
o2$scale




# summary(lm(Y ~ . , data=co2))$sigma
