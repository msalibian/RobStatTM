# test script

# R CMD INSTALL --preclean --clean robustbroli 

rm(list=ls())

# library(robustbroli)
library(robustbase)
library(pyinit)
library(quantreg)
source('R/lmrob2.R')
source('R/DCML.R')
source('R/refineSM.R')

data(coleman)

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
m2 <- lmrob2(Y ~ ., data=coleman) # MMPY
# m1 <- lmrob2(Y ~ ., control=lmrob2.control(candidates="SS", initial='S'), data=coleman)
m0 <- lmrob(Y ~ ., data=coleman, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
c(m0$scale, m2$scale)
all.equal(coef(m1), coef(m0), check.attributes=FALSE)
all.equal(coef(m2), coef(m0), check.attributes=FALSE)
all.equal(m1$cov[,], m0$cov[,], check.attributes=FALSE)
all.equal(m2$cov[,], m0$cov[,], check.attributes=FALSE)

## Default for a very long time:
m2 <- lmrob2(Y ~ . - 1 , data=coleman) # MMPY
m1 <- lmrob2(Y ~ . - 1 , control=lmrob2.control(candidates="SS", initial='S'), data=coleman)
m0 <- lmrob(Y ~ . - 1 , data=coleman, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
c(m0$scale, m1$scale, m2$scale)
all.equal(coef(m1), coef(m0), check.attributes=FALSE)
all.equal(coef(m2), coef(m0), check.attributes=FALSE)
all.equal(m1$cov[,], m0$cov[,], check.attributes=FALSE)
all.equal(m2$cov[,], m0$cov[,], check.attributes=FALSE)

round(summary(lm(Y~.-1, data=coleman))$cov.unscaled * summary(lm(Y~.-1, data=coleman))$sigma^2, 5)




rm(list=ls())

# library(robustbroli)
library(robustbase)
library(pyinit)
library(quantreg)
source('R/lmrob2.R')
source('R/DCML.R')
source('R/refineSM.R')

set.seed(123)
n <- 50
x1 <- rnorm(n)
x2 <- rexp(n, rate=.3)
x3 <- runif(n)
y <- 2*x1 - x2 + 1 + rnorm(n, sd=.5)
ou <- rbinom(n, size=1, prob=.05)
y[ou==1] <- 8
x2[ou==1] <- 2
coef(m2 <- lmrob2(y~x1+x2+x3))
coef(m0 <- lm(y~x1+x2+x3))
coef(m1 <- lmrob(y~x1+x2+x3))

# > m2 <- lmrob2(y~x1+x2+x3, control=lmrob2.control(candidates='PY'))
# (Intercept)         x1         x2        x3
# (Intercept)    7.680582  1.2390906 -1.0037035 -6.585118
# x1             1.239091  1.2299773 -0.1744918 -1.257578
# x2            -1.003703 -0.1744918  0.2207745  0.346452
# x3            -6.585118 -1.2575779  0.3464520 11.460809
# 


rm(list=ls())
source('R/DCML_FINAL_NEW.R')
set.seed(123)
x1 <- rnorm(20)
x2 <- rexp(20, rate=.3)
x3 <- runif(20)
y <- rnorm(20, sd=3.5)
mf <- model.frame(y ~ x1 + x2 + x3)
int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1)
a <- splitFrame(mf, type='f')
Z <- a$x1 
if(int.present) Z <- Z[, -1]
X <- a$x2
options(warn=-1)
smpy.v <- MMPY(y=y, X=X, intercept=FALSE)
DCML_FINAL(X=X,y=y, outMM=smpy.v, intercept=FALSE)$cov
  

round(summary(lm(y~x1+x2))$cov.unscaled * summary(lm(y~x1+x2))$sigma^2, 5)




# With factors!
data(breslow.dat, package='robust')
data(stack.dat, package='robust')
st2 <- stack.dat
for(j in 1:length(st2)) st2[[j]] <- as.double(st2[[j]])
m2 <- lmrob2(Loss ~ ., control=lmrob2.control(corr.b=TRUE), data = st2)
m1 <- lmrob2(Loss ~ ., control=lmrob2.control(candidates="SS", initial='S', corr.b=FALSE), data=st2)
m0 <- lmrob(Loss ~ ., data=st2, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
c(m0$scale, m1$scale, m2$scale)
all.equal(coef(m1), coef(m0), check.attributes=FALSE)
all.equal(coef(m2), coef(m0), check.attributes=FALSE)
all.equal(m1$cov[,], m0$cov[,], check.attributes=FALSE)
all.equal(m2$cov[,], m0$cov[,], check.attributes=FALSE)


sqrt(summary(lm(Loss~., data=st2))$cov[1,1]) * summary(lm(Loss~.,data=st2))$sigma




# synthetic data
n <- 100
p <- 5
lf <- 4
set.seed(123456)
x1 <- rnorm(n)
x2 <- rexp(n)
x3 <- runif(n)
x4 <- runif(n) # rbinom(n, size=4, prob=.7)
f1 <- as.factor(LETTERS[rbinom(n, size=lf, prob=.5) + 1])
f2 <- as.factor(rbinom(n, size=lf, prob=.5) + 1)
set.seed(77)
y2 <- rnorm(n, sd=.5) + 7*x2 # 2*(as.numeric(f2) - 1)
dat2 <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, f1=f1, f2=f2, y=y2)
# y <- rnorm(n, sd=.5) + 0 # 2*(as.numeric(f2) - 1)
# dat <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, f1=f1, f2=f2, y=y)


m2 <- lmrob2(y ~ ., control=lmrob2.control(candidates='PY', initial='SM', refine.PY=500), data=dat2) #)$coef # MMPY
m2$scale

m22 <- lmrob2(y ~ ., control=lmrob2.control(candidates='PY', initial='SM', refine.PY=500), data=dat2) #)$coef # MMPY
m22$scale


set.seed(123)
(m1 <- lmrob2(y ~ ., control=lmrob2.control(candidates='SS', initial='SM'), data=dat))$coef # MMPY
set.seed(123)
(m12 <- lmrob2(y ~ ., control=lmrob2.control(candidates='SS', initial='SM'), data=dat2))$coef # MMPY

set.seed(123)
(m0 <- lmrob(y ~ ., control=lmrob.control(tuning.psi=3.4434, subsampling='simple'),  init='M-S', data=dat))$coef # MMPY
set.seed(123)
(m02 <- lmrob(y ~ ., control=lmrob.control(tuning.psi=3.4434, subsampling='simple'),  init='M-S', data=dat2))$coef # MMPY

c(m0$scale, m1$scale, m2$scale)
c(m02$scale, m12$scale, m22$scale)




all.equal(coef(m1), coef(m0), check.attributes=FALSE)
all.equal(coef(m2), coef(m0), check.attributes=FALSE)


X <- model.matrix(y ~ ., data=dat)
m2 <- old.MMPY(X=X, y=dat$y, intercept=FALSE)



m1 <- lmrob2(y ~ .  , control=lmrob2.control(candidates="SS", initial='S'), data=dat)
m0 <- lmrob(y ~ .  , control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), data=dat)
c(m0$scale, m1$scale, m2$scale)
all.equal(coef(m1), coef(m0), check.attributes=FALSE)
all.equal(coef(m2), coef(m0), check.attributes=FALSE)
all.equal(m1$cov[,], m0$cov[,], check.attributes=FALSE)
all.equal(m2$cov[,], m0$cov[,], check.attributes=FALSE)




# bombs! a single 2-level factor + intercept is trouble!
br2 <- breslow.dat
for(j in 1:length(br2)) if(class(br2[[j]])=='integer') br2[[j]] <- as.double(br2[[j]])
m2 <- lmrob2(Ysum ~ Base + Age + Trt, control=lmrob2.control(initial='SM', corr.b=TRUE), data = br2)
m1 <- lmrob2(Ysum ~ Base + Age + Trt, control=lmrob2.control(candidates="SS", initial='S', corr.b=FALSE), data = br2)
m0 <- lmrob(Ysum ~ Base + Age + Trt, data=br2, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
c(m0$scale, m1$scale, m2$scale)
all.equal(coef(m1), coef(m0), check.attributes=FALSE)
all.equal(coef(m2), coef(m0), check.attributes=FALSE)
all.equal(m1$cov[,], m0$cov[,], check.attributes=FALSE)
all.equal(m2$cov[,], m0$cov[,], check.attributes=FALSE)


co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

# Matias' version of SM+PY
(m2 <- lmrob2(Y ~ ., control=lmrob2.control(candidates='PY', initial='SM'), data=co2))$coef
(m1 <- lmrob2(Y ~ ., control=lmrob2.control(candidates='SS', initial='SM'), data=co2))$coef
(m0 <- lmrob(Y ~ ., control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
c(m0$scale, m1$scale, m2$scale)

(m2 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='PY', initial='SM'), data=co2))$coef
(m1 <- lmrob2(Y ~ .-1, control=lmrob2.control(candidates='SS', initial='SM'), data=co2))$coef
(m0 <- lmrob(Y ~ .-1, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
c(m0$scale, m1$scale, m2$scale)

rm(list=ls())
library(robustbase)
library(pyinit)
library(quantreg)
data(coleman)
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

source('R/DCML_FINAL.R')
mf <- model.frame(Y ~ . -1, data=co2)
a <- splitFrame(mf, type='f') 
Z <- a$x1 # x1 = factors, x2 = continuous, if there's an intercept it's in x1
X <- a$x2
y <- co2$Y
(smpy.victor <- SM_PY(y=y, X=X, Z=Z, intercept=TRUE))
SMPY(mf=mf, y=y, split=a)





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
