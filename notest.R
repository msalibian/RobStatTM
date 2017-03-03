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
(int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1))
a <- splitFrame(mf, type='f') # type = c("f","fi", "fii"))
head(a$x1) # x1 = factors, x2 = continuous, if there's an intercept it's in x1!
head(a$x2)

mf <- model.frame(Y ~ .*educ , data=co2)
a <- splitFrame(mf, type='fi') # type = c("f","fi", "fii"))
head(a$x1) # x1 = factors + interactions
head(a$x2) # x2 = continuous

mf <- model.frame(Y ~ . -1, data=co2)
# is there an intercept?
(int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1))
a <- splitFrame(mf, type='f') # type = c("f","fi", "fii"))
Z <- a$x1
X <- a$x2
old.SMPY(X=X, y=coleman$Y, Z=Z, intercept=FALSE)
svd(Z)$d

# set.seed(123456)
# x1 <- rnorm(50)
# x2 <- rnorm(50)
# x3 <- runif(50)
# dum <- as.factor(LETTERS[rbinom(50, size=7, prob=.3)+1])
# y <- rnorm(50)
# co2 <- data.frame(x1=x1, x2=x2, x3=x3, dum=dum, Y=y)

# SMPY "by hand"
co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

control <- lmrob.control(tuning.chi = 1.5477, bb = 0.5, tuning.psi = 3.4434)
mf <- model.frame(Y ~ . , data=co2)
y <- co2$Y
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
for( i in 1:p) gamma[,i] <- coef(rq(X[,i]~Z-1))
X1 <- X - Z %*% gamma
dee <- .5*(1-((p+1)/n))
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
fi <- coef(rq(y1~Z-1))
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
XX <- cbind(X,Z)
nc <- ncol(XX) + if(int.present) 1 else 0
dee <- .5*(1-(nc/n))
ss <- mscale(u=res, tol=1e-5, delta=dee, tuning.chi=control$tuning.chi)
uu <- list(coef=beta00, scale=ss)
#Compute the MMestimator using lmrob
control <- lmrob.control(tuning.chi = 1.5477, bb = dee, tuning.psi = 3.4434)
outlmrob <- lmrob(y~XX,control=control,init=uu)
outlmrob$coef
uu$scale

# pretty much the same as
control <- lmrob.control(tuning.chi = 1.5477, tuning.psi = 3.4434, init='M-S')
o2 <- lmrob(Y ~ . , control=control, data=co2)
o2$coef
o2$scale


# old.SMPY=function(y,X,Z, intercept=TRUE)
old.SMPY(y=y, X=X, Z=Z, intercept=int.present)


# summary(lm(Y ~ . , data=co2))$sigma
