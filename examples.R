# test script

# R CMD INSTALL --preclean --clean robustbroli

library(RobStatTM)
set.seed(123)
x <- rnorm(5000, sd=2)
mscale(u=x, delta=.5, tuning.chi=1.547645, family='bisquare')
mscale(u=x, delta=.5, tuning.chi=lmrobdet.control()$tuning.chi, family='bisquare')
mscale(u=x, delta=.5, tuning.chi=lmrobdet.control(family='modified.optimal')$tuning.chi, family='modified.optimal')
mscale(u=x, delta=.01, tuning.chi=lmrobdet.control(bb=.01)$tuning.chi, family='bisquare')
mscale(u=x, delta=.01, tuning.chi=lmrobdet.control(family='modified.optimal', bb=.01)$tuning.chi, family='modified.optimal')
mscale(u=x, delta=.01, tuning.chi=lmrobdet.control(family='modified.optimal', efficiency=.95, bb=.01)$tuning.chi, family='modified.optimal')
sd(x)
t.chi <- lmrobdet.control(bb=.25)$tuning.chi
tmp <- mscale(u=x, delta=.25, tuning.chi=t.chi, family='bisquare')
mean( rho(u=x/tmp, family='bisquare', cc=t.chi) )

t.chi <- lmrobdet.control(bb=.25, family='modified.optimal')$tuning.chi
tmp <- mscale(u=x, delta=.25, tuning.chi=t.chi, family='modified.optimal')
mean( rho(u=x/tmp, family='modified.optimal', cc=t.chi) )

t.chi <- lmrobdet.control(family='modified.optimal')$tuning.chi
tmp <- mscale(u=x, delta=.5, tuning.chi=t.chi, family='modified.optimal')
mean( rho(u=x/tmp, family='modified.optimal', cc=t.chi) )

t2.chi <- lmrobdet.control(family='modified.optimal', efficiency=.95)$tuning.chi
tmp2 <- mscale(u=x, delta=.5, tuning.chi=t2.chi, family='modified.optimal')
mean( rho(u=x/tmp2, family='modified.optimal', cc=t2.chi) )



library(RobStatTM)
data(coleman, package='robustbase')
m2 <- lmrobdetMM(Y ~ ., data=coleman)
m4 <- lmrobdetMM(Y ~ ., data=coleman, control=lmrobdet.control(family='modified.optimal'))
m5 <- lmrobdetMM(Y ~ ., data=coleman, control=lmrobdet.control(family='optimal'))
m1 <- lmrobdetDCML(Y ~ ., data=coleman)
m3 <- lmrobdetDCML(Y ~ ., data=coleman, control=lmrobdet.control(efficiency=.999))
m6 <- lmrobdetDCML(Y ~ ., data=coleman, control=lmrobdet.control(family='modified.optimal', efficiency=.999))
m0 <- lm(Y ~ ., data=coleman)

coef(m0)
coef(m6)
coef(m3)
coef(m2)
coef(m4)
coef(m1)

# m0 <- lmrob(Y ~ ., data=coleman, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
step.lmrobdetMM(m2)
step.lmrobdetMM(m2, whole.path=TRUE)
m3 <- lmrobdetMM(Y~sstatus, data=coleman)
step.lmrobdetMM(m3)
m4 <- lmrobdetMM(Y ~ 1, data=coleman)

m2 <- lmrobdetMM(Y ~ ., data=coleman)
m2.ls <- lmrobdetMM(Y ~ ., data=coleman, control=lmrobdet.control(efficiency = .999))
step.lmrobdetMM(m2)
step.lmrobdetMM(m2.ls)
step(lm(Y~., data=coleman), direction='backward') #, scale=2.074296)

data(coleman, package='robustbase')
m2 <- lmrobdetMM(Y ~ ., data=coleman)

# lmrobdet.RFPE(m2)
object <- m2
p <- length(object$coef)
scale <- object$scale
res <- residuals(object)/scale
a <- mean(rho(u=res, family=object$control$tuning.psi))
b <- p * mean(rhoprime(u=res, family=object$control$tuning.psi)^2) #, standardize=TRUE)^2)
d <- mean(rhoprime2(u=res, family=object$control$tuning.psi)) #, standardize=TRUE))
tun <- object$control$tuning.psi$cc
a2 <- mean(rho(u=res, family=object$control$tuning.psi, standardize=TRUE))
b2 <- p * mean(rhoprime(u=res, family=object$control$tuning.psi, standardize=TRUE)^2)
d2 <- mean(rhoprime2(u=res, family=object$control$tuning.psi, standardize=TRUE))
((a+b/d*6/tun^2))
((a2+b2/d2))




## The rho functions and their scaling
tt <- seq(-6, 6, length=500)
uu <- rho(u=tt, family=bisquare(.9), standardize=TRUE)
plot(tt, uu, type='l')

tt <- seq(-6, 6, length=500)
uu <- rhoprime(u=tt, family=bisquare(.9))
plot(tt, uu, type='l')

x0 <- .7
eff <- .85
rhoprime(u=x0, family='bisquare', cc=1.54764, standardize=TRUE)
ep <- 1e-4
(rho(u=x0+ep, family='bisquare', cc=1.54764, standardize=TRUE) -
    rho(u=x0-ep, family='bisquare', cc=1.54764, standardize=TRUE) ) / (2*ep)

x0 <- .7
eff <- .85
rhoprime(u=x0, family='bisquare', cc=1.54764, standardize=FALSE)
ep <- 1e-4
(rho(u=x0+ep, family='bisquare', cc=1.54764, standardize=FALSE) -
    rho(u=x0-ep, family='bisquare', cc=1.54764, standardize=FALSE) ) / (2*ep)


rhoprime2(u=x0, family='bisquare', cc=1.54764, standardize=TRUE)
ep <- 1e-4
(rhoprime(u=x0+ep, family='bisquare', cc=1.54764, standardize=TRUE) -
    rhoprime(u=x0-ep, family='bisquare', cc=1.54764, standardize=TRUE) ) / (2*ep)

rhoprime2(u=x0, family='bisquare', cc=1.54764, standardize=FALSE)
ep <- 1e-4
(rhoprime(u=x0+ep, family='bisquare', cc=1.54764, standardize=FALSE) -
    rhoprime(u=x0-ep, family='bisquare', cc=1.54764, standardize=FALSE) ) / (2*ep)



coef(m2)
coef(m1)
coef(m0)
coef(lm(Y~., data=coleman))

summary(lm(Y~., data=coleman))


summary(m2)
summary(m1)
summary(m0)

m2.null <- lmrobdet(Y ~ . - sstatus -fatherWc, data=coleman)
rob.linear.test(m2, m2.null)

set.seed(123)
x1 <- rnorm(50)
x2 <- rnorm(50)
y <- rnorm(50, sd=.1) + x1 - 2*x2
das <- data.frame(x1=x1, x2=x2, y=y)
tmp <- lmrobdetMM(y ~ . , data=das)
summary(tmp)
summary(lm(y~., data=das))


# fixed designs
set.seed(123)
n <- 50
x1 <- factor(letters[rbinom(n, size=2, prob=.3)+1])
x2 <- factor(c('L1', 'L2', 'L3')[3-rbinom(n, size=2, prob=.3)])
y <- rnorm(n, sd=.7) + 2*(as.numeric(x1)-1) - 3*(as.numeric(x2)-1)
das <- data.frame(x1=x1, x2=x2, y=y)
tmp <- lmrobM(y ~ . , data=das)

# library(quantreg)
# tmp2 <- lmrobM.bak(y~x1 + x2)
# coef(tmp)
# coef(tmp2)

## Default for a very long time:
m2 <- lmrobdet(Y ~ . - 1 , data=coleman) # MMPY
m0 <- lmrob(Y ~ . - 1 , data=coleman, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))

plot(m2)

data(stack.dat, package='robust')
# pyinit can't take integers?
st2 <- stack.dat
for(j in 1:length(st2)) st2[[j]] <- as.double(st2[[j]])
m2 <- lmrobdet(Loss ~ ., control=lmrobdet.control(corr.b=TRUE), data = st2)
m0 <- lmrob(Loss ~ ., data=st2, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
m2$coef
m0$coef
plot(m2)
plot(m0)



# synthetic data
n <- 50
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
ou <- rbinom(n, size=1, prob=.1)
y2[ou==1] <- rnorm(sum(ou), mean=-7, sd=.2)
x3[ou==1] <- runif(sum(ou), min=-1.2, max=-.8)
dat2 <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, f1=f1, f2=f2, y=y2)
# y <- rnorm(n, sd=.5) + 0 # 2*(as.numeric(f2) - 1)
# dat <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, f1=f1, f2=f2, y=y)


m2 <- lmrobdet(y ~ ., control=lmrobdet.control(initial='SM', refine.PY=10), data=dat2) #)$coef # MMPY

# MMPY with factors needs very small prosac, it's not robust anymore!
m22 <- lmrobdet(y ~ ., control=lmrobdet.control(initial='S', refine.PY=10, prosac=.02), data=dat2) #)$coef # MMPY

m1 <- lmrob(y ~ ., control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=dat2)

m0 <- lm(y~., data=dat2)
c(m2$scale, m22$scale, m1$scale)


coef(m0)
coef(m1)
coef(m2)
coef(m22)


data(breslow.dat, package='robust')
br2 <- breslow.dat
for(j in 1:length(br2)) if(class(br2[[j]])=='integer') br2[[j]] <- as.double(br2[[j]])
m2 <- lmrob2(Ysum ~ Base + Age + Trt, control=lmrob2.control(initial='SM'), data = br2)
m0 <- lmrob(Ysum ~ Base + Age + Trt, data=br2, init='M-S', control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
c(m0$scale, m2$scale)
coef(m2)
coef(m0)
coef(lm(Ysum ~ Base + Age + Trt, data=br2))
m2 <- lmrob2(Ysum ~ Base + Age + Trt-1, control=lmrob2.control(initial='SM'), data = br2)
m0 <- lmrob(Ysum ~ Base + Age + Trt-1, data=br2, init='M-S', control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
c(m0$scale, m2$scale)
coef(m2)
coef(m0)
coef(lm(Ysum ~ Base + Age + Trt-1, data=br2))



co2 <- coleman
set.seed(123)
co2$educ <- as.factor(LETTERS[rbinom(nrow(co2), size=2, prob=.3)+1])

# Matias' version of SM+PY
(m2 <- lmrobdet(Y ~ ., control=lmrobdet.control(initial='SM'), data=co2))$coef
(m0 <- lmrob(Y ~ ., control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
c(m0$scale, m2$scale)
coef(m2)
coef(m0)
coef(lm(Y~., data=co2))

(m2 <- lmrobdet(Y ~ .-1, control=lmrobdet.control(initial='SM'), data=co2))$coef
(m0 <- lmrob(Y ~ .-1, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'), init='M-S', data=co2))$coef
c(m0$scale, m2$scale)
coef(m2)
coef(m0)
coef(lm(Y~.-1, data=co2))




library(RobStatTM)
xx <- seq(-3, 3, length=500)
par(mfrow=c(2,2))
plot(xx, rho(u=xx, family='bisquare', cc=1.547, standardize=TRUE), type='l', main='Rho - Standarize: TRUE', xlab='', ylab='')

plot(xx, rho(u=xx, family='bisquare', cc=1.547, standardize=FALSE), type='l', main='Rho - Standarize: FALSE', xlab='', ylab='')

plot(xx, rhoprime(u=xx, family='bisquare', cc=1.547, standardize=TRUE), type='l', main='Rhoprime - Standarize: TRUE', xlab='', ylab='')

plot(xx, rhoprime(u=xx, family='bisquare', cc=1.547, standardize=FALSE), type='l', main='Rhoprime - Standarize: FALSE', xlab='', ylab='')


library(RobStatTM)
effs <- c(.85, .9, .95, .99, .999)
for(j in 1:length(effs))
  print(round(c(effs[j], (c(optimal(effs[j])))[1:3]), 8))
#                    a      lower      upper 
# 0.85000000 0.04358443 0.10991187 2.50257316 
#                    a      lower      upper 
# 0.90000000 0.02789679 0.07009889 2.70366981 
#                    a      lower      upper 
# 0.95000000 0.01317785 0.03305002 3.00333205 
#                    a      lower      upper 
# 0.99000000 0.00243798 0.00611122 3.56932005  
#                    a      lower      upper 
# 0.99900000 0.00024659 0.00061810 4.20099435 
# 


set.seed(123)
n <- 50
p <- 2
x0 <- matrix(rnorm(n * p), n, p)
y <- rbinom(n, size = 1, prob=.75)
# library(robustbase)
# library(RobStatTM)
tmp <- BYlogreg(x0=x0, y=y)
#tmp2 <- BYlogreg(x0=x0, y=y)
tmp.w <-  WBYlogreg(x0=x0, y=y)
tmp.wml <- WMLlogreg(x0=x0, y=y)


data(vaso, package='robustbase')
x <- model.matrix(Y ~ Volume + Rate, data=vaso)
tmp <- BYlogreg(x0=x, y=vaso$Y, intercept=FALSE)
tmp.w <-  WBYlogreg(x0=x, y=vaso$Y, intercept=FALSE)
tmp.wml <- WMLlogreg(x0=x, y=vaso$Y, intercept=FALSE)




