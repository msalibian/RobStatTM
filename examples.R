# test script

# R CMD INSTALL --preclean --clean robustbroli

library(RobustStatistics)
data(coleman, package='robustbase')
m2 <- lmrobdet(Y ~ ., data=coleman)
m1 <- lmrobdetDCML(Y ~ ., data=coleman)
# m0 <- lmrob(Y ~ ., data=coleman, control=lmrob.control(tuning.psi=3.4434, subsampling='simple'))
step.lmrobdet(m2)
step.lmrobdet(m2, whole.path=TRUE)
m3 <- lmrobdet(Y~sstatus, data=coleman)
step.lmrobdet(m3)
m4 <- lmrobdet(Y ~ 1, data=coleman)

m2 <- lmrobdet(Y ~ ., data=coleman)
m2.ls <- lmrobdet(Y ~ ., data=coleman, control=lmrobdet.control(tuning.psi=bisquare(.999), 
                  bb=.02))
step.lmrobdet(m2)
step.lmrobdet(m2.ls)
step(lm(Y~., data=coleman), direction='backward') #, scale=2.074296) 

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

x1 <- rnorm(50)
x2 <- rnorm(50)
y <- rnorm(50, sd=.1) + x1 - 2*x2
das <- data.frame(x1=x1, x2=x2, y=y)
tmp <- lmrobdet(y ~ . , data=das)
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

