library(RobustStatistics)
# source('RFPE-new.R')


set.seed(123)
n <- 200
x1 <- rnorm(n)
x2 <- runif(n)
x3 <- rexp(n, rate=.2)
x4 <- rnorm(n)
x5 <- rnorm(n, sd=2)
x6 <- rbinom(n, size=3, prob=.3)
x7 <- rt(n, df=2)
x8 <- rpois(n, lambda=5)
# f1 <- as.factor( LETTERS[ rbinom(n, size=2, prob=.3) + 1] )
# f2 <- factor( c('Low', 'Med', 'High')[ rbinom(n, size=2, prob=.3) + 1], levels=c('Low', 'Med', 'High'))
y <- rnorm(n, sd=.5) + 2*x1  # - 3*as.numeric(f2)
d <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, x6=x6, x7=x7, x8=x8, y = y)


a <- lmrob(y ~ ., data=d, control=lmrob.control(k.max=200, refine.tol=1e-3)) # , init='M-S')

a <- lmrobdet(y ~ ., data=d, control=lmrobdet.control(refine.tol=1e-3))
# pyinit bombs!
a <- lmrobdet(y ~ 1, data=d, control=lmrobdet.control(refine.tol=1e-3))
b <- lmrob(y ~ 1, data=d, control=lmrob.control(refine.tol=1e-3))


u <- step.lmrobdet(a, whole.path=TRUE, trace=FALSE)
# is not using MM estimates!
u2 <- step.lmrobdet(a)

lmrobdet.RFPE(lmrobdet(y~x1+x2+x7+x8, data=d, control=lmrobdet.control(refine.tol=1e-3))$MM, scale=a$MM$scale)
lmrobdet.RFPE(lmrobdet(y~x1+x8, data=d, control=lmrobdet.control(refine.tol=1e-3))$MM, scale=a$MM$scale)

fo <- as.formula(names(u[4]))
lmrobdet.RFPE(lmrobdet(formula=fo, data=d, control=lmrobdet.control(refine.tol=1e-3))$MM, scale=a$MM$scale)


# keep x4 in all models
my.scope <- list(lower = . ~ x4, upper = . ~ .)
step.lmrob(a, scope = my.scope, whole.path=TRUE, trace=FALSE)


data(Boston, package='MASS')
Bos2 <- Boston[, c('medv', 'crim', 'nox', 'age', 'dis')]
a <- lmrob(medv ~ crim + nox + age + dis, data=Boston, control=lmrob.control(k.max=200, refine.tol=1e-4))

step.lmrob(a)
step.lmrob(a, whole.path=TRUE)

data(stack.dat, package='robust')
a <- lmrob(Loss ~ ., data = stack.dat)

step.lmrob(a)
step.lmrob(a, whole.path=TRUE, trace=FALSE)


## Keep Water.Temp in the model ##
my.scope <- list(lower = . ~ Acid.Conc., upper = . ~ .)
step.lmrob(a, scope = my.scope)
step.lmrob(a, scope = my.scope, whole.path=TRUE, trace=FALSE)


