# ALL SCRIPTS

# EXAMPLE 4.1 shocks.R  shocks.txt
# Figure 4.1 and 4.3

library(quantreg)  # For L1 fit
library(RobStatTM)
data(shock, package='RobStatTM')
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

#LS fit
shockls <- lm(time ~ n.shocks, data = shock)

#LS fit without outliers
shockls124 <- lm(time ~ n.shocks, data = shock, subset = -c(1, 2, 4))

#--------------------
#Figure 4.1
plot(time ~ n.shocks, data=shock, xlab="number of shocks", ylab="average time", pch=19, cex=1.3)
abline(shockls124, lwd=2, col='gray30')
abline(shockls, lwd=2, col='tomato')
text(shock[c(1,2,4),1],shock[c(1,2,4),2]-.4, labels=c("1","2","4"), cex=1.2)
text(1,10.5,"LS", col='tomato', cex=1.3)
text(1,7.2,"LS -", col='gray30', cex=1.3)
#---------------------------------

#L1 fit
shockl1 <- rq(time~n.shocks , data=shock)

# M fit
shockrob <- lmrobM(time ~ n.shocks, data = shock,control=cont)


#--------------------------
#Figure 4.3
plot(time ~ n.shocks, data=shock, xlab="number of shocks", ylab="average time", pch=19, cex=1.3)
abline(shockls124, lwd=2, col='gray30')
abline(shockls, lwd=2, col='tomato')
abline(shockrob, lwd=2, col='green')
abline(shockl1, lwd=2, col='blue')
text(shock[c(1,2,4),1],shock[c(1,2,4),2]-.4, labels=c("1","2","4"), cex=1.3)
text(1,10.5,"LS", cex=1.3)
text(1,7.2,"LS -", cex=1.3)
text(1,8.3,"L1", cex=1.3)
text(1,7.7,"M", cex=1.3)


# EXAMPLE 4.2  oats.R  oats.txt
# Figure 4.2 and 4.4

library(RobStatTM)
data(oats, package='RobStatTM')

## LS regression
oats1LS <- lm(response1 ~ variety+block, data=oats)
oats1LS_var <- lm(response1 ~ block, data=oats)
oats1LS_block <- lm(response1 ~ variety, data=oats)
oats2LS <- lm(response2 ~ variety+block, data=oats)
oats2LS_var <- lm(response2 ~ block, data=oats)
oats2LS_block <- lm(response2 ~ variety, data=oats)
sigma1LS <- summary(oats1LS)$sigma
sigma2LS <- summary(oats2LS)$sigma

cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

## Classical ANOVA  tests
anov1_var  <- anova(oats1LS, oats1LS_var)
anov1_block <- anova(oats1LS, oats1LS_block)
anov2_var <- anova(oats2LS, oats2LS_var)
anov2_block <- anova(oats2LS, oats2LS_block)

## M regressions
oats1M <- lmrobM(response1 ~ variety+block, control=cont, data=oats)
oats1M_var <- lmrobM(response1 ~ block, control=cont, data=oats)
oats1M_block <- lmrobM(response1 ~ variety, control=cont, data=oats)
oats2M <- lmrobM(response2 ~ variety+block, control=cont, data=oats)
oats2M_var <- lmrobM(response2 ~ block, control=cont, data=oats)
oats2M_block <- lmrobM(response2 ~ variety, control=cont, data=oats)
sM2 <- oats2M$scale

## Robust ANOVA  tests
anov1M_var <- rob.linear.test(oats1M, oats1M_var)
anov1M_block <- rob.linear.test(oats1M, oats1M_block)
anov2M_var <- rob.linear.test(oats2M, oats2M_var)
anov2M_block <- rob.linear.test(oats2M, oats2M_block)

plot(oats2LS, which=2)
abline(h=0, lty=2)


tmp <- qqnorm(resid(oats2M)/sM2, ylab="Standardized residuals", pch=19, col='gray30')
qqline(resid(oats2M)/sM2)
abline(h=c(-2.5, 0, 2.5), lty=2)
w <- c(24,36,1,35,20)
text(tmp$x[w] + .1, tmp$y[w] + .1, w)



# EXAMPLE 5.1  mineral.R  mineral.txt
# Figures 5.1 - 5.7

library(quantreg)
library(RobStatTM)

data(mineral, package='RobStatTM')
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

#LS fit
mineralls <- lm(zinc ~  copper, data=mineral)

#L1 fit
minerall1 <- rq(zinc ~  copper, data=mineral)

#LS without outlier fit
minerallssin <- lm(zinc ~ copper, data = mineral[-15,  ])

#Fig 5.1.
plot(zinc ~ copper, data=mineral, pch=19, cex=1.3)
abline(mineralls, lwd=2, col='red')
abline(minerall1, lwd=2, col='blue')
abline(minerallssin, lwd=2, col='green4')
text(c(600,600,600), c(29,55,82), c("LS-15","L1", "LS"), cex=1.3)
text(mineral[15,1], mineral[15,2]-6, "15", cex=1.3)



#Figure 5.2
plot(mineralls, which=2, add.smooth=FALSE, pch=19, id.n=2, cex.id = 1.2)
abline(h=c(2.5, 0, -2.5), lty=2, lwd=2)

#Figure 5.3
sigmaLS <- summary( mineralls)$sigma
plot(mineralls, which=1, add.smooth=FALSE, pch=19, id.n=2, cex.id = 1.2)
abline(h=c(2.5, 0, -2.5)*sigmaLS, lty=2, lwd=2)

#--------------------------------------------


#MM fit
mineralMM <- lmrobdetMM(zinc ~ copper, data = mineral, control=cont)

#Fig 5.4
plot(zinc ~ copper, data=mineral, pch=19, cex=1.3)
abline(minerallssin, lwd=2, col='green4')
abline(mineralMM, lwd=2, col='magenta')
text(c(600,600), c(45,29), c("ROB","LS-15"), cex=1.2)
text(mineral[15,1], mineral[15,2]-6, "15", cex=1.2)
#------------------------------------------------------------


#------------------------------------------------------------
#Fig 5.5
plot(mineralMM, which=4, add.smooth=FALSE, pch=19, cex.id=1.3)

#-----------------------------------------------

#Fig 5.6
plot(mineralMM, which=2, pch=19, id.n=3, cex.id=1.2)
abline(h=c(-2.5, 0, 2.5) * mineralMM$scale, lty=2, lwd=2)

#----------------------------------------------------------

#Fig 5.7
plot(sort(abs(resid(mineralls)))[-53], sort(abs(resid(mineralMM)))[-53], xlab="Least Squares residuals",
     ylab="Robust residuals", pch=19, cex=1.3)
abline(0, 1, lwd=2, col='red')




# EXAMPLE 5.2  wood.R  wood from robustbase
# Figures 5.8 - 5.12

library(RobStatTM)
data(wood, package = "robustbase")
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

#MM fit
woodMM <- lmrobdetMM(y ~ ., data=wood, control=cont)

#LS fit
woodLS <- lm(y ~ ., data=wood)

#-------------------------------------------------------
#Fig 5.8
# Nothing happens in this figure
# text(woodLS$fitted[28]+.03,woodLSresst[28],"28")
sigmaLS <- summary(woodLS)$sigma
plot(woodLS, which=1, add.smooth=FALSE, pch=19, id.n=2, cex.id = 1.2)
abline(h=c(-2.5, 0, 2.5) * sigmaLS, lty=2, lwd=2)

#----------------------------
#Figure 5.9
plot(woodLS, which=2, pch=19)


#-------------------------------------------
# Fig. 5.10
plot(woodMM, which=4, add.smooth=FALSE, pch=19, cex.id=1.3, id.n=4)


#--------------------------
#Figure 5.11
plot(woodMM, which=2, add.smooth=FALSE, pch=19, cex.id=1.3, id.n=4)
abline(h=c(-2.5, 0, 2.5) * woodMM$scale, lty=2)
#----------------------------

#Figure 5.12
wq <- which( abs(resid(woodMM)) > 2.5 * woodMM$scale )
plot(sort(abs(resid(woodLS)[-wq] ) ), sort(abs(resid(woodMM)[-wq])), pch=19, xlim=c(0, .05), ylim=c(0, .05))
abline(0,1)


# EXAMPLE 5.3  step.R

library(RobStatTM)

cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")

set.seed(300)
X <- matrix(rnorm(50*6), 50, 6)
beta <- c(1,1,1,0,0,0)
y <- as.vector(X %*% beta) + 1 + rnorm(50)
y[1:6] <- seq(30, 55, 5)

for (i in 1:6) X[i,] <- c(X[i,1:3],i/2,i/2,i/2)
Z <- cbind(y,X)
Z <- as.data.frame(Z)
obj <- lmrobdetMM(y ~ ., data=Z, control=cont)
out <- step.lmrobdetMM(obj)

obj2 <- lm(y ~ ., data=Z)
out2 <- step(obj2)


# EXAMPLE 5.4  algae.R  algae,txt
# Figures 5.14-5.15

library(RobStatTM)
data(algae, package='RobStatTM')

#Robust fit
cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")
algaerob <- lmrobdetMM(V12 ~ ., data=algae, control=cont)

#LS fit
algaels <- lm(V12 ~ ., data=algae)

#LS fit without outliers
algaelsd <- lm(V12 ~ ., data=algae, subset= -c(36, 77))

#-----------------------------------------------------
#Fig 5.14
plot(algaels, which=2, pch=19)
abline(h=c(-2.5, 0, 2.5), lty=2)

#-------------------------------------------------------
#Fig 5.15
plot(algaerob, which=2, id.n=2, pch=19)
abline(h=c(-2.5, 0, 2.5)*algaerob$scale, lty=2)

#-----------------------------------------------

# EXAMPLE 5.5  ExactFit.R

library(RobStatTM)
library(pyinit)
set.seed(1003)
n <- 100
m <- 50
rr <- rnorm(m)
x1 <- sort(rnorm(n))
x2 <- sort(rr)*2
sig <- 0.1
y1 <- x1 + sig*rnorm(n)   # "good"  data
y2 <- -x2 + sig*rnorm(m)  # outliers
x <- c(x1,x2)
y <- c(y1,y2)
out1 <- lm(y~x) # LSE
out2 <- lmrobdetMM(y~x)  #MM

plot(y ~ x, pch=19, col='gray30')
abline(out1, lwd=3, col='blue3')
abline(out2, lwd=3, col='red3')
text(c(-3.5,3),c(0,2),c("LS","MM"), cex=1.3, col=c('blue3', 'red3'))
#--------------------------------------------------

# EXAMPLE 6.1   biochem.R  biochemdata.txt
# Figures 6.1, 6.2 and Table 6.1
data(biochem, package='RobStatTM')
X <- as.matrix(biochem)
colnames(X) <- c('Phosphate', 'Chloride')
plot(X, pch=19, main='Biochem Data scatterplot')

qqnorm(X[,1], pch=19)
qqline(X[,1], lwd=2, col='gray20')

mu <- colMeans(X)
cov.mat <- var(X)
vv <- diag(cov.mat)
rho <- cov.mat[1,2]
a <- cbind(t(mu), t(vv), rho)
X3 <- X[-3,]    #delete obs. 3 and recompute
mu2 <- colMeans(X3)
cov.mat2 <- var(X3)
vv2 <- diag(cov.mat2)
rho2 <- cov.mat2[1,2]
a2  <- cbind(t(mu2), t(vv2), rho2)
print("Means  Vars,   Correl")
print(rbind(a,a2))
#----------------------------------------------

# EXAMPLE 6.2  wine.R  wine.txt
# Figure 6.3

library(RobStatTM)
data(wine, package='RobStatTM')
X <- as.matrix(wine)
xbar <- colMeans(X)
C <- cov(X)
disC <- mahalanobis(X, xbar, C)
resu <- MultiRobu(X, type="MM")
mu <- resu$mu
V <- resu$V
disM <- mahalanobis(X, mu, V)

resu <- MultiRobu(X, type="Rocke");
mu <- resu$mu
V <- resu$V
disR <- mahalanobis(X,mu,V)

#Figure 6.3
par(mfrow=c(2,2))
plot(disC, xlab="index", ylab="Distances", main = "Classical", cex.main=0.9, pch=19)
plot(qchisq(ppoints(59),13),sort(disC),xlab="chi squared quantiles", ylab="Sorted distances", main ="Classical", cex.main=0.9, pch=19)
lines(sort(disC),sort(disC))
plot(disM, xlab="index", ylab="Distances", main = "Robust", cex.main=0.9, pch=19)
plot(qchisq(ppoints(59),13),sort(disM),xlab="chi squared quantiles", ylab="Sorted distances", main="Robust", cex.main=0.9, pch=19)
lines(sort(disM),sort(disM))
#----------------------------------------------------------

#EXAMPLE 6.3    vehicle.R   vehi4.txt7
# Figure 6.7

library(RobStatTM)
library("rrcov")
data(vehicle, package='RobStatTM')
X <- as.matrix(vehicle)
n <- dim(X)[1]
p <- dim(X)[2]
xbar <- colMeans(X)
C <- cov(X)
disC1 <- mahalanobis(X,xbar,C);  disC=sort(disC1)  #Classical estimator
resu <- MultiRobu(X,"Rocke") #Rocke estimator
muR <- resu$mu
VR <- resu$V
disR1 <- mahalanobis(X, muR, VR)
disR <- sort(disR1)
resu <- CovMcd(X)
muM <- slot(resu, 'center')
VM <- slot(resu, 'cov')
disM1 <- slot(resu, 'mah')
disM <- sort(disM1)
resu <- CovSest(X, method= "bisquare")  #S-estimator
muS <- slot(resu, 'center')
VS <- slot(resu, 'cov')
disS1 <- mahalanobis(X, muS, VS)
disS <-sort(disS1)

# Figure 6.7
par(mfrow=c(1,3))
qua <- qchisq(ppoints(n), p)
plot(qua,disC, xlab="Chi squared quantiles", ylab="Sorted distances", main="Classical", cex.main=0.9, pch=19)
abline(0,1)
plot(qua, disS, xlab="Chi squared quantiles", ylab="Sorted distances", main="S-Bisquare", cex.main=0.9, pch=19)
abline(0,1)
plot(qua, disR, xlab="Chi squared quantiles", ylab="Sorted distances", main="Rocke", cex.main=0.9, pch=19)
abline(0,1)
#---------------------------------------------

# EXAMPLE 6.4   bus.R  busdata.txt
# PCA  Table 6.6 and Figure 6.10

library(RobStatTM)
ClassPCA <- function(X, ncomp) {
  mu <- colMeans(X)
  Xcen <- scale(X, center=mu, scale=FALSE)
  sx <- svd(Xcen, nu=0)
  VCC <- sx$v
  Q <- VCC[,1:ncomp]
  lam <- sx$d^2
  fit <- scale(Xcen%*%Q%*%t(Q), center=-mu, scale=FALSE)
  return( list(fit=fit, evec=VCC, lam=lam) )
}
data(bus, package='RobStatTM')
X0 <- as.matrix(bus)
X1 <- X0[,-9]
ss <- apply(X1, 2, mad)
mu <- apply(X1, 2, median)
X <- scale(X1, center=mu, scale=ss)
n <- dim(X)[1]
p <- dim(X)[2]

#Classical PCA
q <- 3  #compute three components
resC <- ClassPCA(X, q)
lamC <- resC$lam
lamC <- lamC/sum(lamC)
prC <- cumsum(lamC)
nonC <- 1-prC  #proportion of unexaplained variance
fitC <- resC$fit
resiC <- X-fitC
dC <- rowSums(resiC^2)

#Robust PCA
rr <- SMPCA(X, q, 0.99)
propex <- rr$propex
fitM <- rr$fit
resiM <- X-fitM
dM <- rowSums(resiM^2)
alfa <- seq(from=0.1, to=0.9, by=0.1)
qC <- quantile(dC,alfa)
qM <- quantile(dM,alfa)
print(rbind(qC,qM))

dCsor <- sort(dC)
dMsor <- sort(dM)
par(mfrow=c(1,1))
plot(dCsor, dMsor, xlab="Classic", ylab="Robust", log="xy")
abline(0,1)
#---------------------------------------------


# Example 6.5 and 6.6

library(GSE)
library(RobStatTM)
data(wine, package='RobStatTM')
X <- as.matrix(wine)
DM <- dim(X)
n <- DM[1]
p <- DM[2]


#omitted data Figure 6.11
#GSE
set.seed(100)
RR <- matrix(runif(n*p)<.2,n,p)
X2 <- X
for (i in 1:n) {
  for(j in 1:p) {
    if(RR[i,j]) X2[i,j] <- NA
  }
}

qq <- qchisq(.999, p)

out2 <- GSE(X2)
md2 <- getDistAdj(out2)
smd2 <- sort(md2)
v2 <- (md2 > qq)
z2 <- 1:n
z2 <- z2[v2]
#number of outliers
no2 <- length(z2)

#EM estimator
out3 <- CovEM(X2)
md3 <- getDistAdj(out3)
smd3 <- sort(md3)
v3 <- ( md3 > qq )
z3 <- 1:n
z3 <- z3[v3]
#number of outliers
no3 <- length(z3)

#Figure 6.11
abc <- qchisq(ppoints(n), p)
par(mfrow=c(2,2))
plot(md3, xlab="Index", ylab="Adjusted distances", pch=19)
abline(h=qq)
plot(abc, abc, xlab="Chi-square quantiles", ylab="Adjusted distance quantiles", type='n', pch=19)
points(abc,smd3)
abline(0,1)
plot(md2, xlab="Index", ylab="Adjusted distances", pch=19)
abline(h=qq)
plot(abc, smd2, xlab="Chi square quantile", ylab="Adjusted distance quantiles", pch=19)
abline(0,1)


#----------------------------------------------------------
#Analysis with independent contamination Figure 6.12

#MM
set.seed(100)
out <- MultiRobu(X, type="MM")
md <- out$dist
smd <- sort(md)
v <- ( md> qq )
z <- 1:n
z <- z[v]
# number of outliers
no <- length(z)

#TSGS
out4 <- TSGS(X, method="bisquare", init="emve", filter="UBF-DDC")
mu4 <- getLocation(out4)
Sigma4 <- getScatter(out4)
winef4 <- getFiltDat(out4)
md4 <- mahalanobis(X, mu4, Sigma4)
smd4 <- sort(md4)
v4 <- ( md4 > qq )
z4 <- 1:n
z4 <- z4[v4]
# number of outliers
no4 <- length(z4)


#Figure 6.12
abc <- qchisq(ppoints(n),p)
par(mfrow=c(2,2))
plot(md, xlab="Index", ylab="Adjusted distances", pch=19)
abline(h = qq) #lines(c(0,60),c(qq,qq))
plot(abc, smd, xlab="Chi square quantiles", ylab="Adjusted distance quantiles", pch=19)
abline(0,1)
plot(md4, xlab="Index", ylab="Adjusted distances", pch=19)
abline(h=qq)
plot(abc, smd4, xlab="Chi square quantile", ylab="Adjusted distance quantiles", pch=19)
abline(0,1)

#---------------------------------------------------


# Example 6.7
### AUTISM
library(RobStatTM)
library(robustvarComp)
library(nlme)
data(autism, package='WWGbook')
autism <- autism[complete.cases(autism),]
completi <- table(autism$childid)==5
completi <- names(completi[completi])
indici <- as.vector(unlist(sapply(completi, function(x) which(autism$childid==x))))
ind <- rep(FALSE, nrow(autism))
ind[indici] <- TRUE
autism <- subset(autism, subset=ind) ## complete cases 41
attach(autism)
sicdegp.f <- factor(sicdegp)
age.f <- factor(age)
age.2 <- age - 2
sicdegp2 <- sicdegp
sicdegp2[sicdegp == 3] <- 0
sicdegp2[sicdegp == 2] <- 2
sicdegp2[sicdegp == 1] <- 1
sicdegp2.f <- factor(sicdegp2)
autism.updated <- subset(data.frame(autism, sicdegp2.f, age.2), !is.na(vsae))
autism.grouped <- groupedData(vsae ~ age.2 | childid, data=autism.updated, order.groups = FALSE)
p <- 5
n <- 41
z1 <- rep(1, p)
z2 <- c(0, 1, 3, 7, 11)
z3 <- z2^2
K <- list()
K[[1]] <- tcrossprod(z1,z1)
K[[2]] <- tcrossprod(z2,z2)
K[[3]] <- tcrossprod(z3,z3)
K[[4]] <- tcrossprod(z1,z2) + tcrossprod(z2,z1)
K[[5]] <- tcrossprod(z1,z3) + tcrossprod(z3,z1)
K[[6]] <- tcrossprod(z3,z2) + tcrossprod(z2,z3)
names(K) <- c("Int", "age", "age2", "Int:age", "Int:age2", "age:age2")

groups <- cbind(rep(1:p, each=n), rep((1:n), p))

## Composite Tau
AutismCompositeTau <- varComprob(vsae ~ age.2 + I(age.2^2)
                                 + sicdegp2.f + age.2:sicdegp2.f + I(age.2^2):sicdegp2.f,
                                 groups = groups, data = autism.grouped, varcov = K,
                                 control=varComprob.control(lower=c(0.01,0.01,0.01,-Inf,-Inf,-Inf)))

summary(AutismCompositeTau)

## Classic S
AutismS <- varComprob(vsae ~ age.2 + I(age.2^2)
                      + sicdegp2.f + age.2:sicdegp2.f + I(age.2^2):sicdegp2.f,
                      groups = groups, data = autism.grouped, varcov = K,
                      control=varComprob.control(method="S", psi="rocke", cov.init="covOGK", lower=c(0.01,0.01,0.01,-Inf,-Inf,-Inf)))
summary(AutismS)

# Example 7.1
# Leukemia

#loading packages
# library(robust)
# library(robustbase)
library(RobStatTM)
# library(MASS)
data(leuk.dat, package='robust')

Xleuk <-as.matrix( leuk.dat[, 1:2] )
yleuk <- leuk.dat$y
# weighted M fit
leukWBY <- logregWBY(Xleuk, yleuk, intercept=1)
pr1 <- as.vector( leukWBY$fitted.values )

# ML fit
leukML <- glm(formula = y ~ wbc + ag, family = binomial, data = leuk.dat)

#--------------------------------
#Figure 7.4

dev1 <- abs(leukWBY$residual.deviances)
dev2 <- abs(resid(leukML, type='deviance'))
n <- length(dev1)
ord1 <- order(dev1)
sdev1 <- sort(dev1) # dev1[ord1]
sdev2 <- sort(dev2) #[ord2] # typo in leukemia.R

plot(ppoints(n), sdev1, type="b",pch=1,xlab="quantiles", ylab= " deviance residuals")
lines(ppoints(n), sdev2, type="b",pch=2)
xuu <- ppoints(n)[n]
text(xuu - .03, max(sdev1) + .1, ord1[n])
text(xuu, max(sdev2) + .3, ord1[n])
legend(x="topleft",legend=c("weighted M","maximum likelihood"), pch=c(1,2))


#other estimates
#M fit
leukBY <- logregBY(Xleuk, yleuk, intercept=1)

#cubif fit
library(robust)
leukCUBIF <- glmRob(formula = y ~ wbc + ag, family = binomial, data = leuk.dat, fit.method = "cubif")

#weighted ML fit
leukWML <- logregWML(Xleuk, yleuk, intercept=1)


# Example 7.2
# Figure 7.5
# loading packages
library(RobStatTM)
# library(robustbase)
# library(robust)
# skin=read.table("skin.txt", header=T)

Xskin <- as.matrix( skin[, 1:2] )
yskin <- skin$vasoconst
#weighted M fit
skinWBY <- logregWBY(Xskin, yskin, intercept=1)

# ML fit
skinML <- glm(formula=vasoconst~logVOL+logRATE, family = binomial, data = skin)

#Figure 7.5
dev1 <- abs(skinWBY$residual.deviances)
dev2 <- abs(resid(skinML, type='deviance'))
sdev1 <- sort(dev1)
sdev2 <- sort(dev2)
tt <- c(18,4)
uu <- order(tt)
xuu <- ppoints(39)[38:39]
yuu1 <- sdev1[38:39]
yuu2 <- sdev2[38:39]
tx <- c("18","4")

plot(ppoints(39), sdev1, type="b", pch=1, xlab="quantiles", ylab= "deviance residuals")
lines(ppoints(39), sdev2, type="b", pch=2)
text(xuu, yuu1 + .1, tt)
text(xuu, yuu2 + .1, tt)
legend(x="topleft", legend=c("weighted M","maximum likelihood"), pch=c(1,2))


# other estimates
# M fit
skinBY <- logregBY(Xskin, yskin, intercept=1)

#cubif fit
library(robust)
skinCUBIF <- glmRob(formula =vasoconst~logVOL+logRATE, family = binomial, data = skin, fit.method = "cubif")

# weighted ML fit
skinWML <- logregWML(Xskin, yskin, intercept=1)

# Example 7.3
# Breslow data

library(RobStatTM)
data(breslow.dat, package='robust')


#CUBIF Estimator
yy <- breslow.dat[, 10]
xx1 <- breslow.dat[, 11]
xx2 <- breslow.dat[, 12]
xx3 <- breslow.dat[, 8]=="progabide"
xx4 <- xx2*xx3

XX <- cbind(rep(1,59), xx1, xx2, xx3, xx4)
colnames(XX) <- c("intercept","Age10","Base4","Progabide","interac.Base4-Progabide")

# # Does not run
# ufact <- 1.1
# ctrl <-  cubinf.control(ufact=ufact)
# epiCUBIF <- cubinf(XX, yy, family=poisson(), null.dev = FALSE, control=ctrl)
# epiCUBIF_coefficients <- epiCUBIF$coefficients
# epiCUBIF_ste <- sqrt(diag(epiCUBIF$cov))
# epiCUBIF_deviances <- sign(yy-epiCUBIF$fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiCUBIF$fitted)+epiCUBIF$fitted))

#ML estimator
epiMV <- glm(yy~xx1+xx2+xx3+xx4,family=poisson)
epiMV_coefficients <- epiMV$coefficients
tt <- summary(epiMV)
epiMV_ste <- sqrt(diag(tt$cov.unscaled))
epiMV_deviances <- tt$deviance.resid

#RQL estimator
library(robustbase)
epiRQL <- glmrob(yy~xx1+xx2+xx3+xx4,family=poisson)
epiRQL_coefficients <- epiRQL$coefficients
epiRQL_ste <- sqrt(diag(epiRQL$cov))
epiRQL_deviances <- residuals(epiRQL)


#MT estimator
epiMT <- glmrob(yy~	xx1+xx2+xx3+xx4, family=poisson, method="MT")
epiMT_coefficients <- epiMT$coefficients
epiMT_ste <- sqrt(diag(epiMT$cov))
epiMT_deviances <- sign(yy-epiMT$fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiMT$fitted)+epiMT$fitted))


#epiMP
epiMP_coefficients <- c(2.0078,    0.0707,    0.1346,   -0.4898  ,  0.0476)
epiMP_fitted <- exp(XX%*%epiMP_coefficients)
epiMP_deviances <- sign(yy-epiMP_fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiMP_fitted)+epiMP_fitted))

#boxplot Figure 7.6
par(mfrow=c(1,2), cex=0.6)
devC <- cbind(abs(epiMV_deviances),abs(epiCUBIF_deviances),abs(epiMT_deviances),abs(epiRQL_deviances),abs(epiMP_deviances))
boxplot(devC,names=c("ML","CUBIF","MT","QL","MP"),ylab="Absolute Deviance Residuals",main="with outliers")
boxplot(devC,names=c("ML","CUBIF","MT","QL","MP"),ylab="Absolute Deviance Residuals",main="without outliers",outline=FALSE)

# Figure 7.7

sdev1 <- sort(abs(epiMT_deviances))
sdev2 <- sort(abs(epiMV_deviances))

plot(ppoints(59)[1:48],sdev1[1:48], type="b",pch=1,xlab="quantiles", ylab= "deviance residuals")
lines(ppoints(59)[1:48],sdev2[1:48], type="b",pch=2)
legend(x="topleft",legend=c("MT","ML"), pch=c(1,2))













