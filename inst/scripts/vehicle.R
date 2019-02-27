# vehicle.R
# EXAMPLE 6.3
# Figure 6.7

library(RobStatTM)
data(vehicle)
X <- as.matrix(vehicle)
n <- dim(X)[1]
p <- dim(X)[2]
xbar <- colMeans(X)
C <- cov(X)
disC1 <- mahalanobis(X,xbar,C);  disC=sort(disC1)  #Classical estimator
resu <- covRobRocke(X) #Rocke estimator
muR <- resu$mu
VR <- resu$V
disR1 <- mahalanobis(X, muR, VR)
disR <- sort(disR1)
resu <- rrcov::CovMcd(X)
muM <- slot(resu, 'center')
VM <- slot(resu, 'cov')
disM1 <- slot(resu, 'mah')
disM <- sort(disM1)
resu <- rrcov::CovSest(X, method= "bisquare")  #S-estimator
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
par(mfrow=c(1,1))








