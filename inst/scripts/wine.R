# wine.R
# EXAMPLE 6.2
# Figure 6.3

library(RobStatTM)
data(wine)
X <- as.matrix(wine)
xbar <- colMeans(X)
C <- cov(X)
disC <- mahalanobis(X, xbar, C)
resu <- covRobMM(X)
mu <- resu$mu
V <- resu$V
disM <- mahalanobis(X, mu, V)

resu <- covRobRocke(X);
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
par(mfrow=c(1,1))












