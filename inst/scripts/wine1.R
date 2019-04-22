# wine1.R
# Example 6.5 and 6.6

# Must install GSE

library(RobStatTM)
data(wine)
X <- as.matrix(wine)
DM <- dim(X)
n <- DM[1]
p <- DM[2]


#omitted data Figure 6.11
#GSE
set.seed(2400)
RR <- matrix(runif(n*p)<.2,n,p)
X2 <- X
for (i in 1:n) {
  for(j in 1:p) {
    if(RR[i,j]) X2[i,j] <- NA
  }
}

qq <- qchisq(.999, p)

out2 <- GSE::GSE(X2)
md2 <- GSE::getDistAdj(out2)
smd2 <- sort(md2)
v2 <- (md2 > qq)
z2 <- 1:n
z2 <- z2[v2]
#number of outliers
no2 <- length(z2)

#EM estimator
out3 <- GSE::CovEM(X2)
md3 <- GSE::getDistAdj(out3)
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
par(mfrow=c(1,1))

# NOTE:  The difference between the plots made by the script above and
# the plots in Figure 6.11 of the book are due to the use of a different
# seed for the generation of random numbers.


#----------------------------------------------------------
#Analysis with independent contamination Figure 6.12

#MM
set.seed(2400)
out <- covRobMM(X)
md <- out$dist
smd <- sort(md)
v <- ( md> qq )
z <- 1:n
z <- z[v]
# number of outliers
no <- length(z)

#TSGS
out4 <- GSE::TSGS(X, method="bisquare", init="emve", filter="UBF-DDC")
mu4 <- GSE::getLocation(out4)
Sigma4 <- GSE::getScatter(out4)
winef4 <- GSE::getFiltDat(out4)
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
par(mfrow=c(1,1))
