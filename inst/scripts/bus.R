# bus.R
# EXAMPLE 6.4    
# Figure 6.10
# Table 6.6 

library(RobStatTM)
data(bus)
X0 <- as.matrix(bus)
X1 <- X0[,-9]
ss <- apply(X1, 2, mad)
mu <- apply(X1, 2, median)
X <- scale(X1, center=mu, scale=ss)
n <- dim(X)[1]
p <- dim(X)[2]

#Classical PCA
q <- 3  #compute three components
resC <- prcomp(X) 
prC <- as.vector(summary(resC)$importance['Cumulative Proportion', ] )
nonC <- 1 - prC  #proportion of unexaplained variance
Xcent <- scale(X, center=resC$center, scale=FALSE)
fitC <- scale(Xcent %*% resC$rotation[, 1:q] %*% t(resC$rotation[, 1:q]), center=-resC$center, scale=FALSE) #resC$fit
resiC <- X - fitC
dC <- rowSums(resiC^2)

#Robust PCA
rr <- pcaRobS(X, q, 0.99)
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
















