# biochem.R
# EXAMPLE 6.1   
# Figures 6.1, 6.2
# Table 6.1

data(biochem, package='RobStatTM')
X <- as.matrix(biochem)
colnames(X) <- c('Phosphate', 'Chloride')
plot(X, pch=19, main='Biochem Data scatterplot')
text(.95,4.4,"3",col = "black",cex = 1.2)

qqnorm(X[,1], pch=19)
qqline(X[,1], lwd=2, col='gray20')

mu <- colMeans(X)
cov.mat <- var(X)
vv <- diag(cov.mat)
rho <- cor(X)[1,2]
a <- cbind(t(mu), t(vv), rho)
X3 <- X[-3,]    #delete obs. 3 and recompute
mu2 <- colMeans(X3)
cov.mat2 <- var(X3)
vv2 <- diag(cov.mat2)
rho2 <- cor(X3)[1,2]
a2  <- cbind(t(mu2), t(vv2), rho2)
print("Means  Vars,   Correl")
print(round(rbind(a,a2),2))












