# step.R
# EXAMPLE 5.3  

# NOTE:  The sequence of models in Table 5.2 of the book is correct,
# but the RFPE values are wrong, and the ones computed below are correct.

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















