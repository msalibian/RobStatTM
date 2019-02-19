# ExactFit.R
# EXAMPLE 5.5  

library(RobStatTM)
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
















