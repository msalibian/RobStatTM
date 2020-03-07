# leukemia.R
# Example 7.1
# Figure 7.4
# Table 7.1


library(RobStatTM)

# Must install robust
data(leuk.dat) 

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
ufact <- 1.1
ctrl <-  robcbi::cubinf.control(ufact=ufact)
yy <- leuk.dat$y
XX <- cbind(rep(1, length(yy)), leuk.dat$wbc, leuk.dat$ag)
leukCUBIF <- robcbi::cubinf(XX, yy, family=binomial(), null.dev = FALSE, control=ctrl)


#weighted ML fit
leukWML <- logregWML(Xleuk, yleuk, intercept=1)
















