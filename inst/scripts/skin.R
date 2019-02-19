# skin.R
# Example 7.2
# Figure 7.5

library(RobStatTM)
data(skin)

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
skinCUBIF <- robust::glmRob(formula =vasoconst~logVOL+logRATE, family = binomial, data = skin, fit.method = "cubif")

# weighted ML fit
skinWML <- logregWML(Xskin, yskin, intercept=1)















