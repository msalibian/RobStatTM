# epilepsy.R
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

 ufact <- 1.1
 ctrl <-  robcbi::cubinf.control(ufact=ufact)
 epiCUBIF <- robcbi::cubinf(XX, yy, family=poisson(), null.dev = FALSE, control=ctrl)
 epiCUBIF_coefficients <- epiCUBIF$coefficients
 epiCUBIF_ste <- sqrt(diag(epiCUBIF$cov))
 epiCUBIF_deviances <- sign(yy-epiCUBIF$fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiCUBIF$fitted)+epiCUBIF$fitted))

#ML estimator
epiMV <- glm(yy~xx1+xx2+xx3+xx4,family=poisson)
epiMV_coefficients <- epiMV$coefficients
tt <- summary(epiMV)
epiMV_ste <- sqrt(diag(tt$cov.unscaled))
epiMV_deviances <- tt$deviance.resid

# NOTE 1:  The MLE values in Table 7.3 of the book are incorrect.

#RQL estimator
epiRQL <- robustbase::glmrob(yy~xx1+xx2+xx3+xx4,family=poisson)
epiRQL_coefficients <- epiRQL$coefficients
epiRQL_ste <- sqrt(diag(epiRQL$cov))
epiRQL_deviances <- residuals(epiRQL)


#MT estimator
epiMT <- robustbase::glmrob(yy~xx1+xx2+xx3+xx4, family=poisson, method="MT")
epiMT_coefficients <- epiMT$coefficients
epiMT_ste <- sqrt(diag(epiMT$cov))
epiMT_deviances <- sign(yy-epiMT$fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiMT$fitted)+epiMT$fitted))


#epiMP
# NOTE 2: There is no R code for the MP estimator, so the following was computed with MATLAB
epiMP_coefficients <- c(2.0078,    0.0707,    0.1346,   -0.4898  ,  0.0476)
epiMP_fitted <- exp(XX%*%epiMP_coefficients)
epiMP_deviances <- sign(yy-epiMP_fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiMP_fitted)+epiMP_fitted))

#boxplot Figure 7.6
par(mfrow=c(1,2), cex=0.6)
# devC <- cbind(abs(epiMV_deviances),abs(epiCUBIF_deviances),abs(epiMT_deviances),abs(epiRQL_deviances),abs(epiMP_deviances))
# boxplot(devC,names=c("ML","CUBIF","MT","QL","MP"),ylab="Absolute Deviance Residuals",main="with outliers")
# boxplot(devC,names=c("ML","CUBIF","MT","QL","MP"),ylab="Absolute Deviance Residuals",main="without outliers",outline=FALSE)
devC <- cbind(abs(epiMV_deviances),abs(epiMT_deviances),abs(epiRQL_deviances),abs(epiMP_deviances))
boxplot(devC,names=c("ML","MT","QL","MP"),ylab="Absolute Deviance Residuals",main="with outliers")
boxplot(devC,names=c("ML","MT","QL","MP"),ylab="Absolute Deviance Residuals",main="without outliers",outline=FALSE)
par(mfrow = c(1,1))

# Figure 7.7
sdev1 <- sort(abs(epiMT_deviances))
sdev2 <- sort(abs(epiMV_deviances))

plot(ppoints(59)[1:48],sdev1[1:48], type="b",pch=1,xlab="quantiles", ylab= "deviance residuals")
lines(ppoints(59)[1:48],sdev2[1:48], type="b",pch=2)
legend(x="topleft",legend=c("MT","ML"), pch=c(1,2))













