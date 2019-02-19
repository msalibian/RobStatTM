

library(robustbase)
library(robust)
data(breslow.dat)
 

#CUBIF Estimator
yy<-breslow.dat[,10]
xx1<-breslow.dat[,11]
xx2<-breslow.dat[,12]
xx3<-breslow.dat[,8]=="progabide"
xx4<-xx2*xx3

XX<-cbind(rep(1,59),xx1,xx2,xx3,xx4)
colnames(XX)=c("intercept","Age10","Base4","Progabide","interac.Base4-Progabide")
ufact = 1.1
ctrl = cubinf.control(ufact=ufact)
epiCUBIF=cubinf(XX,yy,family=poisson(),null.dev = FALSE, control=ctrl)
epiCUBIF_coefficients =epiCUBIF$coefficients
epiCUBIF_ste=sqrt(diag(epiCUBIF$cov))
epiCUBIF_deviances=sign(yy-epiCUBIF$fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiCUBIF$fitted)+epiCUBIF$fitted))

#ML estimator
epiMV<-glm(yy~xx1+xx2+xx3+xx4,family=poisson)
epiMV_coefficients=epiMV$coefficients
tt=summary(epiMV)
epiMV_ste=sqrt(diag(tt$cov.unscaled))
epiMV_deviances=tt$deviance.resid

#RQL estimator
epiRQL=glmrob(yy~xx1+xx2+xx3+xx4,family=poisson) 
epiRQL_coefficients=epiRQL$coefficients
epiRQL_ste=sqrt(diag(epiRQL$cov))
epiRQL_deviances=residuals(epiRQL)


#MT estimator
epiMT<-glmrob(yy~	xx1+xx2+xx3+xx4,family=poisson,method="MT") 
epiMT_coefficients=epiMT$coefficients
epiMT_ste=sqrt(diag(epiMT$cov))
epiMT_deviances= sign(yy-epiMT$fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiMT$fitted)+epiMT$fitted))


#epiMP
epiMP_coefficients= c(2.0078,    0.0707,    0.1346,   -0.4898  ,  0.0476)
epiMP_fitted=exp(XX%*%epiMP_coefficients)
epiMP_deviances= sign(yy-epiMP_fitted)*sqrt(2*(yy*log( pmax(yy,1))-yy-yy*log(epiMP_fitted)+epiMP_fitted))

#boxplot Figure 7.6
pdf(file='Epilepsy abs dev boxplot.pdf', bg='transparent')
par(mfrow=c(1,2),cex=0.6)
devC<-cbind(abs(epiMV_deviances),abs(epiCUBIF_deviances),abs(epiMT_deviances),abs(epiRQL_deviances),abs(epiMP_deviances))
boxplot(devC,names=c("ML","CUBIF","MT","QL","MP"),ylab="Absolute Deviance Residuals",main="with outliers")
boxplot(devC,names=c("ML","CUBIF","MT","QL","MP"),ylab="Absolute Deviance Residuals",main="without outliers",outline=FALSE)
dev.off()

#Quantiles plot Figure 7.7
#pdf(file='Epilepsy abs dev quantiles plot.pdf', bg='transparent') 

 

sdev1=sort(abs(epiMT_deviances)) 
sdev2=sort(abs(epiMV_deviances)) 
 
 
#plot(ppoints(59)[1:57],sdev1[1:57]
#, type="b",pch=1,xlab="quantiles", ylab= " deviance residuals")
#lines(ppoints(59)[1:57],sdev2[1:57], type="b",pch=2)
 
#title("Epilepsy: quantiles smaller than 57/59")
 
#legend(x="topleft",legend=c("MT","ML"), pch=c(1,2))

plot(ppoints(59)[1:48],sdev1[1:48]
, type="b",pch=1,xlab="quantiles", ylab= " deviance residuals")
lines(ppoints(59)[1:48],sdev2[1:48], type="b",pch=2)
 
legend(x="topleft",legend=c("MT","ML"), pch=c(1,2))

#dev.off()

 


 

