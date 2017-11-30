
#This function compute a weighted likelihood estimator for the logistic model, where the weights penalize high leverage observations. In this version the weights are zero or one
#INPUT
#x0: covariable pxn matrix, p is the number of covariables, n the number of observations
#y:dimension n response vector
#intercept: 1 if the model include intercept, 0 if it does not include intercept

#OUTPUT
#$xweights n-dimensional  vector of zeros and ones used to compute the weighted maimum likelihood estimator
#$coefficients:    (p+intercept)- dimensional vector of     regression coefficients
#$standard.deviation:   vector of  the standard deviations of the regression coefficient estimators
#$fitted.values: a n dimensional  vector with the probabilities of success 
#$cov: a (p+intercept)x(p+intercept) covariance matrix of  the regression coefficients
#$residual.deviances: a n-dimensional  vector with the residual deviances



  WMLlogreg=function (x0, y, intercept = 1)
  {  
  ttx=colnames(x0)  
  if (!is.numeric(y)) 
        y <- as.numeric(y)
  if (!is.null(dim(y))) {
  if (ncol(y) != 1) 
            stop("y is not onedimensional")
        y <- as.vector(y)
    }
    n <- length(y)
  #  if (is.data.frame(x0)) {
  #      x0 <- data.matrix(x0)
  #  }
  #else if (!is.matrix(x0)) {
  #    x0 <- matrix(x0, length(x0), 1, dimnames = list(names(x0), 
  #       deparse(substitute(x0))))
  # }
   x0<-as.matrix(x0)
   if (nrow(x0) != n) 
        stop("Number of observations in x and y not equal")
   na.x <- !is.finite(rowSums(x0))
   na.y <- !is.finite(y)
   ok <- !(na.x | na.y)
   if (!all(ok)) {
   x0 <- x0[ok, , drop = FALSE]
   y <- y[ok]
    }
   n <- nrow(x0)  
   if (n == 0) 
   stop("All observations have missing values!")
   p<-ncol(x0)
   family <- binomial()
   p<-dim(x0)[2]
   zz<-rep(0,p)
   for (i in 1:p)
   {zz[i]<-sum((x0[,i]==0)|(x0[,i]==1))}
   tt<-which(zz!=n)	
   p1<-length(tt)
   x0=as.matrix(x0,n,p)
   x00<-x0[,tt]
   if(p1==1)
   {rdx<-abs(x00-median(x00))/mad(x00)
   wrd<-(rdx<=qnorm(.9875))}
   if(p1>1)
   { 
   mcdx<-covMcd(x00,alpha=.75)	
   rdx<-mahalanobis(x00,center=mcdx$center,cov=mcdx$cov)		
   vc<-qchisq(0.975,p)
   wrd<-(rdx<=vc)}
   if(p1==0)
   {wrd=1:n}
   
   x000=x0[wrd,]
   colnames(x000)=ttx
   y00=y[wrd] 

   if (intercept==1) 
   {
   out=glm(y00~x000, family = family)
    x <- cbind(Intercept = 1, x0)}
   if (intercept==0) 
   {out<-glm(y00~x000-1, family = family)
   x<-x0}
   out1<-summary(out)
   fitted.linear<-(x%*%out$coeff)
   fitted.linear<-pmin(fitted.linear,1e2)
   fitted.values<-exp(fitted.linear)/(1+exp(fitted.linear))
   residual.deviances<-sign(y-fitted.values)*sqrt(-2*(y*log(fitted.values)+(1-y)*log(1-fitted.values)))
   result<-list(xweights=wrd, coefficients=out$coeff,standard.deviation=sqrt(diag(out1$cov.uns)),fitted.values=t(fitted.values),cov=out1$cov.unscaled, residual.deviances= t(residual.deviances))
   result
 }

    
 
     
     