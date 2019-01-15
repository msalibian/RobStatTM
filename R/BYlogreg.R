#' Bianco and Yohai estimator for logistic regression
#'
#' This function computes the M-estimator proposed by Bianco and Yohai for
#' logistic regression. By default, an intercept term is included and p
#' parameters are estimated. Modified by Yohai (2018) to take as initial estimator
#' a weighted ML estimator with weights derived from the MCD estimator.
#' For more details we refer to Croux, C., and Haesbroeck, G. (2002),
#' "Implementing the Bianco and Yohai estimator for Logistic Regression"
#'
#' @export BYlogreg logregBY
#' @aliases BYlogreg logregBY
#' @rdname BYlogreg
#'
#' @param x0 matrix of explanatory variables;
#' @param y vector of binomial responses (0 or 1);
#' @param intercept 1 or 0 indicating if an intercept is included or or not
#' @param const tuning constant used in the computation of the estimator (default=0.5);
#' @param kmax maximum number of iterations before convergence (default=1000);
#' @param maxhalf max number of step-halving (default=10).
#'
#' @return A list with the following components:
#' \item{coefficients}{estimates for the regression coefficients}
#' \item{standard.deviation}{standard deviations of the coefficients}
#' \item{fitted.values}{fitted values}
#' \item{residual.deviances}{residual deviances}
#' \item{components}{logical value indicating whether convergence was achieved}
#' \item{objective}{value of the objective function at the minimum}
#'
#' @author Christophe Croux, Gentiane Haesbroeck, Victor Yohai
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(skin)
#' Xskin <- as.matrix( skin[, 1:2] )
#' yskin <- skin$vasoconst
#' skinBY <- logregBY(Xskin, yskin, intercept=1)
#' skinBY$coeff
#' skinBY$standard.deviation
#'
logregBY <- BYlogreg <- function(x0,y, intercept=1, const=0.5,kmax=1000,maxhalf=10)
{
  sigmamin=0.0001
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
    mcdx <- robustbase::covMcd(x00,alpha=.75)
    rdx <- mahalanobis(x00,center=mcdx$center,cov=mcdx$cov)
    vc<-qchisq(0.975,p)
    wrd<-(rdx<=vc)}
  if(p1==0)
  {wrd=1:n}

  if (intercept==1)
  {out<-glm(y[wrd]~x0[wrd,], family = family)
  x <- cbind(Intercept = 1, x0)}
  if (intercept==0)
  {out<-glm(y[wrd]~x0[wrd,]-1, family = family)
  x<-x0}
  gstart<-out$coeff
  p=ncol(x)
  sigmastart=1/sqrt(sum(gstart^2))
  xistart=gstart*sigmastart
  stscores=x %*% xistart


  #Initial value for the objective function

  oldobj=mean(phiBY3(stscores/sigmastart,y,const))
  kstep=1
  jhalf=1

  while (( kstep < kmax) & (jhalf<maxhalf))
  {

    optimsig=optimize(sigmaBY3,lower=0,upper=10^3,y=y,s=stscores,c3=const)
    sigma1=optimsig$minimum


    if (sigma1<sigmamin)
    {
      print("Explosion");kstep=kmax
    }
    else
    {


      gamma1=xistart/sigma1
      scores=stscores/sigma1
      newobj=mean(phiBY3(scores,y,const))
      oldobj=newobj

      gradBY3=apply((derphiBY3(scores,y,const)%*%t(rep(1,p)))*x,2,mean)

      h=-gradBY3+(sum(gradBY3*xistart) *xistart)
      finalstep=h/sqrt(sum(h^2))
      xi1=xistart+finalstep
      xi1=xi1/(sum(xi1^2))
      scores1=(x%*%xi1)/sigma1
      newobj=mean(phiBY3(scores1,y,const))

      ####stephalving

      hstep=1
      jhalf=1
      while  ((jhalf <=maxhalf) & (newobj>oldobj))
      {
        hstep=hstep/2
        xi1=xistart+finalstep*hstep
        xi1=xi1/sqrt(sum(xi1^2))
        scores1=x%*%xi1/sigma1
        newobj=mean(phiBY3(scores1,y,const))
        jhalf=jhalf+1
      }

      if ((jhalf==maxhalf+1) & (newobj>oldobj))

      { #print("Convergence Achieved")
      } else
      { jhalf=1
      xistart=xi1
      oldobj=newobj
      stscores=x%*% xi1
      kstep=kstep+1
      }

    }

  }


  if (kstep == kmax)
  {print("No convergence")
    resultat=list(convergence=F,objective=0,coef=t(rep(NA,p)))

    resultat}

  else
  {
    gammaest=xistart/sigma1
    stander=sterby3(x,y,const,gammaest)

    fitted.linear<- x%*%gammaest
    fitted.linear<-pmin(fitted.linear,1e2)
    fitted.values<-exp(fitted.linear)/(1+exp(fitted.linear))
    residual.deviances<-sign(y-fitted.values)*sqrt(-2*(y*log(fitted.values)+(1-y)*log(1-fitted.values)))
    result <-list(convergence=TRUE,objective=oldobj, coefficients=gammaest, standard.deviation=stander,fitted.values=t(fitted.values),  residual.deviances= t(residual.deviances))



    result
  }




}




###############################################################
###############################################################
#Functions needed for the computation of estimator of Bianco and Yohai

phiBY3=function(s,y,c3)
{
  s=as.double(s)
  dev=log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
  res=rhoBY3(dev,c3)+GBY3Fs(s,c3)+GBY3Fsm(s,c3)
  res

}


rhoBY3 <- function(t,c3)
{
  (t*exp(-sqrt(c3))*as.numeric(t <= c3))+
    (((exp(-sqrt(c3))*(2+(2*sqrt(c3))+c3))-(2*exp(-sqrt(t))*(1+sqrt(t))))*as.numeric(t >c3))
}

psiBY3 <- function(t,c3)
{
  (exp(-sqrt(c3))*as.numeric(t <= c3))+(exp(-sqrt(t))*as.numeric(t >c3))

}


derpsiBY3 <- function(t,c3)
{
  (0*as.numeric(t <= c3))+(-(exp(-sqrt(t))/(2*sqrt(t)))*as.numeric(t >c3))

}

sigmaBY3<-function(sigma,s,y,c3)
{
  mean(phiBY3(s/sigma,y,c3))
}



derphiBY3=function(s,y,c3)
{
  Fs=	exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))

  ds=Fs*(1-Fs)
  dev=log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
  Gprim1=log(1+exp(-abs(s)))+abs(s)*(s<0)
  Gprim2=log(1+exp(-abs(s)))+abs(s)*(s>0)
  -psiBY3(dev,c3)*(y-Fs)+((psiBY3(Gprim1,c3)-psiBY3(Gprim2,c3))*ds)
}

der2phiBY3=function(s,y,c3)
{
  s=as.double(s)
  Fs=	exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
  ds=Fs*(1-Fs)
  dev=log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
  Gprim1=log(1+exp(-abs(s)))+abs(s)*(s<0)
  Gprim2=log(1+exp(-abs(s)))+abs(s)*(s>0)
  der2=(derpsiBY3(dev,c3)*(Fs-y)^2)+(ds*psiBY3(dev,c3))
  der2=der2+(ds*(1-2*Fs)*(psiBY3(Gprim1,c3)-psiBY3(Gprim2,c3)))
  der2=der2-(ds*((derpsiBY3(Gprim1,c3)*(1-Fs))+(derpsiBY3(Gprim2,c3)*Fs)))
  der2
}


GBY3Fs <- function(s,c3)
{
  Fs=	exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
  resGinf=exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fs))))-1)
  resGinf=(resGinf+(Fs*exp(-sqrt(-log(Fs)))))*as.numeric(s <= -log(exp(c3)-1))
  resGsup=((Fs*exp(-sqrt(c3)))+(exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1)))*as.numeric(s > -log(exp(c3)-1))
  resG=resGinf+resGsup
  resG
}



GBY3Fsm <- function(s,c3)
{
  Fsm=exp(-(log(1+exp(-abs(s)))+abs(s)*(s>0)))
  resGinf=exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fsm))))-1)
  resGinf=(resGinf+(Fsm*exp(-sqrt(-log(Fsm)))))*as.numeric(s >= log(exp(c3)-1))
  resGsup=((Fsm*exp(-sqrt(c3)))+(exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1)))*as.numeric(s < log(exp(c3)-1))
  resG=resGinf+resGsup
  resG
}

sterby3 <- function(z,y,const,estim)

{
  n=dim(z)[[1]]
  p=dim(z)[[2]]


  argum=z %*% estim

  matM=matrix(data=0,nrow=p,ncol=p)
  for (i in 1:n)
  {
    matM=matM+(der2phiBY3(argum[i],y[i],const) * (z[i,] %*% t(z[i,])))
  }

  matM=matM/n
  matMinv=solve(matM)

  IFsquar=matrix(data=0,nrow=p,ncol=p)
  for (i in 1:n)
  {
    IFsquar=IFsquar+((derphiBY3(argum[i],y[i],const))^2 * (z[i,] %*% t(z[i,])))
  }

  IFsquar=IFsquar/n
  asvBY=matMinv %*% IFsquar %*% t(matMinv)
  sqrt(diag(asvBY))/sqrt(n)
}
