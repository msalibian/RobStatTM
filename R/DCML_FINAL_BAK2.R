

require(robustbase) # require(lmrob)
require(pyinit) # require(penseinit)
require(quantreg)


mscale=function(u, ep,delta)
#M scaleof a sample u  
#ep accuracy
#delta breakdown point (right side)
{s0=median(abs(u))/.6745
h=1
it=0
while((h>ep)&(it<50))
{       it=it+1
s1=(1/delta)*(s0^2)*mean(rho.bi((u/s0),1.547))
s1=s1^(1/2)
h=abs(s1-s0)/s0
s0=s1}
s0}



covdcml=function(resLS,res0,C,sig0,t0,p,n)
##Computation of the asymptotic covariance matrix of the DCML estimator
{t0=1-t0
c=mean(psibs(res0/sig0,3.44)*resLS)
a1=mean(psibs(res0/sig0,3.44)^2)
b=mean(psibspri(res0/sig0,3.44))
deltasca=0.5*(1-(p/n))
a2=mscale(resLS,.00001,deltasca)^2
tuti=t0^2*sig0^2*a1/b^2 + a2*(1-t0)^2 +2*t0*(1-t0)*sig0*c/b
V=tuti*solve(C)
V
}



  psibs=function(t,c)  
# bisquare psi function
 {r=t/c
 gg=r*(1-r^2)^2*(abs(r)<=1)
 gg
}
        
 
 psibspri=function (t,c)   
#Derivative of the bisquare psi function
{ r=t/c
 gg=(1-r^2)*(1-5*r^2)*(abs(r)<=1)/c
gg
}
        
 MMPY=function(X,y,  intercept=TRUE)
#Compute an MM-estimator taking as initial Pe?a Yohai
#INPUT
#X nxp matrix, where n is the number of observations and p the number of  columns, the 1's of the intercept are not included
#y vector of dimension  n with the responses
#intercept FALSE the regression does not includes intercept, TRUE, the regression includes intercept
#OUTPUT
#outMM output of the MM estimator (lmrob) with 85% of efficiency and PY as initial
{
cont1=lmrob.control(tuning.psi=3.44)
n=nrow(X)
p=ncol(X)
if(intercept==TRUE )
{ p=p+1}
dee=.5*(1-(p/n))
a <- pyinit(intercept,X=X, y=y, deltaesc=dee, cc.scale=1.547, prosac=.5,clean.method="threshold", C.res = 2, prop=.2, py.nit = 20, en.tol=1e-5)
betapy2=a$initCoef[,1]
sspy2=a$objF[1]
uu=list(coeff=betapy2,scale=sspy2)
 if(intercept==TRUE)
    { outMM=lmrob(y~X, control=cont1,init=uu)}else
    { outMM=lmrob(y~X-1, control=cont1,init=uu)}
outMM
}


DCML_FINAL=function(X,y, outMM, intercept=TRUE)
#INPUT
#X nxp matrix, where n is the number of observations and p the number of  columns, the 1's of the intercept are not included
#y vector of dimension  n with the responses
#intercept FALSE the regression does not includes intercept, TRUE, the regression includes intercept
#outMM output of MMPY or SM_PY
#OUTPUT
#coef vector of coefficients, the first element is the intercept when there is one
#cov covariance matrix of the coefficients
#resid vector of residuals
#weight  vector with weights that the MMestimator assigns to every observation
#sigma standad error of errot term
 

{
XX=X
n=nrow(X)
p=ncol(X)
if(intercept==TRUE )
    {XX=cbind(rep(1,n),X)
    p=p+1}
cont1=lmrob.control(tuning.psi=3.44)
if(intercept==TRUE)
   {out3=lm(y~X)}else
   {out3=lm(y~X-1)}
betaLS=out3$coeff
resLS=out3$resid
dee=.5*(1-(p/n))
beta0=outMM$coeff
weight=outMM$rweights
residuos=outMM$resid
sigma=mscale(residuos,.00001,dee)
##Begin the computation of the DCML
deltas=.3*p/n 
C=t(XX)%*%diag(weight)%*%XX/sum(weight) 	
d=(beta0-betaLS)%*%C%*%(beta0-betaLS)
d=d/sigma^2
t0=min(1,sqrt(deltas/d))
beta1=t0*betaLS+(1-t0)*beta0
V=covdcml (resLS,residuos,C,sigma,t0,p,n)/n	
resid=y-XX%*%beta1
sigma=mscale(resid, .00001,dee) 
out=list(coef=beta1, cov=V, resid=resid,   sigma=sigma )
out
}
SM_PY= function(y,X,Z, intercept=TRUE)

{
cont1=lmrob.control(tuning.psi=3.44) 
n=nrow(X)
q=ncol(Z)
p=ncol(X)
hh1=1+intercept
hh2=p+intercept
hh3=q+intercept
 gamma=matrix(0, (q+intercept),p)
 
  
for( i in 1:p)
{gamma[,i]=rq(X[,i]~Z)$coeff}
X1=X-Z%*%gamma[hh1:hh3,]
#aa=c(dim(Z),dim(X),dim(X1),dim(y))
 

dee=.5*(1-((p+1)/n))
out= pyinit(intercept=intercept,X=X1, y=y, deltaesc=dee, cc.scale=1.547, prosac=.5,clean.method="threshold", C.res = 2, prop=.2, py.nit = 20, en.tol=1e-5)
  betapy=out$initCoef[,1]
     sspy=out$objF[1]
     uu=list(coeff=betapy,scale=sspy)
     out0=lmrob(y~X1, control=cont1,init=uu) 
     beta=out0$coeff

y1=y-X1%*%beta[hh1:hh2]


fi=rq(y1~Z)$coeff
oo=NULL
if(intercept==TRUE)
oo=fi[1]
#print(length(beta))
#print(dim(gamma[hh1:hh3,]))
tt=gamma[hh1:hh3,]
if(p==1)
tt=matrix(tt,q,p)
beta00=c(oo,beta[hh1:hh2],fi[hh1:hh3]-  tt%*%beta[hh1:hh2])
res=y1-fi[1]-Z%*%fi[hh1:hh3]
res=as.vector(res)
 dee=.5*(1-((p+q+intercept)/n))
ss=mscale(res,.0001,dee)
uu=list(coeff=beta00,scale=ss)
XX=cbind(X,Z)
beta0=lmrob(y~XX,control=cont1,init=uu)$coeff
beta0
}


rho.bi=function(x,k)
{
# bicuadratic rho function
n=length(x)
dif=(abs(x)-k<0)
ro1=vector("numeric",n)
ro1[dif]=(x[dif]^2)/2-(x[dif]^4)/(2*(k^2))+(x[dif]^6)/(6*(k^4))
ro1[!dif]=k^2/6
ro1=ro1*6/(k^2)
ro1
}





 
 
