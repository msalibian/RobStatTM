#source("MScale.R")

MLocDis<- function(x,psi="Bis",eff=0.9, maxit=50, tol=1.e-4) {
  #M-estimators of location and dispersion
  #psi= "Bis" (quare) (default) of "Hub" (er)
  #eff=desired asymptotic efficiency. Options: 0.9 (default), 0.85 or 0.95
  #Output: mu=location estimator, std.mu= estimated standard deviation of mu, 
  # disper= M-dispersion estimator
  ########
  if (psi=="Bis") {kpsi=1
  } else  if (psi=="Hub") {kpsi=2
  } else {print(c(psi, " No such psi")); kpsi=0
  }
  kBis=c(3.44, 3.88, 4.685)
  kHub=c(0.732, 0.981, 1.34)
  kk=rbind(kBis, kHub)
  efis=c(0.85, 0.90, 0.95)
  if (is.element(eff, efis)) {keff=match(eff,efis); 
  } else {print(c(eff, " No such eff")); keff=0}
if (kpsi>0 & keff>0) {
  ktun=kk[kpsi, keff]
  mu0=median(x); sig0=mad(x)
  if (sig0<1.e-10) {mu=0; sigma=0
  } else { #initialize
    dife=1.e10; iter=0
    while (dife>tol & iter<maxit) {
      iter=iter+1
      resi=(x-mu0)/sig0; ww=wfun(resi/ktun, kpsi)
      mu=sum(ww*x)/sum(ww)
      dife=abs(mu-mu0)/sig0; mu0=mu
    } # end while
  } # end if sig
} #end if k
  rek=resi/ktun; pp=psif(rek, kpsi)*ktun
  n=length(x)
  a=mean(pp^2); b=mean(psipri(rek, kpsi))
  sigmu=sig0^2 *a/(n*b^2)
  sigmu=sqrt(sigmu)
  scat=mscale(x-mu)
  resu=list(mu=mu, std.mu=sigmu, disper=scat)
  return(resu)
} # end function

wfun<- function(x,k) { #weight function
  if (k==1) ww=(1-x^2)^2 *(abs(x)<=1)
  else  ww=(abs(x)<=1)+(abs(x)>1)/(abs(x)+1.e-20)
  return(ww)
}

psif<-function(x,k) {return(x*wfun(x,k))}

psipri<-function(x,k) {
  if (k==1) pp=	(((1 - (x^2))^2) - 4 * (x^2) * (1 - (x^2))) * (abs(x) < 1)
  else pp=(abs(x)<=1)
  return(pp)
}