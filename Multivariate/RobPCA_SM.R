library("rrcov")
library("rARPACK")
options(digits=4)


SMPCA<-function(X, ncomp, desprop=0.9, deltasca=0.5, maxit=100) {
# Robust principal components based on minimization of "residual" M-scale (Maronna 2005)
# INPUT
# X=data matrix
# ncomp=desired (maximum) number of components
# desprop= desired (minimum) proportion of unexplained variability (default= 0.9)
# maxit=maximum number of iterations (default= 100)
# deltasca= "delta" parameter of scale M-estimator (default=0.5)  
# OUTPUT
# q=actual number of components
# propex=actual proportion  of unexplained variability
#  q=min(ncomp,q1), where q1 is the minimum number of components which makes propex>=desprop
# eigvec=eigenvectors (pxq)
# fit= (nxp) rank-q approximation to X
# repre=(nxq) representation of data in R^q (scores)
# propSPC (p-vector): cumulative explained variance from initial SPC


n=dim(X)[1];  p=dim(X)[2]
tol=1.e-4; 
#Inicial
sp=PcaLocantore(X, k=p)
mu0=sp@center; lamL=sp@eigenvalues
lamL=lamL/sum(lamL);  propSPC=cumsum(lamL)
q1=sum(propSPC<desprop)+1
q=min(c(q1,ncomp))

QL=sp@loadings[,1:q]; Xcen=scale(X, center=mu0, scale=FALSE)
fitL=scale(Xcen%*%QL%*%t(QL), center=-mu0,scale=FALSE) 
#sigini es para prop. de var. explicada
rr=colSums(t(Xcen)^2); 
sigini=mscale(sqrt(rr),delta=deltasca,  tuning.chi = 1)^2 
#escala inicial
rr=colSums(t(X-fitL)^2); 
sig0=mscale(sqrt(rr),delta=deltasca,  tuning.chi = 1)^2 
 ww=Wf(rr/sig0);  B0=QL
            
 iter=0; del=100;  
while (iter<maxit & abs(del)>tol) {
	iter=iter+1;     
#update mu	
	mu=colSums(X*ww)/sum(ww)
	Xcen=scale(X, center=mu, scale=FALSE)
#update B
	C=t(Xcen)%*%diag(ww)%*%Xcen
	B=eigs_sym(C, q)$vectors
	fit=Xcen%*%B%*%t(B)
	res=Xcen-fit   #residuals
	rr=colSums(t(res)^2)
	sig=mscale(sqrt(rr),delta=deltasca,  tuning.chi = 1) ^2
	del1=1-sig/sig0;  #print(c(iter,del1))
	U=diag(q)-abs(t(B)%*%B0)
	del2=mean(abs(U)); del=max(c(del1,del2))
	sig0=sig; B0=B
	ww=Wf(rr/sig);               
	repre=Xcen%*%B  #represnt. in R^2
	#	print(c(iter, del1, del2))
} #endo
propex=1-sig/sigini  #prop. var. expli
fit=scale(fit, center=-mu, scale=FALSE)
resu=list(eigvec=B, fit=fit, repre=repre, propex=propex,propSPC=propSPC,mu=mu)
return(resu)
}

Wf<-function(r){ return((1-r)^2*(r<=1))}    #Bisquare weights for r=resid^2