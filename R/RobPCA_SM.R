#' Robust principal components
#'
#' This function computes robust principal components based on the minimization of
#' the "residual" M-scale.
#'
#' @export SMPCA pcaRobS
#' @aliases SMPCA pcaRobS
#' @rdname SMPCA
#'
#' @param X a data matrix with observations in rows.
#' @param ncomp desired (maximum) number of components
#' @param desprop desired (minimum) proportion of explained variability (default = 0.9)
#' @param maxit maximum number of iterations (default= 100)
#' @param deltasca "delta" parameter of the scale M-estimator (default=0.5)
#'
#' @return A list with the following components:
#' \item{q}{The actual number of principal components}
#' \item{propex}{The actual proportion of unexplained variability}
#' \item{eigvec}{Eigenvectors, in a \code{p x q} matrix}
#' \item{fit}{an \code{n x p} matrix with the rank-q approximation to \code{X}}
#' \item{repre}{An \code{n x q} matrix with representation of data in R^q (scores)}
#' \item{propSPC}{A vector of length \code{p} with the cumulative explained variance from initial SPC}
#'
#' @author Ricardo Maronna, \email{rmaronna@retina.ar}, based on original code
#' by D. Pen~a and J. Prieto
#'
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(bus)
#' X0 <- as.matrix(bus)
#' X1 <- X0[,-9]
#' ss <- apply(X1, 2, mad)
#' mu <- apply(X1, 2, median)
#' X <- scale(X1, center=mu, scale=ss)
#' q <- 3  #compute three components
#' rr <- pcaRobS(X, q, 0.99)
#' round(rr$eigvec, 3)
#'
pcaRobS <- SMPCA <-function(X, ncomp, desprop=0.9, deltasca=0.5, maxit=100) {

  Wf <- function(r){ return((1-r)^2*(r<=1)) }    #Bisquare weights for r=resid^2

  n=dim(X)[1];  p=dim(X)[2]
  tol=1.e-4;
  #Inicial
  sp= rrcov::PcaLocantore(X, k=p)
  mu0=sp@center; lamL=sp@eigenvalues
  lamL=lamL/sum(lamL);  propSPC=cumsum(lamL)
  q1=sum(propSPC<desprop)+1
  q=min(c(q1,ncomp))

  QL=sp@loadings[,1:q]; Xcen=scale(X, center=mu0, scale=FALSE)
  fitL=scale(Xcen%*%QL%*%t(QL), center=-mu0,scale=FALSE)
  #sigini es para prop. de var. explicada
  rr=colSums(t(Xcen)^2);
  sigini=mscale(sqrt(rr), delta=deltasca,  tuning.chi = 1, family='bisquare')^2
  #escala inicial
  rr=colSums(t(X-fitL)^2);
  sig0=mscale(sqrt(rr), delta=deltasca,  tuning.chi = 1, family='bisquare')^2
  ww=Wf(rr/sig0);  B0=QL

  iter=0; del=100;
  while (iter<maxit & abs(del)>tol) {
    iter=iter+1;
    #update mu
    mu=colSums(X*ww)/sum(ww)
    Xcen=scale(X, center=mu, scale=FALSE)
    #update B
    C=t(Xcen)%*%diag(ww)%*%Xcen
    B <- svd(C, nu=q, nv=q)$u 	# B=rARPACK::eigs_sym(C, q)$vectors
    fit=Xcen%*%B%*%t(B)
    res=Xcen-fit   #residuals
    rr=colSums(t(res)^2)
    sig=mscale(sqrt(rr),delta=deltasca,  tuning.chi = 1, family='bisquare') ^2
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
  resu=list(eigvec=B, fit=fit, repre=repre, propex=propex, propSPC=propSPC, mu=mu, q=q)
  return(resu)
}
