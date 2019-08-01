#' Robust multivariate location and scatter estimators
#'
#' This function computes robust estimators for multivariate location and scatter.
#'
#' This function computes robust estimators for multivariate location and scatter.
#' The default behaviour (\code{type = "auto"}) computes a "Rocke" estimator
#' (as implemented in \code{\link{covRobRocke}}) if the number
#' of variables is greater than or equal to 10, and an MM-estimator with a
#' SHR rho function (as implemented in \code{\link{covRobMM}}) otherwise.
#'
#' @export Multirobu covRob
#' @aliases Multirobu covRob
#' @rdname Multirobu
#'
#' @param X a data matrix with observations in rows.
#' @param type a string indicating which estimator to compute. Valid options
#' are "Rocke" for Rocke's S-estimator, "MM" for an MM-estimator with a
#' SHR rho function, or "auto" (default) which selects "Rocke" if the number
#' of variables is greater than or equal to 10, and "MM" otherwise.
#' @param maxit Maximum number of iterations, defaults to 50.
#' @param tol Tolerance for convergence, defaults to 1e-4.
#'
#' @return A list with class \dQuote{covClassic} with the following components:
#' \item{mu}{The location estimate}
#' \item{V}{The scatter matrix estimate, scaled for consistency at the normal distribution}
#' \item{center}{The location estimate. Same as \code{mu} above.}
#' \item{cov}{The scatter matrix estimate, scaled for consistency at the normal distribution. Same as \code{V} above.}
#' \item{dist}{Robust Mahalanobis distances}
#'
#' @author Ricardo Maronna, \email{rmaronna@retina.ar}
#'
#' @seealso \code{\link{covRobRocke}}, \code{\link{covRobMM}}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(bus)
#' X0 <- as.matrix(bus)
#' X1 <- X0[,-9]
#' tmp <- covRob(X1)
#' round(tmp$cov[1:10, 1:10], 3)
#' tmp$mu
#'
covRob <- Multirobu <- function(X, type="auto", maxit=50, tol=1e-4)  {
if (type=="auto") {
  p=dim(X)[2]
  if (p<10) {type="MM"
  } else {type="Rocke"}
}

 if (type=="Rocke") {
   resu=RockeMulti(X, maxit=maxit, tol=tol)
 } else {resu=MMultiSHR(X, maxit=maxit, tolpar=tol)  #MM
 }
  mu=resu$mu; V=resu$V

  # Feed list into object and give class
  z <- list(mu=mu, V=V, dist=mahalanobis(X,mu,V), cov=V, center=mu)
  class(z) <- c("covRob")
  return(z)
}

#' Rocke's robust multivariate location and scatter estimator
#'
#' This function computes Rocke's robust estimator for multivariate location and scatter.
#'
#' This function computes Rocke's robust estimator for multivariate location and scatter.
#'
#' @export RockeMulti covRobRocke
#' @aliases RockeMulti covRobRocke
#' @rdname RockeMulti
#'
#' @param X a data matrix with observations in rows.
#' @param initial A character indicating the initial estimator. Valid options are 'K' (default)
#' for the Pena-Prieto 'KSD' estimate, and 'mve' for the Minimum Volume Ellipsoid.
#' @param maxsteps Maximum number of steps for the line search section of the algorithm.
#' @param propmin Regulates the proportion of weights computed from the initial estimator that
#' will be different from zero. The number of observations with initial non-zero weights will
#' be at least p (the number of columns of X) times propmin.
#' @param qs Tuning paramater for Rocke's loss functions.
#' @param maxit Maximum number of iterations.
#' @param tol Tolerance to decide converngence.#'
#'
#' @return A list with class \dQuote{covRob} containing the following elements:
#' \item{mu}{The location estimate}
#' \item{V}{The scatter matrix estimate, scaled for consistency at the normal distribution}
#' \item{center}{The location estimate. Same as \code{mu} above.}
#' \item{cov}{The scatter matrix estimate, scaled for consistency at the normal distribution. Same as \code{V} above.}
#' \item{dista}{Robust Mahalanobis distances}
#' \item{w}{weights}
#' \item{gamma}{Final value of the constant gamma that regulates the efficiency}
#'
#' @author Ricardo Maronna, \email{rmaronna@retina.ar}
#'
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(bus)
#' X0 <- as.matrix(bus)
#' X1 <- X0[,-9]
#' tmp <- covRobRocke(X1)
#' round(tmp$cov[1:10, 1:10], 3)
#' tmp$mu
#'
covRobRocke <- RockeMulti <- function(X, initial='K', maxsteps=5, propmin=2, qs=2, maxit=50, tol=1e-4)
{
  d <- dim(X)
  n <- d[1]
  p <- d[2]


  gamma0 <- consRocke(p=p, n=n, initial )$gamma # tuning constant
  if(initial=='K')
{  out=KurtSDNew(X)
   mu0=out$center; V0=out$cova
   V0=V0/(det(V0)^(1/p))
dista0=mahalanobis(X,mu0,V0)
dista=dista0}

if(initial=='mve')
{out=fastmve(X)

mu0=out$center
    V0=out$cov
   V0=V0/(det(V0)^(1/p))
dista0=mahalanobis(X,mu0,V0)
dista=dista0}


  delta <- (1-p/n)/2 # max breakdown
  #gamma0 <- consRocke(p,n,'K')$gamma
  sig <- MScalRocke(x=dista, gamma=gamma0, q=qs, delta=delta) #Inicializar
  # %Buscar gama que asegure que al menos p*propmin elementos tengan w>0
  didi <- dista / sig
  dife <- sort( abs( didi - 1) )
  gg <- min( dife[ (1:n) >= (propmin*p) ] )
  gamma <- max(gg, gamma0)
#print(gamma)
  sig0 <- MScalRocke(x=dista, gamma=gamma, delta=delta, q=qs)

  iter <- 0
  difpar <- difsig <- +Inf
  while( ( ( (difsig > tol) | (difpar > tol) ) &
           (iter < maxit) ) & (difsig > 0) ) {
    iter <- iter + 1
    w <- WRoTru(tt=dista/sig, gamma=gamma, q=qs)
    mu <- colMeans( X * w ) / mean(w) # as.vector( t(w) %*% X ) / sum(w)
    Xcen <- scale(X, center=mu, scale=FALSE)
    V <- t(Xcen) %*% (Xcen * w) / n;
    V <- V / ( det(V)^(1/p) )
    dista <- mahalanobis(x=X, center=mu, cov=V)

    sig <- MScalRocke(x=dista, gamma=gamma, delta=delta, q=qs)
    # %Si no desciende, hacer Line search
    step <- 0
    delgrad <- 1
    while( (sig > sig0) & (step < maxsteps) ) {
      delgrad <- delgrad / 2
      step <- step + 1
      mu <- delgrad * mu + (1 - delgrad)*mu0
      V <- delgrad*V + (1-delgrad)*V0
      V <- V / ( det(V)^(1/p) )
      dista <- mahalanobis(x=X, center=mu, cov=V)

      sig <- MScalRocke(x=dista, gamma=gamma, delta=delta, q=qs)
    }
    dif1 <- as.vector( t(mu - mu0) %*% solve(V0, mu-mu0) ) / p
    dif2 <- max(abs(solve(V0, V)-diag(p)))
    difpar <- max(dif1, dif2)
    difsig <- 1 - sig/sig0
    mu0 <- mu
    V0 <- V
    sig0 <- sig
  }
  tmp <- scalemat(V0=V0, dis=dista, weight='X')
  V <- tmp$V
  ff <- tmp$ff
  dista <- dista/ff

  # GSB: Feed list into object and give class
  z <- list(mu=mu, V=V, sig=sig, dista=dista, w=w, gamma=gamma, cov=V, center=mu)
  class(z) <- c("covRob")
  return(z)
}

consRocke <- function(p, n, initial) {
  if(initial=='M') {
    beta <- c(-5.4358, -0.50303, 0.4214)
  } else {
    beta <- c(-6.1357, -1.0078, 0.81564)
  }
  if( p >= 15 ) {
    a <- c(1, log(p), log(n))
    alpha <- exp( sum( beta*a ) )
    gamma <- qchisq(1-alpha, df=p)/p - 1
    gamma <- min(gamma, 1)
  } else {
    gamma <- 1
    alpha <- 1e-6
  }
  return(list(gamma=gamma, alpha=alpha))
}


WRoTru <- function(tt, gamma, q) {
  ss <- (tt - 1) / gamma
  w <- 1 - ss^q
  w[ abs(ss) > 1 ] <- 0
  return(w)
}


rhorotru <- function(tt, gamma, q) {
  u <- (tt - 1) / gamma
  y <- ( (u/(2*q)*(q+1-u^q) +0.5) )
  y[ u >= 1 ] <- 1
  y[ u < (-1) ] <- 0
  return(y)
}



MScalRocke <- function(x, gamma, q, delta = 0.5, tol=1e-5)
{
  # sigma= solucion de ave{rhorocke1(x/sigma)}=delta
  n <- length(x)
  y <- sort(abs(x))
  n1 <- floor(n * (1-delta) )
  n2 <- ceiling(n * (1 - delta) / (1 - delta/2) )
  qq <- y[c(n1, n2)]
  u <- 1 + gamma*(delta-1)/2 #asegura rho(u)<delta/2
  sigin <- c(qq[1]/(1+gamma), qq[2]/u)
  if( qq[1] >= 1) {
    tolera <- tol
  } else {
    tolera <- tol * qq[1]
  }
  if( mean(x==0) > (1-delta) ) {
    sig <- 0
  } else {
    sig <- uniroot(f=averho, interval=sigin, x=x,
                   gamma=gamma, delta=delta, q=q, tol=tolera)$root
  } # solucion de ave{rhorocke1(x/sigma)}=delta
  return(sig)
}

averho <- function(sig, x, gamma, delta, q)
  return( mean( rhorotru(x/sig, gamma, q) ) - delta )


scalemat <- function(V0, dis, weight='X')
{
  p <- dim(V0)[1]
  if( weight == 'M') {
    sig <- M_Scale(x=sqrt(dis), normz=0)^2
    cc <- 4.8421*p-2.5786 #ajuste empirico
  } else {
    sig <- median(dis)
    cc <- qchisq(0.5, df=p)
  }
  ff <- sig/cc
  return(list(ff=ff, V=V0*ff))
}


M_Scale <- function(x, normz=1, delta=0.5, tol=1e-5)
{
  n <- length(x)
  y <- sort(abs(x))
  n1 <- floor(n*(1-delta))
  n2 <- ceiling(n*(1-delta)/(1-delta/2));
  qq <- y[c(n1, n2)]
  u <- rhoinv(delta/2)
  sigin <- c(qq[1],  qq[2]/u) # intervalo inicial
  if (qq[1]>=1) {
    tolera=tol
  } else {
    tolera = tol * qq[1]
  }
  #tol. relativa o absol. si sigma> o < 1
  if( mean(x==0) >= (1-delta) ) {
    sig <- 0
  } else {
    sig <- uniroot(f=averho.uni, interval=sigin, x=x,
                   delta=delta, tol=tolera)$root
  }
  if( normz > 0) sig <- sig / 1.56
  return(sig)
}



rhobisq <- function(x) {
  r <- 1 - (1-x^2)^3
  r[ abs(x) > 1 ] <- 1
  return(r)
}

averho.uni <- function(sig, x, delta)
  return( mean( rhobisq(x/sig) ) - delta )

rhoinv <- function(x)
  return(sqrt(1-(1-x)^(1/3)))

#' MM robust multivariate location and scatter estimator
#'
#' This function computes an MM robust estimator for multivariate location and scatter with the "SHR" loss function.
#'
#' This function computes an MM robust estimator for multivariate location and scatter with the "SHR" loss function.
#'
#' @export MMultiSHR covRobMM
#' @aliases MMultiSHR covRobMM
#' @rdname MMultiSHR
#'
#' @param X a data matrix with observations in rows.
#' @param maxit Maximum number of iterations.
#' @param tolpar Tolerance to decide converngence.#'
#'
#' @return A list with class \dQuote{covRob} containing the following elements
#' \item{mu}{The location estimate}
#' \item{V}{The scatter matrix estimate, scaled for consistency at the normal distribution}
#' \item{center}{The location estimate. Same as \code{mu} above.}
#' \item{cov}{The scatter matrix estimate, scaled for consistency at the normal distribution. Same as \code{V} above.}
#' \item{dista}{Robust Mahalanobis distances}
#'
#' @author Ricardo Maronna, \email{rmaronna@retina.ar}
#'
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(bus)
#' X0 <- as.matrix(bus)
#' X1 <- X0[,-9]
#' tmp <- covRobMM(X1)
#' round(tmp$cov[1:10, 1:10], 3)
#' tmp$mu
#'
covRobMM <- MMultiSHR <- function(X, maxit=50, tolpar=1e-4) {
  d <- dim(X)
  n <- d[1]; p <- d[2]
  delta <- 0.5*(1-p/n) #max. breakdown
  cc <- consMMKur(p,n)
  const <- cc[1]
  out=KurtSDNew(X)
  mu0=out$center; V0=out$cova
  V0=V0/(det(V0)^(1/p))
  distin=mahalanobis(X,mu0,V0)

  sigma <- MscalSHR(t=distin, delta=delta)*const;
  iter <- 0
  difpar <- +Inf
  V0 <- diag(p)
  mu0 <- rep(0, p)
  dista <- distin
  while ((difpar>tolpar) & (iter<maxit)) {
    iter <- iter+1
    w <- weightsSHR(dista/sigma)
    mu <- colMeans(X * w) / mean(w) # mean(repmat(w,1,p).*X)/mean(w)
    Xcen <- scale(X, center=mu, scale=FALSE)
    V <- t(Xcen) %*% (Xcen * w) # Xcen'*(repmat(w,1,p).*Xcen);
    V <- V/(det(V)^(1/p)) # set det=1
    # V0in <- solve(V0);
    dista <- mahalanobis(x=X,center=mu, cov=V) # %rdis=sqrt(dista);
    dif1 <- t(mu - mu0) %*% solve(V0, mu-mu0) # (mu-mu0)*V0in*(mu-mu0)'
    dif2 <-  max(abs(c(solve(V0, V) - diag(p)))) # max(abs(V0in*V-eye(p)));
    difpar <- max(c(dif1, dif2)) #errores relativos de parametros
    mu0 <- mu
    V0 <- V
  }
  tmp <- scalemat(V0=V0, dis=dista, weight='X');
  dista <- dista/tmp$ff

  # GSB: Feed list into object and give class
  z <- list(V=tmp$V, mu=mu0, dista=dista, w=w, center=mu0, cov=tmp$V)
  class(z) <- c("covRob")
  return(z)
}


MscalSHR <- function(t, delta=0.5, sig0=median(t), niter=50, tol=1e-4) {
  if(mean(t<1e-16) >= 0.5) { sig <- 0
  } else { # fixed point
    sig1 <- sig0
    y1 <- meanrhoSHR(sig1, t, delta)
    # make meanrho(sig11, t) <= sig1
    while( (y1 > sig1) & ( (sig1/sig0) < 1000) ) {
      sig1 <- 2*sig1
      y1 <- meanrhoSHR(sig1, t, delta)
    }
    if( (sig1/sig0) >= 1000) { warning('non-convex function')
      sig <- sig0 } else {
        iter <- 0
        sig2 <- y1
        while( iter <= niter) {
          iter <- iter+1
          y2 <- meanrhoSHR(sig2, t, delta)
          den <- sig2-y2+y1-sig1 # secante
          if( abs(den) < tol*sig1) {
            sig3 <- sig2
            iter <- niter+1 } else {
              sig3 <- (y1*sig2-sig1*y2)/den
            }
          sig1 <- sig2
          sig2 <- sig3
          y1 <- y2
        }
        sig <- sig2
      }
  }
  return(sig)
}


meanrhoSHR <- function(sig, d, delta) {
  return( sig * mean( rhoSHR( d/sig )) /delta )
}

rhoSHR <- function(d)  { # SHR # Optima, squared distances (d=x^2)
  G1 <- (-1.944)
  G2 <- 1.728
  G3 <- (-0.312)
  G4 <- 0.016
  u <- (d > 9.0)
  v <- (d < 4)
  w <- (1-u)*(1-v)
  z <- v*d/2 + w*((G4/8)*(d^4) + (G3/6)*(d^3) +
                    (G2/4)*(d^2)+ (G1/2)*d + 1.792) + 3.25*u
  z <- z/3.25
  return(z)
}


weightsSHR <- function(d){ # derivative of SHR rho
  G1 <- (-1.944)
  G2 <- 1.728
  G3 <- (-0.312)
  G4 <- 0.016
  u <- (d > 9.0)
  v <- (d < 4)
  w <- (1-u)*(1-v)
  z <- v/2+ w*((G4/2)*(d^3) + (G3/2)*(d^2) +(G2/2)*d+ (G1/2));
  w <- z/3.25
  return(w)
}

consMMKur <- function(p, n) {
  # Constante para eficiewncia 90% de MM, partiendo de KSD (Pe?a-Prieto)
  # cc(1): rho "SHR" ("Optima"), cc(2): Bisquare
  x <- c(1/p, p/n, 1)
  beta <- c(4.5041, -1.1117, 0.61161, 2.5724, -0.7855, 0.716)
  beta <- matrix(beta, byrow=TRUE, ncol=3)
  return( t(x) %*% t(beta) )
}




mahdist <- function(x, center=c(0,0), cov)
{
x <- as.matrix(x)
if(any(is.na(cov)))
	stop("Missing values in covariance matrix not allowed")
if(any(is.na(center)))
	stop("Missing values in center vector not allowed")
dx <- dim(x)
p <- dx[2]
if(length(center) != p)
	stop("center is the wrong length")
# produce the inverse of the covariance matrix
if(!is.qr(cov)) cov <- qr(cov)
dc <- dim(cov$qr)
if(p != dc[1] || p != dc[2])
	stop("Covariance matrix is the wrong size")
rank.cov <- cov$rank
flag.rank <- (rank.cov < p)
mahdist <- NA
if(!flag.rank)  {
	#cov.inv <- rep(c(1, rep(0, p)), length = p * p)
	#dim(cov.inv) <- c(p, p)
	#cov.inv <- qr.coef(cov, cov.inv)
	n <- dx[1]
 	mahdist <- mahalanobis(x, center=center, cov=cov)
 	if(length(dn <- dimnames(x)[[1]]))
  		names(mahdist) <- dn
}
ans <- list(mahdist,flag.rank)
names(ans) <- c("mahdist","flag.rank")
ans
}
##############################3
desceRocke <- function(X, gamma0, muini, Vini,
                       maxsteps=5, propmin=2,
                       qs=2, maxit=50, tol=1e-4)
## Iterations to minimize Rocke's scale strting from initial vector muini and matrix Vini
{
  d <- dim(X)
  n <- d[1]
  p <- d[2]
  delta <- (1-p/n)/2 # max breakdown
  mu0 <- muini
  V0 <- Vini
  dista <- dista0 <- mahalanobis(x=X,center=mu0,cov=V0);
  gamma0 <- consRocke(p,n,'K')$gamma
  sig <- MScalRocke(x=dista, gamma=gamma0, q=qs, delta=delta) #Inicializar
  # %Buscar gama que asegure que al menos p*propmin elementos tengan w>0
  didi <- dista / sig
  dife <- sort( abs( didi - 1) )
  gg <- min( dife[ (1:n) >= (propmin*p) ] )
  gamma <- max(gg, gamma0)
  sig0 <- MScalRocke(x=dista, gamma=gamma, delta=delta, q=qs)

  iter <- 0
  difpar <- difsig <- +Inf
  while( ( ( (difsig > tol) | (difpar > tol) ) &
           (iter < maxit) ) & (difsig > 0) ) {
    iter <- iter + 1
    w <- WRoTru(tt=dista/sig, gamma=gamma, q=qs)
    mu <- colMeans( X * w ) / mean(w) # as.vector( t(w) %*% X ) / sum(w)
    Xcen <- scale(X, center=mu, scale=FALSE)
    V <- t(Xcen) %*% (Xcen * w) / n;
    V <- V / ( det(V)^(1/p) )
    dista <- mahalanobis(x=X, center=mu, cov=V)
    sig <- MScalRocke(x=dista, gamma=gamma, delta=delta, q=qs)
    # %Si no desciende, hacer Line search
    step <- 0
    delgrad <- 1
    while( (sig > sig0) & (step < maxsteps) ) {
      delgrad <- delgrad / 2
      step <- step + 1
      mu <- delgrad * mu + (1 - delgrad)*mu0
      V <- delgrad*V + (1-delgrad)*V0
      V <- V / ( det(V)^(1/p) )
      dista <- mahalanobis(x=X, center=mu, cov=V)
      sig <- MScalRocke(x=dista, gamma=gamma, delta=delta, q=qs)
    }
    dif1 <- as.vector( t(mu - mu0) %*% solve(V0, mu-mu0) ) / p
    dif2 <- max(abs(solve(V0, V)-diag(p)))
    difpar <- max(dif1, dif2)
    difsig <- 1 - sig/sig0
    mu0 <- mu
    V0 <- V
    sig0 <- sig
  }
  tmp <- scalemat(V0=V0, dis=dista, weight='M')
  V <- tmp$V
  ff <- tmp$ff
  dista <- dista/ff
  return(list(mu=mu, V=V, sig=sig, dista=dista, w=w, gamma=gamma))
}



consRocke <- function(p, n, initial) {
  if(initial=='M') {
    beta <- c(-5.4358, -0.50303, 0.4214)
  } else {
    beta <- c(-6.1357, -1.0078, 0.81564)
  }
  if( p >= 15 ) {
    a <- c(1, log(p), log(n))
    alpha <- exp( sum( beta*a ) )
    gamma <- qchisq(1-alpha, df=p)/p - 1
    gamma <- min(gamma, 1)
  } else {
    gamma <- 1
    alpha <- 1e-6
  }
  return(list(gamma=gamma, alpha=alpha))
}


WRoTru <- function(tt, gamma, q) {
  ss <- (tt - 1) / gamma
  w <- 1 - ss^q
  w[ abs(ss) > 1 ] <- 0
  return(w)
}


rhorotru <- function(tt, gamma, q) {
  u <- (tt - 1) / gamma
  y <- ( (u/(2*q)*(q+1-u^q) +0.5) )
  y[ u >= 1 ] <- 1
  y[ u < (-1) ] <- 0
  return(y)
}



MScalRocke <- function(x, gamma, q, delta = 0.5, tol=1e-5)
{
  # sigma= solucion de ave{rhorocke1(x/sigma)}=delta
  n <- length(x)
  y <- sort(abs(x))
  n1 <- floor(n * (1-delta) )
  n2 <- ceiling(n * (1 - delta) / (1 - delta/2) )
  qq <- y[c(n1, n2)]
  u <- 1 + gamma*(delta-1)/2 #asegura rho(u)<delta/2
  sigin <- c(qq[1]/(1+gamma), qq[2]/u)
  if( qq[1] >= 1) {
    tolera <- tol
  } else {
    tolera <- tol * qq[1]
  }
  if( mean(x==0) > (1-delta) ) {
    sig <- 0
  } else {
    sig <- uniroot(f=averho, interval=sigin, x=x,
                   gamma=gamma, delta=delta, q=q, tol=tolera)$root
  } # solucion de ave{rhorocke1(x/sigma)}=delta
  return(sig)
}

averho <- function(sig, x, gamma, delta, q)
  return( mean( rhorotru(x/sig, gamma, q) ) - delta )


scalemat <- function(V0, dis, weight='M')
{
  p <- dim(V0)[1]
  if( weight == 'M') {
    sig <- M_Scale(x=sqrt(dis), normz=0)^2
    cc <- 4.8421*p-2.5786 #ajuste empirico
  } else {
    sig <- median(dis)
    cc <- qchisq(0.5, df=p)
  }
  ff <- sig/cc
  return(list(ff=ff, V=V0*ff))
}


M_Scale <- function(x, normz=1, delta=0.5, tol=1e-5)
{
  n <- length(x)
  y <- sort(abs(x))
  n1 <- floor(n*(1-delta))
  n2 <- ceiling(n*(1-delta)/(1-delta/2));
  qq <- y[c(n1, n2)]
  u <- rhoinv(delta/2)
  sigin <- c(qq[1],  qq[2]/u) # intervalo inicial
  if (qq[1]>=1) {
    tolera=tol
  } else {
    tolera = tol * qq[1]
  }
  #tol. relativa o absol. si sigma> o < 1
  if( mean(x==0) >= (1-delta) ) {
    sig <- 0
  } else {
    sig <- uniroot(f=averho.uni, interval=sigin, x=x,
                   delta=delta, tol=tolera)$root
  }
  if( normz > 0) sig <- sig / 1.56
  return(sig)
}



rhobisq <- function(x) {
  r <- 1 - (1-x^2)^3
  r[ abs(x) > 1 ] <- 1
  return(r)
}

averho.uni <- function(sig, x, delta)
  return( mean( rhobisq(x/sig) ) - delta )

rhoinv <- function(x)
  return(sqrt(1-(1-x)^(1/3)))


#' Classical Covariance Estimation
#'
#' Compute an estimate of the covariance/correlation matrix and location
#' vector using classical methods.
#'
#' Its main intention is to return an object compatible to that
#' produced by \code{\link{covRob}}, but fit using classical methods.
#'
#' @param data a numeric matrix or data frame containing the data.
#' @param corr a logical flag.  If \code{corr = TRUE} then the estimated correlation matrix is computed.
#' @param center a logical flag or a numeric vector of length \code{p} (where \code{p} is the number of columns of \code{x}) specifying the center.  If \code{center = TRUE} then the center is estimated.  Otherwise the center is taken to be 0.
#' @param distance a logical flag.  If \code{distance = TRUE} the Mahalanobis distances are computed.
#' @param na.action a function to filter missing data.  The default \code{na.fail} produces an error if missing values are present.  An alternative is \code{na.omit} which deletes observations that contain one or more missing values.
#' @param unbiased a logical flag. If \code{TRUE} the unbiased estimator is returned (computed with denominator equal to \code{n-1}), else the MLE (computed with denominator equal to \code{n}) is returned.
#'
#' @return a list with class \dQuote{covClassic} containing the following elements:
#' \item{call}{an image of the call that produced the object with all the arguments named.}
#' \item{cov}{a numeric matrix containing the estimate of the covariance/correlation matrix.}
#' \item{center}{a numeric vector containing the estimate of the location vector.}
#' \item{dist}{a numeric vector containing the squared Mahalanobis distances. Only present if \code{distance = TRUE} in the \code{call}.}
#' \item{corr}{a logical flag.  If \code{corr = TRUE} then \code{cov}
#' contains an estimate of the correlation matrix of \code{x}.}
#'
#' @note Originally, and in S-PLUS, this function was called \code{cov}; it has
#' been renamed, as that did mask the function in the standard package \pkg{stats}.
#'
#'
#' @examples
#' data(wine)
#' round( covClassic(wine)$cov, 2)
#'
#' @export
covClassic <- function(data, corr = FALSE, center = TRUE, distance = TRUE,
                       na.action = na.fail, unbiased = TRUE)
{
  the.call <- match.call(expand.dots = FALSE)

  data <- na.action(data)
  if(!is.matrix(data))
    data <- as.matrix(data)

  n <- nrow(data)
  p <- ncol(data)
  dn <- dimnames(data)
  dimnames(data) <- NULL
  rowNames <- dn[[1]]
  if(is.null(rowNames)) rowNames <- 1:n
  colNames <- dn[[2]]
  if(is.null(colNames)) colNames <- paste("V", 1:p, sep = "")

  if(length(center) != p && is.logical(center))
    center <- if(center) apply(data, 2, mean) else numeric(p)

  data <- sweep(data, 2, center)

  covmat <- crossprod(data) / (if(unbiased) (n - 1) else n)

  if(distance)
    dist <- mahalanobis(data, rep(0, p), covmat)

  if(corr) {
    std <- sqrt(diag(covmat))
    covmat <- covmat / (std %o% std)
  }

  dimnames(covmat) <- list(colNames, colNames)
  names(center) <- colNames

  if(distance)
    names(dist) <- rowNames

  ans <- list(call = the.call, cov = covmat, center = center, dist = dist, corr = corr)
  oldClass(ans) <- c("covClassic")
  ans
}

# The functions below are similar to those in the Robust package

#' @export
summary.covRob <- function (object, ...)
{
    evals <- eigen(object$cov, symmetric = TRUE, only.values = TRUE)$values
    names(evals) <- paste("Eval.", 1:length(evals))
    object$evals <- evals
    object <- object[c("cov", "center", "evals", "dist")]
    oldClass(object) <- "summary.covRob"
    object
}

print.summary.covRob <- function (x, digits = max(3, getOption("digits") - 3), print.distance = TRUE,
    ...)
{
    cat("Robust Estimates of Covariance: \n")
    print(x$cov, digits = digits, ...)

    cat("\n Robust Estimates of Location: \n")
    print(x$center, digits = digits, ...)

    cat("\nEigenvalues: \n")
    print(x$evals, digits = digits, ...)

    if (print.distance && !is.null(x$dist)) {
        cat("\nSquared Mahalanobis Distances: \n")
        print(x$dist, digits = digits, ...)
    }

    invisible(x)
}

# These were taken from the robust package

#' @export
summary.covClassic <- function (object, ...)
{
    evals <- eigen(object$cov, symmetric = TRUE, only.values = TRUE)$values
    names(evals) <- paste("Eval.", 1:length(evals))
    object$evals <- evals

    object <- object[c("call", "cov", "center", "evals", "dist",
        "corr")]

    oldClass(object) <- "summary.covClassic"
    object
}

print.summary.covClassic <- function (x, digits = max(3, getOption("digits") - 3), print.distance = TRUE,
    ...)
{
    if (x$corr)
        cat("\nClassical Estimate of Correlation: \n")
    else cat("\nClassical Estimate of Covariance: \n")
    print(x$cov, digits = digits, ...)

    cat("\nClassical Estimate of Location: \n")
    print(x$center, digits = digits, ...)

    cat("\nEigenvalues: \n")
    print(x$evals, digits = digits, ...)

    if (print.distance && !is.null(x$dist)) {
        cat("\nSquared Mahalanobis Distances: \n")
        print(x$dist, digits = digits, ...)
    }

    invisible(x)
}
