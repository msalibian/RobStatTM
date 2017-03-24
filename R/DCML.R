
# require(lmrob)
# require(penseinit)
# require(quantreg)
 

mscale <- function(u, tol, delta=.5, max.it=100, tuning.chi=1.5477) {
  # M-scale of a sample u  
  # tol: accuracy
  # delta: breakdown point (right side)
  # Initial 
  s0 <- median(abs(u))/.6745
  err <- tol + 1
  it <- 0
  while( (err > tol) & ( it < max.it) ) {
    it <- it+1
    s1 <- sqrt( s0^2 * mean(rho((u/s0), tuning.chi)) / delta )
    err <- abs(s1-s0)/s0
    s0 <- s1
  }
  return(s0)
}



cov.dcml <- function(res.LS, res.R, CC, sig.R, t0, p, n, control) {
  ##Computation of the asymptotic covariance matrix of the DCML estimator
  t0 <- 1-t0
  rr <- res.R/sig.R
  tpr <- rhoprime(r=rr, cc=control$tuning.psi)
  c0 <- mean( tpr * res.LS )
  a1 <- mean( tpr^2 )
  b0 <- mean(rhoprime2(r=rr, cc=control$tuning.psi))
  dee <- control$bb
  if(control$corr.b) dee <- dee * (1 - p/n)
  a2 <- mscale(u=res.LS, tol=control$mscale.tol, delta=dee, tuning.chi=control$tuning.chi)
  tt <- t0^2*sig.R^2*a1/b0^2 + a2^2*(1-t0)^2 +2*t0*(1-t0)*sig.R*c0/b0
  V <- tt*solve(CC)
  return(V)
}


rho <- function(u, cc=1.5477) {
  w <- as.numeric( abs(u) <= cc )
  v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
  v <- v*6/cc^2
  return(v)
}

rhoint <- function(e)
  return(integrate(function(a, cc) rho(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value)


find.tuning.chi <- function(delta, low=.5, upp=10) {
  return( uniroot( function(e) (rhoint(e)-delta), lower=low, upper=upp)$root )
}


rhoprime <- function(r, cc) { 
  # bisquare rho' = psi function
  r <- r / cc
  gg <- r*(1-r^2)^2*as.numeric(abs(r)<=1)
  return(gg)
}



rhoprime2 <- function(r, cc) {
  #Derivative of the bisquare psi function
  r <- r/cc
  gg <- (1-r^2)*(1-5*r^2)*as.numeric(abs(r)<=1)/cc
  return(gg)
}

# effi <- function(e) {
#   a <- integrate(function(a, cc) (rhoprime(a, cc)^2)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
#   b <- integrate(function(a, cc) rhoprime2(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
#   return( 1/(a/b^2) ) # efficiency!
# }
# 
# uniroot( function(e) (effi(e)-.85), lower=3, upper=4)$root


        
MMPY <- function(X, y, control, mf, corr.b=control$corr.b) {
   # This function will be called from lmrob, so control will be valid
   # X will already contain a column of ones if needed
   # Compute an MM-estimator taking as initial Pe?a Yohai
   # INPUT
   # X nxp matrix, where n is the number of observations and p the number of  columns
   # y vector of dimension  n with the responses
   #
   # OUTPUT
   # outMM output of the MM estimator (lmrob) with 85% of efficiency and PY as initial
   n <- nrow(X)
   p <- ncol(X)
   dee <- control$bb
   if(corr.b) dee <- dee * (1-(p/n))
   a <- pyinit(X=X, y=y, intercept=FALSE, deltaesc=dee, 
               cc.scale=control$tuning.chi, 
               prosac=control$prosac, clean.method=control$clean.method, 
               C.res = control$C.res, prop=control$prop, 
               py.nit = control$py.nit, en.tol=control$en.tol, 
               mscale.maxit = control$mscale.maxit, mscale.tol = control$mscale.tol,
               mscale.rho.fun=control$mscale.rho.fun)
   # refine the PY candidates to get something closer to an S-estimator for y ~ X1
   kk <- dim(a$initCoef)[2]
   best.ss <- +Inf
   for(i in 1:kk) {
     tmp <- refine.sm(x=X, y=y, initial.beta=a$initCoef[,i], 
                      initial.scale=a$objF[i], 
                      k=control$refine.PY, conv=1, b=dee, cc=control$tuning.chi, step='S')
     if(tmp$scale.rw < best.ss) {
       best.ss <- tmp$scale.rw # initial$objF[1]
       betapy <- tmp$beta.rw # initial$initCoef[,1]
     }
   }
   S.init <- list(coef=betapy, scale=best.ss)
   control$method <- 'M'
   control$cov <- ".vcov.w"
   control$subsampling <- 'simple' 
   # # lmrob() does the above when is.list(init)==TRUE, in particular:
   outMM <- lmrob.fit(X, y, control, init=S.init, mf=mf)
   return(outMM)
}



DCML <- function(x, y, z, z0, control) {
  # x: design matrix
  # z: robust fit
  # z0: LS fit
  beta.R <- as.vector(z$coefficients)
  beta.LS <- as.vector(z0$coefficients)
  p <- length(beta.R)
  n <- length(z$residuals)
  dee <- control$bb
  if(control$corr.b) dee <- dee*(1-p/n)
  si.dcml <- mscale(u=z$residuals, tol = control$mscale.tol, delta=dee, tuning.chi=control$tuning.chi)
  deltas <- .3*p/n 
  CC <- t(x * z$rweights) %*% x / sum(z$rweights) 	
  # print(all.equal(CC, t(x) %*% diag(z$rweights) %*% x / sum(z$rweights)))
  d <- as.numeric(crossprod(beta.R-beta.LS, CC %*% (beta.R-beta.LS)))/si.dcml^2
  t0 <- min(1, sqrt(deltas/d))
  beta.dcml <- t0*coef(z0) +(1-t0)*coef(z)
  V.dcml <- cov.dcml(res.LS=z0$residuals, res.R=z$residuals, CC=CC, 
                     sig.R=si.dcml, t0=t0, p=p, n=n, control=control) / n
  re.dcml <- as.vector(y - x %*% beta.dcml)
  si.dcml.final <- mscale(u=re.dcml, tol = control$mscale.tol, delta=dee, tuning.chi=control$tuning.chi)
  return(list(coefficients=beta.dcml, cov=V.dcml, residuals=re.dcml, scale=si.dcml.final, t0=t0))
}


SMPY <- function(mf, y, control, split, corr.b=control$corr.b) { 
  if(missing(control)) 
    control <- lmrob2.control(tuning.chi = 1.5477, bb = 0.5, tuning.psi = 3.4434)
  # int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1)
  if(missing(split)) {
    split <- splitFrame(mf, type=control$split.type) 
  }
  # step 1 - build design matrices, x1 = factors, x2 = continuous, intercept is in x1
  Z <- split$x1  
  X <- split$x2
  n <- nrow(X)
  q <- ncol(Z)
  p <- ncol(X)
  dee <- control$bb
  gamma <- matrix(NA, q, p)
  # Eliminate Z from ea. column of X with L1 regression, result goes in X1
  for( i in 1:p) gamma[,i] <- coef( lmrob.lar(x=Z, y=X[,i], control = control, mf = NULL) ) # coef(rq(X[,i]~Z-1))
  X1 <- X - Z %*% gamma
  # Eliminate Z from y, L1 regression, result goes in y1
  tmp0 <- lmrob.lar(x=Z, y=y, control = control, mf = NULL)
  y1 <- as.vector(tmp0$residuals)
  # Now regress y1 on X1, find PY candidates
  if(corr.b) dee <- dee*(1-p/n)
  initial <- pyinit(intercept=FALSE, X=X1, y=y1, 
                    deltaesc=dee, cc.scale=control$tuning.chi, 
                    prosac=control$prosac, clean.method=control$clean.method, 
                    C.res = control$C.res, prop=control$prop, 
                    py.nit = control$py.nit, en.tol=control$en.tol, 
                    mscale.maxit = control$mscale.maxit, mscale.tol = control$mscale.tol,
                    mscale.rho.fun=control$mscale.rho.fun)
  # choose best candidates including factors into consideration!
  # recompute scales adjusting for Z
  dee <- control$bb
  if(corr.b) dee <- dee*(1-(ncol(X1)+ncol(Z))/n)
  kk <- dim(initial$initCoef)[2]
  # do the first, then iterate over the rest looking for a better one
  betapy <- initial$initCoef[,1]
  r <- as.vector(y1 - X1 %*% betapy)
  best.tmp <- lmrob.lar(x=Z, y=r, control = control, mf = NULL)
  sspy <- mscale(u=best.tmp$residuals, tol=control$mscale.tol, delta=dee, tuning.chi=control$tuning.chi)
  for(i in 2:kk) {
    r <- as.vector(y1 - X1 %*% initial$initCoef[,i])
    tmp <- lmrob.lar(x=Z, y=r, control = control, mf = NULL)
    s.cand <- mscale(u=tmp$residuals, tol=control$mscale.tol, delta=dee, tuning.chi=control$tuning.chi)
    if( s.cand < sspy ) {
      sspy <- s.cand
      betapy <- initial$initCoef[,i]
      best.tmp <- tmp
    }
  }
  tmp <- best.tmp
  # re-express estimators in terms of X, Z
  gammapy <- as.vector(coef(tmp0) - gamma %*% betapy)
  gammapy <- gammapy + tmp$coef # gamma.cand
  res <- tmp$residuals
  sih <- sspy
  # Run a few IRWLS iterations, adjusting with Z
  max.it <- control$refine.PY
  for(i in 1:max.it) {
    weights <- f.w( tmp$residuals/sih, cc=control$tuning.chi)
    xw <- X * sqrt(weights)
    yw <- y * sqrt(weights)
    beta <- our.solve( t(xw) %*% xw ,t(xw) %*% yw )
    r1 <- as.vector(y - X %*% beta)
    tmp <- lmrob.lar(x=Z, y=r1, control = control, mf = NULL)
    sih <- mscale(u=tmp$residuals, tol=control$mscale.tol, delta=dee, tuning.chi=control$tuning.chi)
    if(sih < sspy) {
      sspy <- sih
      betapy <- beta
      gammapy <- tmp$coeff
      res <- tmp$residuals
    }
  }  
  beta00 <- c(betapy, gammapy)
  ss <- sspy
  XX <- model.matrix(attr(mf, 'terms'), mf)
  ii <- charmatch('(Intercept)', colnames(cbind(X, Z)))
  if(!is.na(ii)) beta00 <- c(beta00[ii], beta00[-ii])
  uu <- list(coef=beta00, scale=ss, residuals=res)
  # Compute the MMestimator using lmrob, starting from this initial
  # and associated residual scale estimate
  control$method <- 'M' 
  control$cov <- ".vcov.w"
  control$subsampling <- 'simple' 
  # lmrob() sets the above when is.list(init)==TRUE
  outlmrob <- lmrob.fit(XX, y, control, init=uu, mf=mf)
  return(outlmrob) #, init.SMPY=uu)) 
}
