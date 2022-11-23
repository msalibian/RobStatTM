#' M-scale estimator
#'
#' This function computes an M-scale, which is a robust
#' scale (spread) estimator.
#' M-estimators of scale are a robust alternative to
#' the sample standard deviation. Given a vector of
#' residuals \code{r}, the M-scale estimator \code{s}
#' solves the non-linear equation \code{mean(rho(r/s, cc))=b},
#' where \code{b} and \code{cc} are user-chosen tuning constants.
#' In this package the function \code{rho} is one of
#' Tukey's bisquare family.
#' The breakdown point of the estimator is \code{min(b, 1-b)},
#' so the optimal choice for \code{b} is 0.5. To obtain a
#' consistent estimator the constant
#' \code{cc} should be chosen such that E(rho(Z, cc)) = b, where
#' Z is a standard normal random variable.
#'
#' The iterative algorithm starts from the scaled median of
#' the absolute values of the input vector, and then
#' cycles through the equation s^2 = s^2 * mean(rho(r/s, cc)) / b.
#'
#' @export mscale scaleM
#' @aliases mscale scaleM
#' @rdname mscale
#'
#' @param u vector of residuals
#' @param delta the right hand side of the M-scale equation
#' @param family string specifying the name of the family of loss function to be used (current valid
#' options are "bisquare", "opt" and "mopt").
#' @param tuning.chi the tuning object for the rho function as returned
#' by \code{\link{lmrobdet.control}}, \link{bisquare}, \link{mopt} or \link{opt}.
#' It should correspond to the family of rho functions specified in the argument \code{family}.
#' @param tol relative tolerance for convergence
#' @param max.it maximum number of iterations allowed
#' @param tolerancezero smallest (in absolute value) non-zero value accepted as a scale. Defaults to \code{.Machine$double.eps}
#'
#' @return The scale estimate value at the last iteration or at convergence.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @examples
#' set.seed(123)
#' r <- rnorm(150, sd=1.5)
#' mscale(r)
#' sd(r)
#' # 10% of outliers, sd of good points is 1.5
#' set.seed(123)
#' r2 <- c(rnorm(135, sd=1.5), rnorm(15, mean=-5, sd=.5))
#' mscale(r2)
#' sd(r2)
#'
scaleM <- mscale <- function(u, delta=0.5, tuning.chi=1.547645, family ="bisquare",
                             max.it=100, tol=1e-6, tolerancezero=.Machine$double.eps) {
  # M-scale of a sample u
  # tol: accuracy
  # delta: breakdown point (right side)
  # Initial
  s0 <- median(abs(u))/.6745
  if(s0 < tolerancezero) return(0)
  err <- tol + 1
  it <- 0
  while( (err > tol) && ( it < max.it) ) {
    it <- it+1
    s1 <- sqrt( s0^2 * mean(rho(u/s0, family = family, cc = tuning.chi)) / delta )
    err <- abs(s1-s0)/s0
    s0 <- s1
  }
  return(s0)
}


#' Approximate covariance matrix of the DCML regression estimator.
#'
#' The estimated covariance matrix of the DCML regression estimator.
#' This function is used internally and not meant to be used
#' directly.
#'
#' @param res.LS vector of residuals from the least squares fit
#' @param res.R vector of residuals from the robust regression fit
#' @param CC estimated covariance matrix of the robust regression estimator
#' @param sig.R robust estimate of the scale of the residuals
#' @param t0 mixing parameter
#' @param p,n the dimensions of the problem, needed for the finite
#' sample correction of the tuning constant of the M-scale
#' @param control a list of control parameters as returned by \code{\link{lmrobdet.control}}
#'
#' @return The covariance matrix estimate.
#'
#' @rdname cov.dcml
#' @author Victor Yohai, \email{victoryohai@gmail.com}
#'
#' @export
cov.dcml <- function(res.LS, res.R, CC, sig.R, t0, p, n, control) {
  ##Computation of the asymptotic covariance matrix of the DCML estimator
  t0 <- 1-t0
  rr <- res.R/sig.R
  #tpr <- rhoprime(r=rr, cc=control$tuning.psi)
  tpr <- rhoprime(rr, family = control$family, cc = control$tuning.psi)
  c0 <- mean( tpr * res.LS )
  a1 <- mean( tpr^2 )
  #b0 <- mean(rhoprime2(r=rr, cc=control$tuning.psi))
  b0 <- mean(rhoprime2(rr, family = control$family, cc = control$tuning.psi))
  dee <- control$bb
  if(control$corr.b) dee <- dee * (1 - p/n)
  a2 <- mscale(u=res.LS, tol=control$mscale_tol, delta=dee, tuning.chi=control$tuning.chi, family=control$family, max.it = control$mscale_maxit)
  tt <- t0^2*sig.R^2*a1/b0^2 + a2^2*(1-t0)^2 +2*t0*(1-t0)*sig.R*c0/b0
  V <- tt*solve(CC)
  return(V)
}

#' MM regression estimator using Pen~a-Yohai candidates
#'
#' This function computes MM-regression estimator using Pen~a-Yohai
#' candidates for the initial S-estimator. This function is used
#' internally by \code{\link{lmrobdetMM}}, and not meant to be used
#' directly.
#'
#' @param X design matrix
#' @param y response vector
#' @param control a list of control parameters as returned by \code{\link{lmrobdet.control}}
#' @param mf model frame
#'
#' @return an \code{\link{lmrob}} object witht the M-estimator
#' obtained starting from the S-estimator computed with the
#' Pen~a-Yohai initial candidates. The properties of the final
#' estimator (efficiency, etc.) are determined by the tuning constants in
#' the argument \code{control}.
#'
#' @rdname MMPY
#' @author Victor Yohai, \email{victoryohai@gmail.com}, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \code{\link{DCML}}, \code{\link{MMPY}}, \code{\link{SMPY}}
#'
#' @export
MMPY <- function(X, y, control, mf) {
  # INPUT
  # X nxp matrix, where n is the number of observations and p the number of  columns
  # y vector of dimension  n with the responses
  #
  # OUTPUT
  # outMM output of the MM estimator (lmrob) with 85% of efficiency and PY as initial
  n <- nrow(X)
  p <- ncol(X)
  dee <- control$bb
  if(control$corr.b) dee <- dee * (1-(p/n))
  a <- pyinit::pyinit(x=X, y=y, intercept=FALSE, delta=dee,
                      cc=1.54764,
                      psc_keep=control$psc_keep*(1-(p/n)), resid_keep_method=control$resid_keep_method,
                      resid_keep_thresh = control$resid_keep_thresh, resid_keep_prop=control$resid_keep_prop,
                      maxit = control$py_maxit, eps=control$py_eps,
                      mscale_maxit = control$mscale_maxit, mscale_tol = control$mscale_tol,
                      mscale_rho_fun="bisquare")
  # refine the PY candidates to get something closer to an S-estimator for y ~ X1
  kk <- dim(a$coefficients)[2]
  best.ss <- +Inf
  iters.py <- vector('numeric', kk)
  for(i in 1:kk) {
    tmp <- refine.sm(x=X, y=y, initial.beta=a$coefficients[,i], initial.scale=a$objective[i],
                     k=control$refine.PY, conv=1, b=dee, family=control$family, cc=control$tuning.chi, step='S',
                     tol=control$refine.S.py)
    iters.py[i] <- tmp$iterations
    if(tmp$scale.rw < best.ss) {
      best.ss <- tmp$scale.rw # initial$objF[1]
      betapy <- tmp$beta.rw # initial$initCoef[,1]
    }
  }
  S.init <- list(coef=betapy, scale=best.ss)
  orig.control <- control

  control$psi <- control$family # tuning.psi$name
  # control$tuning.psi <- control$tuning.psi # $cc

  control$method <- 'M'
  control$cov <- ".vcov.w"
  control$subsampling <- 'simple'

  # # lmrob() does the above when is.list(init)==TRUE, in particular:
  if(S.init$scale > 0) {
    outMM <- lmrob.fit(X, y, control, init=S.init, mf=mf)
    outMM$control <- orig.control
    coefnames <- names(coef(outMM))
    residnames <- names(resid(outMM))
    outMM$iters.py <- iters.py
    # if weights were all zero, then return the S estimator
    # S.resid <- as.vector(y - X %*% S.init$coef) / S.init$scale
    # ws <- rhoprime(S.resid, family=orig.control$family, cc=orig.control$tuning.psi)
  } else {
    outMM <- list(coefficients = S.init$coef, scale = S.init$scale)
    outMM$fitted.values <- as.vector(X %*% S.init$coef) # y - b$y # y = fitted + res
    outMM$residuals <- setNames(y - outMM$fitted.values, rownames(X))
    names(outMM$coefficients) <- colnames(X)
    outMM$iters.py <- iters.py
    ## robustness weights
    outMM$rweights <- lmrob.rweights(outMM$residuals, 
                                     outMM$scale, 
                                     control$tuning.chi, 
                                     control$psi)
    outMM$converged <- FALSE
  }
  if(all( outMM$rweights == 0 ) ) {
    outMM$coefficients <- as.vector(S.init$coef)
    outMM$scale <- S.init$scale
    names(outMM$coefficients) <- coefnames
    outMM$residuals <- as.vector(y - X %*% S.init$coef)
    names(outMM$residuals) <- residnames
  }
  return(outMM)
}


#' DCML regression estimator
#'
#' This function computes the DCML regression estimator. This function is used
#' internally by \code{\link{lmrobdetDCML}}, and not meant to be used
#' directly.
#'
#' @param x design matrix
#' @param y response vector
#' @param z robust fit as returned by \code{\link{MMPY}} or \code{\link{SMPY}}
#' @param z0 least squares fit as returned by \code{\link{lm.fit}}
#' @param control a list of control parameters as returned by \code{\link{lmrobdet.control}}
#'
#' @return a list with the following components
#' \item{coefficients}{the vector of regression coefficients}
#' \item{cov}{the estimated covariance matrix of the DCML regression estimator}
#' \item{residuals}{the vector of regression residuals from the DCML fit}
#' \item{scale}{a robust residual (M-)scale estimate}
#' \item{t0}{the mixing proportion between the least squares and robust regression estimators}
#'
#' @rdname DCML
#' @author Victor Yohai, \email{victoryohai@gmail.com}, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \code{\link{DCML}}, \code{\link{MMPY}}, \code{\link{SMPY}}
#'
#' @export
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
  si.dcml <- mscale(u=z$residuals, tol = control$mscale_tol, delta=dee, tuning.chi=control$tuning.chi, family=control$family, max.it = control$mscale_maxit)
  deltas <- .3*p/n
  CC <- t(x * z$rweights) %*% x / sum(z$rweights)
  # print(all.equal(CC, t(x) %*% diag(z$rweights) %*% x / sum(z$rweights)))
  d <- as.numeric(crossprod(beta.R-beta.LS, CC %*% (beta.R-beta.LS)))/si.dcml^2
  t0 <- min(1, sqrt(deltas/d))
  beta.dcml <- t0*coef(z0) +(1-t0)*coef(z)
  V.dcml <- cov.dcml(res.LS=z0$residuals, res.R=z$residuals, CC=CC,
                     sig.R=si.dcml, t0=t0, p=p, n=n, control=control) / n
  fi.dcml <- as.vector( x %*% beta.dcml )
  re.dcml <- as.vector(y - fi.dcml )
  si.dcml.final <- mscale(u=re.dcml, tol = control$mscale_tol, delta=dee, tuning.chi=control$tuning.chi, family=control$family, max.it = control$mscale_maxit)
  return(list(coefficients=beta.dcml, cov=V.dcml, residuals=re.dcml, fitted.values = fi.dcml,
              scale=si.dcml.final, t0=t0))
}

#' SM regression estimator using Pen~a-Yohai candidates
#'
#' This function computes a robust regression estimator when there
#' are categorical / dummy explanatory variables. It uses Pen~a-Yohai
#' candidates for the S-estimator. This function is used
#' internally by \code{\link{lmrobdetMM}}, and not meant to be used
#' directly.
#'
#' @param mf model frame
#' @param y response vector
#' @param control a list of control parameters as returned by \code{\link{lmrobdet.control}}
#' @param split a list as returned by \code{\link{splitFrame}} containing the continuous and
#' dummy components of the design matrix
#'
#' @return an \code{\link{lmrob}} object witht the M-estimator
#' obtained starting from the MS-estimator computed with the
#' Pen~a-Yohai initial candidates. The properties of the final
#' estimator (efficiency, etc.) are determined by the tuning constants in
#' the argument \code{control}.
#'
#' @rdname SMPY
#' @author Victor Yohai, \email{victoryohai@gmail.com}, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \code{\link{DCML}}, \code{\link{MMPY}}, \code{\link{SMPY}}
#'
#' @export
SMPY <- function(mf, y, control, split) {
  if(missing(control))
    control <- lmrobdet.control(tuning.chi = 1.5477, bb = 0.5, tuning.psi = 3.4434)
  # int.present <- (attr(attr(mf, 'terms'), 'intercept') == 1)
  if(missing(split)) {
    split <- splitFrame(mf, type=control$split.type)
  }
  f.w <- function(u, family, cc)
    Mwgt(x = u, cc = cc, psi = family)
  # step 1 - build design matrices, x1 = factors, x2 = continuous, intercept is in x1
  Z <- split$x1
  X <- split$x2
  n <- nrow(X)
  q <- ncol(Z)
  p <- ncol(X)
  # if( p == 0 ) there are only factors (no continuous), just use L1
  dee <- control$bb
  gamma <- matrix(NA, q, p)
  # Eliminate Z from ea. column of X with L1 regression, result goes in X1
  for( i in 1:p) gamma[,i] <- coef( lmrob.lar(x=Z, y=X[,i], control = control, mf = NULL) ) # coef(rq(X[,i]~Z-1))
  X1 <- X - Z %*% gamma
  # Eliminate Z from y, L1 regression, result goes in y1
  tmp0 <- lmrob.lar(x=Z, y=y, control = control, mf = NULL)
  y1 <- as.vector(tmp0$residuals)
  # Now regress y1 on X1, find PY candidates
  if(control$corr.b) dee <- dee*(1-p/n)
  initial <- pyinit::pyinit(intercept=FALSE, x=X1, y=y1,
                            delta=dee, cc=1.54764,
                            psc_keep=control$psc_keep*(1-(p/n)), resid_keep_method=control$resid_keep_method,
                            resid_keep_thresh = control$resid_keep_thresh, resid_keep_prop=control$resid_keep_prop,
                            maxit = control$py_maxit, eps=control$py_eps,
                            mscale_maxit = control$mscale_maxit, mscale_tol = control$mscale_tol,
                            mscale_rho_fun="bisquare")
  # choose best candidates including factors into consideration!
  # recompute scales adjusting for Z
  dee <- control$bb
  if(control$corr.b) dee <- dee*(1-(ncol(X1)+ncol(Z))/n)
  kk <- dim(initial$coefficients)[2]
  # do the first, then iterate over the rest looking for a better one
  betapy <- initial$coefficients[,1]
  r <- as.vector(y1 - X1 %*% betapy)
  best.tmp <- lmrob.lar(x=Z, y=r, control = control, mf = NULL)
  sspy <- mscale(u=best.tmp$residuals, tol=control$mscale_tol, delta=dee, tuning.chi=control$tuning.chi, family=control$family, max.it = control$mscale_maxit)
  for(i in 2:kk) {
    r <- as.vector(y1 - X1 %*% initial$coefficients[,i])
    tmp <- lmrob.lar(x=Z, y=r, control = control, mf = NULL)
    s.cand <- mscale(u=tmp$residuals, tol=control$mscale_tol, delta=dee, tuning.chi=control$tuning.chi, family=control$family, max.it = control$mscale_maxit)
    if( s.cand < sspy ) {
      sspy <- s.cand
      betapy <- initial$coefficients[,i]
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
    weights <- f.w( tmp$residuals/sih, family=control$family, control$tuning.chi)
    xw <- X * sqrt(weights)
    yw <- y * sqrt(weights)
    beta <- our.solve( t(xw) %*% xw ,t(xw) %*% yw )
    r1 <- as.vector(y - X %*% beta)
    tmp <- lmrob.lar(x=Z, y=r1, control = control, mf = NULL)
    sih <- mscale(u=tmp$residuals, tol=control$mscale_tol, delta=dee, tuning.chi=control$tuning.chi, family=control$family, max.it = control$mscale_maxit)
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
  orig.control <- control
  control$method <- 'M'
  control$cov <- ".vcov.w"
  control$subsampling <- 'simple'
  # lmrob() sets the above when is.list(init)==TRUE
  # control$psi <- control$tuning.psi$name
  # control$tuning.psi <- control$tuning.psi$cc
  control$psi <- control$family # tuning.psi$name
  outlmrob <- lmrob.fit(XX, y, control, init=uu, mf=mf)
  outlmrob$control <- orig.control
  return(outlmrob) #, init.SMPY=uu))
}


