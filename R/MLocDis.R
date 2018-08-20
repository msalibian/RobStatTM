#' Robust univariate location and scale M-estimators
#'
#' This function computes M-estimators for location and scale.
#'
#' This function computes M-estimators for location and scale.
#'
#' @export MLocDis locScaleM
#' @aliases MLocDis locScaleM
#' @rdname MLocDis
#'
#' @param x a vector of univariate observations
#' @param psi a string indicating which score function to use. Valid options are "Bis" for
#' bi-square and "Hub" for a Huber-type.
#' @param eff desired asymptotic efficiency. Valid options are 0.9 (default), 0.85 and 0.95.
#' @param maxit maximum number of iterations allowed.
#' @param tol tolerance to decide convergence of the iterative algorithm.
#'
#' @return A list with the following components:
#' \item{mu}{The location estimator}
#' \item{std.mu}{Estimated standard deviation of the location estimator \code{mu}}
#' \item{disper}{M-scale/dispersion estimator}
#'
#' @author Ricardo Maronna, \email{rmaronna@retina.ar}
#'
#' @references \url{http://thebook}
#'
locScaleM <- MLocDis <- function(x, family = "bisquare", eff = 0.9, maxit = 50, tol = 1.e-4)) {
  cc <- get(family)(eff)

  mu0 <- median(x)
  sig0 <- mad(x)

  if(sig0 < 1.e-10) {
    mu <- 0.0
    sigma <- 0.0
  }
  
  else {
    dife <- 1.e10
    iter <- 0

    while(dife > tol && iter < maxit) {
      iter <- iter + 1
      resi <- (x - mu0) / sig0
      ww <- RobStatTM:::Mwgt(resi, cc, family)
      mu <- sum(ww * x) / sum(ww)
      dife <- abs(mu - mu0) / sig0
      mu0 <- mu
    }
  }

  pp <- rhoprime(resi, family, cc)
  n <- length(x)
  a <- mean(pp^2)
  b <- mean(rhoprime2(resi, family, cc))
  sigmu <- sqrt(sig0^2 * a / (n * b^2))
  
  f <- function(u, family, cc) {
    cc["c"] <- u
    integrate(function(x, fam, cc) rho(x, family, cc) * dnorm(x), -Inf, Inf, fam = famly, cc = cc)$value - 0.5
  }
  
  cc["c"] <- uniroot(f, c(0.01, 10), family = family, cc = cc)$root
  scat <- mscale(x - mu, delta = 0.5, tuning.chi = cc, family = family)

  list(mu = mu, std.mu = sigmu, disper = scat)
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
