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
#' @param psi a string indicating which score function to use. Valid options are "bisquare", "huber",
#' "opt" and "mopt".
#' @param eff desired asymptotic efficiency. Valid options are 0.85, 0.9 and 0.95 (default) when
#' \code{psi} = "bisquare" or "huber", and 0.85, 0.9, 0.95 (default) and 0.99 when
#' \code{psi} = "opt" or "mopt".
#' @param maxit maximum number of iterations allowed.
#' @param tol tolerance to decide convergence of the iterative algorithm.
#' @param na.rm	a logical value indicating whether \code{NA} values should be stripped before
#' the computation proceeds. Defaults to \code{FALSE}
#'
#' @return A list with the following components:
#' \item{mu}{The location estimate}
#' \item{std.mu}{Estimated standard deviation of the location estimator \code{mu}}
#' \item{disper}{M-scale/dispersion estimate}
#'
#' @author Ricardo Maronna, \email{rmaronna@retina.ar}
#'
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' set.seed(123)
#' r <- rnorm(150, sd=1.5)
#' locScaleM(r)
#' # 10% of outliers, sd of good points is 1.5
#' set.seed(123)
#' r2 <- c(rnorm(135, sd=1.5), rnorm(15, mean=-10, sd=.5))
#' locScaleM(r2)
#'
locScaleM <- MLocDis <- function(x, psi="mopt", eff=0.95, maxit=50, tol=1.e-4, na.rm = FALSE) {
  kpsi <- switch(psi, bisquare = 1, huber = 2, opt = 3, mopt = 4, 5)
  # if (psi=="bisquare") kpsi=1
  # if (psi=="huber") kpsi=2
  # } else {print(c(psi, " No such psi")); kpsi=0
  # }
  # Next 6 lines taken from mean.default()
  if (!is.numeric(x)) {
    warning("argument is not numeric: returning NA")
    return(NA_real_)
  }
  if (na.rm)
    x <- x[!is.na(x)]
  if(kpsi == 5) stop(paste0(psi, ' - No such rho function'))
  if(kpsi %in% c(1, 2)) { # Start of Ricardo's code
    kBis=c(3.44, 3.88, 4.685)
    kHub=c(0.732, 0.981, 1.34)
    kk=rbind(kBis, kHub)
    efis=c(0.85, 0.90, 0.95)
    if (is.element(eff, efis)) {keff=match(eff,efis);
    } else {print(c(eff, " No such eff")); keff=0
    }
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
    scat <- mscale(u=x-mu, delta=.5, tuning.chi=1.56, family='bisquare')
    resu <- list(mu=mu, std.mu=sigmu, disper=scat)
  } else { #start of Kjell's code
    family <- psi
    efis <- c(0.85, 0.90, 0.95, 0.99)
    keff <- match(eff,efis);
    if(!is.na(keff)) {
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
          ww <- Mwgt(resi, cc, family) #RobStatTM:::Mwgt(resi, cc, family)
          mu <- sum(ww * x) / sum(ww)
          dife <- abs(mu - mu0) / sig0
          mu0 <- mu
        }
      }
      pp <- rhoprime(resi, family, cc)
      n <- length(x)
      a <- mean(pp^2)
      b <- mean(rhoprime2(resi, family, cc))
      sigmu <- sqrt(sig0^2 * a / (n*b^2))
      f <- function(u, family, cc) {
        cc["c"] <- u
        integrate(function(x, fam, cc) rho(x, fam, cc) * dnorm(x), -Inf, Inf, fam = family, cc = cc)$value - 0.5
      }
      cc["c"] <- uniroot(f, c(0.01, 10), family = family, cc = cc)$root
      scat <- mscale(x - mu, delta = 0.5, tuning.chi = cc, family = family)
      resu <- list(mu = mu, std.mu = sigmu, disper = scat)
    } else { print(c(eff, " No such eff"))
      resu <- NA
    }
  }
  return(resu)
} # end function

wfun <- function(x,k) { #weight function
  if (k==1) ww=(1-x^2)^2 *(abs(x)<=1)
  else  ww=(abs(x)<=1)+(abs(x)>1)/(abs(x)+1.e-20)
  return(ww)
}

psif <- function(x,k) {return(x*wfun(x,k))}

psipri <- function(x,k) {
  if (k==1) pp=	(((1 - (x^2))^2) - 4 * (x^2) * (1 - (x^2))) * (abs(x) < 1)
  else pp=(abs(x)<=1)
  return(pp)
}
