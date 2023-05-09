DNORM1 <- dnorm(1.0)

psiSupportFromTuningConst <- function(a, family.name)
{
  family.name <- match.arg(family.name, choices = FAMILY.NAMES)

  switch(family.name,
    bisquare = c(0, a[1]),
    huber = c(0, +Inf), 
    mopt = psiSupportFromTuningConst_modopt(a[1]),
    opt = psiSupportFromTuningConst_optimal(a[1]),
    moptv0 = psiSupportFromTuningConst_modoptv0(a[1]),
    optv0 = psiSupportFromTuningConst_optimalv0(a[1])
  )
}


computeNuFromFamily <- function(family, cc, psiSupport = NULL, densFun = dnorm)
{
  if(is.null(psiSupport))
    psiSupport <- psiSupportFromTuningConst(cc, family)

  integrand.top <- function(x, family, cc, densFun)
    rhoprime(x, family, cc)^2 * densFun(x)

  nu.top <- integrate(integrand.top, psiSupport[1], psiSupport[2], family = family, cc = cc, densFun = densFun)$value
  # nu.top <- integrate(integrand.top, psiSupport[1], +Inf, family = family, cc = cc, densFun = densFun)$value

  integrand.bottom <- function(x, family, cc, densFun)
    rhoprime2(x, family, cc) * densFun(x)

  nu.bottom <- integrate(integrand.bottom, psiSupport[1], psiSupport[2], family = family, cc= cc, densFun = densFun)$value^2
  # nu.bottom <- integrate(integrand.bottom, psiSupport[1], +Inf, family = family, cc= cc, densFun = densFun)$value^2

  0.5 * nu.top / nu.bottom
}


computeGaussianEfficiencyFromFamily <- function(family, cc, psiSupport = NULL)
  1.0 / computeNuFromFamily(family, cc, psiSupport = psiSupport)


findTuningConstFromGaussianEfficiency <- function(e, family.name)
{
  family.name <- match.arg(family.name, choices = FAMILY.NAMES)

  switch(family.name,
    opt = findTuningConstFromGaussianEfficiency_optimal(e),
    mopt = findTuningConstFromGaussianEfficiency_modopt(e),
    optv0 = findTuningConstFromGaussianEfficiency_optimalv0(e),
    moptv0 = findTuningConstFromGaussianEfficiency_modoptv0(e),
    {
      obj <- function(a, family, e) {
        cc <-  c('c' = a)
        # names(tuning.psi) <- 'c'
        e - computeGaussianEfficiencyFromFamily(family, cc)
      }
      c('c' = uniroot(obj, interval = c(1e-2, 1e2), e = e, family=family.name, check.conv = TRUE, tol = 1e-8)$root)
    }
  )
}


########## bisquare ##########

findBreakdownPoint <- function(family, cc)
{
  f <- function(x, family, cc)
    rho(x, family, cc) * dnorm(x)

  integrate(f, lower = -Inf, upper = Inf, family = family, cc = cc)$value
}


adjustTuningVectorForBreakdownPoint <- function(family, cc, breakdown.point = 0.5)
{
  g <- function(v, family, cc, breakdown.point) {
    # family$cc["c"] <- v
    if( (family == 'opt') | (family == 'mopt') )
      cc["c2"] <- v
    else
      cc["c"] <- v
    findBreakdownPoint(family, cc=cc) - breakdown.point
  }

  # family$cc["c"] <-
  tmp <- uniroot(g, c(1e-5, 1e5), family = family, cc=cc, breakdown.point = breakdown.point, tol = 1e-8)$root
  if( (family == 'opt') | (family == 'mopt') )
    cc["c2"] <- tmp
  else
    cc["c"] <- tmp
  return(cc)
}



########## optimal ##########

Psi_optimalv0 <- function(x, a)
  0.5 * x^2 - a * pi * .Call(R_erfi, x / sqrt(2), PACKAGE='RobStatTM')

Psi_optimal <- function(x, a)
  0.5 * x^2 - a * pi * .Call(R_erfi, x / sqrt(2), PACKAGE='RobStatTM')


A_MAX <- 0.2419707

g <- function(x, a) x * dnorm(x) - a


psiSupportFromTuningConst_optimalv0 <- function(a, tol = 1e-8)
{
  if(a > 0.0 && a < A_MAX) {
    lower <- uniroot(f = g, interval = c(0.0, 1.0), a = a, check.conv = TRUE, tol = tol)

    if(lower$f.root < 0.0) {
      x <- seq(lower$root, lower$root + lower$estim.prec, length.out = 100)
      y <- g(x, a)
      i <- min(which(y >= 0.0))
      if(length(i)) {
        lower$root <- x[i]
        lower$f.root <- y[i]
      }
      else {
        warning("psi is negative at computed endpoint")
      }
    }

    upper <- uniroot(f = g, interval = c(1.0, 5.0), a = a, extendInt = "downX", check.conv = TRUE, tol = tol)

    if(upper$f.root < 0.0) {
      x <- seq(upper$root - upper$estim.prec, upper$root, length.out = 100)
      y <- g(x, a)
      i <- max(which(y >= 0.0))
      if(length(i)) {
        upper$root <- x[i]
        upper$f.root <- y[i]
      }
      else {
        warning("psi is negative at computed endpoint")
      }
    }

    return(c(lower$root, upper$root))
  }

  if(a == 0.0)
    return(c(0.0, Inf))
  
  c(NA_real_, NA_real_)
}

psiSupportFromTuningConst_optimal <- function(a, tol = 1e-8)
{
  if(a > 0.0 && a < A_MAX) {
    lower <- uniroot(f = g, interval = c(0.0, 1.0), a = a, check.conv = TRUE, tol = tol)
    
    if(lower$f.root < 0.0) {
      x <- seq(lower$root, lower$root + lower$estim.prec, length.out = 100)
      y <- g(x, a)
      i <- min(which(y >= 0.0))
      if(length(i)) {
        lower$root <- x[i]
        lower$f.root <- y[i]
      }
      else {
        warning("psi is negative at computed endpoint")
      }
    }
    
    upper <- uniroot(f = g, interval = c(1.0, 5.0), a = a, extendInt = "downX", check.conv = TRUE, tol = tol)
    
    if(upper$f.root < 0.0) {
      x <- seq(upper$root - upper$estim.prec, upper$root, length.out = 100)
      y <- g(x, a)
      i <- max(which(y >= 0.0))
      if(length(i)) {
        upper$root <- x[i]
        upper$f.root <- y[i]
      }
      else {
        warning("psi is negative at computed endpoint")
      }
    }
    
    return(c(lower$root, upper$root))
  }
  
  if(a == 0.0)
    return(c(0.0, Inf))
  
  c(NA_real_, NA_real_)
}



findTuningConstFromGaussianEfficiency_optimal <- function(e, interval = c(1e-6, 0.2255))
{
  param <- optv0(e) #findTuningConstFromGaussianEfficiency_optimalv0(e=eff)
  a <- as.vector(param[2])
  b <- as.vector(param[3])
  x <- seq(a, b, .001)
  y <- rhoprime(u=x, family='optv0', cc=param, standardize=TRUE)
  x1 <- x
  x3 <- x^3
  x5 <- x^5 
  x7 <- x^7
  d <- lm(y~ x1+x3+x5+x7, x=TRUE)
  beta <- coef(d)
  R <- rbind(c(1, a, a^3, a^5, a^7), c(1, b, b^3, b^5, b^7))
  X <- d$x
  V <- solve(t(X)%*%X)
  ccoef <- as.vector( beta - V %*% t(R) %*% solve( R %*% V %*% t(R)) %*% R %*% beta )
  u1 <- (ccoef[1]*a+(ccoef[2]*(a^2)/2) +(ccoef[3]*(a^4)/4)+(ccoef[4]*(a^6)/6)+(ccoef[5]*(a^8)/8)) 
  u2 <- (ccoef[1]*b+(ccoef[2]*(b^2)/2) +(ccoef[3]*(b^4)/4)+(ccoef[4]*(b^6)/6)+(ccoef[5]*(b^8)/8)) 
  # out <- list(coef=ccoef, a=a, b=b, u1=u1, u2=u2)
  out <- c(param, c(ccoef), a, b, 1.0, u1, u2)
  names(out) <- c(names(param), c('p1', 'p2', 'p3', 'p4', 'p5', 'lop', 'upp', 'c2', 'u1', 'u2'))
  out
  #   obj <- function(a, e) {
  #   sup <- psiSupportFromTuningConst(a, "opt")
  #   # fam <- list(name = "optimal", cc = c(a, sup, 1.0, NA, NA))
  #   fam <- 'opt'
  #   cc <- c(a, sup, 1.0, NA, NA)
  #   e - computeGaussianEfficiencyFromFamily(fam, cc, psiSupport = sup)
  # }
  # 
  # uniroot(obj, interval = interval, e = e, check.conv = TRUE, tol = 1e-8)$root
}

findTuningConstFromGaussianEfficiency_optimalv0 <- function(e, interval = c(1e-6, 0.2255))
{
  obj <- function(a, e) {
    sup <- psiSupportFromTuningConst(a, "optv0")
    # fam <- list(name = "optimal", cc = c(a, sup, 1.0, NA, NA))
    fam <- 'optv0'
    cc <- c(a, sup, 1.0, NA, NA)
    e - computeGaussianEfficiencyFromFamily(fam, cc, psiSupport = sup)
  }
  
  uniroot(obj, interval = interval, e = e, check.conv = TRUE, tol = 1e-8)$root
}


########## modopt ##########

# Same g as optimal
# g <- function(x, a) x * dnorm(x) - a

psiSupportFromTuningConst_modopt <- function(a, tol = 1e-8)
{
  if(a > 0.0 && a < A_MAX) {
    upper <- uniroot(f = g, interval = c(1.0, 5.0), a = a, extendInt = "downX", check.conv = TRUE, tol = tol)

    if(upper$f.root < 0.0) {
      x <- seq(upper$root - upper$estim.prec, upper$root, length.out = 100)
      y <- g(x, a)
      i <- max(which(y >= 0.0))
      if(length(i)) {
        upper$root <- x[i]
        upper$f.root <- y[i]
      }
      else {
        warning("psi is negative at computed endpoint")
      }
    }
    return(c(0.0, upper$root))
  }
  
  if(a == 0.0)
    return(c(0.0, Inf))
  
  c(NA_real_, NA_real_)
}

psiSupportFromTuningConst_modoptv0 <- function(a, tol = 1e-8)
{
  if(a > 0.0 && a < A_MAX) {
    upper <- uniroot(f = g, interval = c(1.0, 5.0), a = a, extendInt = "downX", check.conv = TRUE, tol = tol)
    
    if(upper$f.root < 0.0) {
      x <- seq(upper$root - upper$estim.prec, upper$root, length.out = 100)
      y <- g(x, a)
      i <- max(which(y >= 0.0))
      if(length(i)) {
        upper$root <- x[i]
        upper$f.root <- y[i]
      }
      else {
        warning("psi is negative at computed endpoint")
      }
    }
    return(c(0.0, upper$root))
  }
  
  if(a == 0.0)
    return(c(0.0, Inf))
  
  c(NA_real_, NA_real_)
}


findTuningConstFromGaussianEfficiency_modopt <- function(e, interval = c(1e-6, 0.2255))
{
  param <- moptv0(e) #findTuningConstFromGaussianEfficiency_optimalv0(e=eff)
  a <- 0.0
  b <- as.vector(param[3])
  x <- seq(a, b, .001)
  y <- rhoprime(u=x, family='moptv0', cc=param, standardize=TRUE)
  x1 <- x
  x3 <- x^3
  x5 <- x^5 
  x7 <- x^7
  d <- lm(y~ x1+x3+x5+x7, x=TRUE)
  beta <- coef(d)
  R <- rbind(c(1, a, a^3, a^5, a^7), c(1, b, b^3, b^5, b^7))
  X <- d$x
  V <- solve(t(X)%*%X)
  ccoef <- as.vector( beta - V %*% t(R) %*% solve( R %*% V %*% t(R)) %*% R %*% beta )
  u1 <- (ccoef[1]*a+(ccoef[2]*(a^2)/2) +(ccoef[3]*(a^4)/4)+(ccoef[4]*(a^6)/6)+(ccoef[5]*(a^8)/8)) 
  u2 <- (ccoef[1]*b+(ccoef[2]*(b^2)/2) +(ccoef[3]*(b^4)/4)+(ccoef[4]*(b^6)/6)+(ccoef[5]*(b^8)/8)) 
  # out <- list(coef=ccoef, a=a, b=b, u1=u1, u2=u2)
  out <- c(param, c(ccoef), a, b, 1.0, u1, u2)
  names(out) <- c(names(param), c('p1', 'p2', 'p3', 'p4', 'p5', 'lop', 'upp', 'c2', 'u1', 'u2'))
  out
}

findTuningConstFromGaussianEfficiency_modoptv0 <- function(e, interval = c(1e-6, 0.2255))
{
  obj <- function(a, e) {
    sup <- psiSupportFromTuningConst(a, "moptv0")
    # fam <- list(name = "modopt", cc = c(a, DNORM1 / (DNORM1 - a), sup[2], 1.0, NA, NA))
    fam <- 'moptv0'
    cc <- c(a, DNORM1 / (DNORM1 - a), sup[2], 1.0, NA, NA)
    e - computeGaussianEfficiencyFromFamily(fam, cc, psiSupport = sup)
  }
  
  uniroot(obj, interval = interval, e = e, check.conv = TRUE, tol = 1e-8)$root
}


#' #' @export
#' findTuningConstFromEfficiencyBisquare <- function(efficiency) {
#'   family <- 'bisquare'
#'   effi <- function(cc, family) {
#'     a <- integrate(function(a, family, cc) (rhoprime(a, family, cc)^2)*dnorm(a), cc=cc, family=family, lower=-Inf, upper=+Inf)$value
#'     b <- integrate(function(a, family, cc) rhoprime2(a, family, cc)*dnorm(a), cc=cc, family=family, lower=-Inf, upper=+Inf)$value
#'     return( 1/(a/b^2) )
#'   }
#'   return( uniroot( function(cc, eff, family) (effi(cc, family)-eff), eff=efficiency, family=family, lower=.1, upper=1e3)$root )
#' }
#'
#' #' @export
#' findTuningConstFromBDPBisquare <- function(bdp) {
#'   family <- 'bisquare'
#'   obj <- function(cc, family) {
#'     a <- integrate(function(a, family, cc) rho(a, family, cc)*dnorm(a), cc=cc, family=family, lower=-Inf, upper=+Inf)$value
#'     return( a )
#'   }
#'   return( uniroot( function(cc, bdp, family) (obj(cc, family)-bdp), bdp=bdp, family=family, lower=.1, upper=1e3)$root )
#' }


