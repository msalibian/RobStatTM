DNORM1 <- dnorm(1.0)

psiSupportFromTuningConst <- function(a, family.name)
{
  family.name <- match.arg(family.name, choices = FAMILY.NAMES)
  
  switch(family.name,
    bisquare = c(0, a[1]),
    modopt = psiSupportFromTuningConst_modopt(a[1]),
    optimal = psiSupportFromTuningConst_optimal(a[1])
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
    optimal = findTuningConstFromGaussianEfficiency_optimal(e),
    modopt = findTuningConstFromGaussianEfficiency_modopt(e),
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
    cc["c"] <- v
    findBreakdownPoint(family, cc=cc) - breakdown.point
  }

  # family$cc["c"] <-
  cc["c"] <- uniroot(g, c(1e-5, 1e5), family = family, cc=cc, breakdown.point = breakdown.point, tol = 1e-8)$root
  return(cc)
}



########## optimal ##########

Psi_optimal <- function(x, a)
  0.5 * x^2 - a * pi * .Call(R_erfi, x / sqrt(2), PACKAGE='RobStatTM')


psiSupportFromTuningConst_optimal <- function(a)
{
  u <- 0.77957 - 0.33349 * log(a)
  # fam <- list(name = "optimal", cc = c(a, 1e-9, Inf, 1.0, NA, NA))
  fam <- 'optimal'
  cc <- c(a, 1e-9, Inf, 1.0, NA, NA)
  lower <- uniroot(f = rhoprime, interval = c(1e-6, u), family = fam, cc = cc, check.conv = TRUE, tol = 1e-8)$root
  upper <- uniroot(f = rhoprime, interval = c(u, 1.5*u), family = fam, cc = cc, check.conv = TRUE, tol = 1e-8)$root
  c(lower, upper)
}


findTuningConstFromGaussianEfficiency_optimal <- function(e, interval = c(1e-6, 0.2255))
{
  obj <- function(a, e) {
    sup <- psiSupportFromTuningConst(a, "optimal")
    # fam <- list(name = "optimal", cc = c(a, sup, 1.0, NA, NA))
    fam <- 'optimal'
    cc <- c(a, sup, 1.0, NA, NA)
    e - computeGaussianEfficiencyFromFamily(fam, cc, psiSupport = sup)
  }

  uniroot(obj, interval = interval, e = e, check.conv = TRUE, tol = 1e-8)$root
}

########## modopt ##########


psiSupportFromTuningConst_modopt <- function(a)
{
  u <- 0.77957 - 0.33349 * log(a)
  # fam <- list(name = "modopt", cc = c(a, 1e-9, Inf, 1.0, 0.0, 0.0))
  fam <- 'modopt'
  cc <- c(a, 1e-9, Inf, 1.0, 0.0, 0.0)
  upper <- uniroot(f = rhoprime, interval = c(u, 1.5*u), family = fam, cc = cc, check.conv = TRUE, tol = 1e-8)$root
  c(0.0, upper)
}


findTuningConstFromGaussianEfficiency_modopt <- function(e, interval = c(1e-6, 0.2255))
{
  obj <- function(a, e) {
    sup <- psiSupportFromTuningConst(a, "modopt")
    # fam <- list(name = "modopt", cc = c(a, DNORM1 / (DNORM1 - a), sup[2], 1.0, NA, NA))
    fam <- 'modopt'
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


