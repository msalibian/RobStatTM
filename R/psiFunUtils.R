DNORM1 <- dnorm(1.0)

psiSupportFromTuningConst <- function(a, family.name)
{
  family.name <- match.arg(family.name, choices = FAMILY.NAMES)
  
  switch(family.name,
    bisquare = c(0, a[1]),
    modified.optimal = psiSupportFromTuningConst_modified.optimal(a[1]),
    optimal = psiSupportFromTuningConst_optimal(a[1])
  )
}


computeNuFromFamily <- function(family, psiSupport = NULL, densFun = dnorm)
{
  if(is.null(psiSupport))
    psiSupport <- psiSupportFromTuningConst(family$cc, family$name)

  integrand.top <- function(x, family, densFun)
    rhoprime(x, family)^2 * densFun(x)
  
  nu.top <- integrate(integrand.top, psiSupport[1], psiSupport[2], family = family, densFun = densFun)$value
  
  integrand.bottom <- function(x, family, densFun)
    rhoprime2(x, family) * densFun(x)
  
  nu.bottom <- integrate(integrand.bottom, psiSupport[1], psiSupport[2], family = family, densFun = densFun)$value^2
  
  0.5 * nu.top / nu.bottom
}


computeGaussianEfficiencyFromFamily <- function(family, psiSupport = NULL)
  1.0 / computeNuFromFamily(family, psiSupport = psiSupport)


findTuningConstFromGaussianEfficiency <- function(e, family.name)
{
  family.name <- match.arg(family.name, choices = FAMILY.NAMES)

  switch(family.name,
    optimal = findTuningConstFromGaussianEfficiency_optimal(e),
    modified.optimal = findTuningConstFromGaussianEfficiency_modified.optimal(e),
    {
      obj <- function(a, e) {
        fam <- list(name = family.name, cc = a)
        e - computeGaussianEfficiencyFromFamily(fam)
      }
      
      uniroot(obj, interval = c(0.1, 15.0), e = e, check.conv = TRUE, tol = 1e-8)$root
    }
  )
}


########## bisquare ##########

findBreakdownPoint <- function(family)
{
  f <- function(x, family)
    rho(x, family = family) * dnorm(x)

  integrate(f, lower = -Inf, upper = Inf, family = family)$value
}


adjustTuningVectorForBreakdownPoint <- function(family, breakdown.point = 0.5)
{
  g <- function(v, family, breakdown.point) {
    family$cc["c"] <- v
    findBreakdownPoint(family) - breakdown.point
  }

  family$cc["c"] <- uniroot(g, c(0.1, 15.0), family = family, breakdown.point = breakdown.point, tol = 1e-8)$root
  family
}



########## optimal ##########

Psi_optimal <- function(x, a)
  0.5 * x^2 - a * pi * .Call(R_erfi, x / sqrt(2))


psiSupportFromTuningConst_optimal <- function(a)
{
  u <- 0.77957 - 0.33349 * log(a)
  fam <- list(name = "optimal", cc = c(a, 1e-9, Inf, 1.0, NA, NA))
  lower <- uniroot(f = rhoprime, interval = c(1e-6, u), family = fam, check.conv = TRUE, tol = 1e-8)$root
  upper <- uniroot(f = rhoprime, interval = c(u, 1.5*u), family = fam, check.conv = TRUE, tol = 1e-8)$root
  c(lower, upper)
}


findTuningConstFromGaussianEfficiency_optimal <- function(e, interval = c(1e-6, 0.2255))
{
  obj <- function(a, e) {
    sup <- psiSupportFromTuningConst(a, "optimal")
    fam <- list(name = "optimal", cc = c(a, sup, 1.0, NA, NA))
    e - computeGaussianEfficiencyFromFamily(fam, psiSupport = sup)
  }

  uniroot(obj, interval = interval, e = e, check.conv = TRUE, tol = 1e-8)$root
}

########## modified.optimal ##########


psiSupportFromTuningConst_modified.optimal <- function(a)
{
  u <- 0.77957 - 0.33349 * log(a)
  fam <- list(name = "modified.optimal", cc = c(a, 1e-9, Inf, 1.0, 0.0, 0.0))
  upper <- uniroot(f = rhoprime, interval = c(u, 1.5*u), family = fam, check.conv = TRUE, tol = 1e-8)$root
  c(0.0, upper)
}


findTuningConstFromGaussianEfficiency_modified.optimal <- function(e, interval = c(1e-6, 0.2255))
{
  obj <- function(a, e) {
    sup <- psiSupportFromTuningConst(a, "modified.optimal")
    fam <- list(name = "modified.optimal", cc = c(a, DNORM1 / (DNORM1 - a), sup[2], 1.0, NA, NA))
    e - computeGaussianEfficiencyFromFamily(fam, psiSupport = sup)
  }

  uniroot(obj, interval = interval, e = e, check.conv = TRUE, tol = 1e-8)$root
}
