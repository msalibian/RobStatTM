# The supported rho families. 

#FAMILY.NAMES <- c("bisquare", "ggw", "hampel", "huber", "lqq", "modified.optimal", "optimal", "welsh")
FAMILY.NAMES <- c("bisquare", "modified.optimal", "optimal")


#' @export
bisquare <- function(efficiency, breakdown.point)
{
  if(missing(breakdown.point))
    cc <- c("c" = findTuningConstFromGaussianEfficiency(efficiency, "bisquare"))
  else
    cc <- c("c" = 1.5477)

  list(name = "bisquare", cc = cc)
}

#' @export
optimal <- function(e)
{
  a <- findTuningConstFromGaussianEfficiency(e, "optimal")
  cc <- c(a, psiSupportFromTuningConst(a, "optimal"), 1.0)
  cc[5] <- Psi_optimal(cc[2], cc[1])
  cc[6] <- Psi_optimal(cc[3], cc[1]) - cc[5]
  names(cc) <- c("a", "lower", "upper", "c", "Psi(lower)", "rho(Inf)")
  list(name = "optimal", cc = cc)
}


#' @export
modified.optimal <- function(e)
{
  a <- findTuningConstFromGaussianEfficiency(e, "modified.optimal")
  cc <- c(a, DNORM1 / (DNORM1 - a), psiSupportFromTuningConst(a, "modified.optimal")[2], 1.0)
  cc[5] <- Psi_optimal(1.0, cc[1])
  cc[6] <- (0.5 + cc[2] * (Psi_optimal(cc[3], cc[1]) - cc[5]))
  names(cc) <- c("a", "normConst", "upper", "c", "Psi(1)", "rho(Inf)")
  list(name = "modified.optimal", cc = cc)
}



#' Tukey bisquare rho function
#'
#' @param u point or vector at which rho is to be evaluated
#' @param cc tuning parameter
#'
#' @return The value of \code{rho_cc} at \code{u}
#'
#' @rdname rho
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @export
rho <- function(u, family, standardize = TRUE)
{
  family.name <- match.arg(family$name, choices = FAMILY.NAMES)

  if(standardize)
    Mchi(u, cc = family$cc, psi = family.name, deriv = 0)
  else
    Mpsi(u, cc = family$cc, psi = family.name, deriv = -1)
}


#' The first derivative of Tukeys bisquare rho function
#'
#' @param r scalar or vector at which the derivative of rho is to be evaluated
#' @param cc tuning parameter
#'
#' @return The value of the first derivative \code{rho_cc} evaluated at \code{r}
#'
#' @rdname rhoprime
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @export
rhoprime <- function(u, family, standardize = FALSE)
{
  family.name <- match.arg(family$name, choices = FAMILY.NAMES)

  if(standardize)
    Mchi(u, cc = family$cc, psi = family.name, deriv = 1)
  else
    Mpsi(u, cc = family$cc, psi = family.name, deriv = 0)
}


#' The second derivative of Tukey bisquare rho function
#'
#' @param r scalar or vector at which the second derivative of rho is to be evaluated
#' @param cc tuning parameter
#'
#' @return The value of the second derivative of \code{rho_cc} evaluated at \code{r}
#'
#' @rdname rhoprime2
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @export
rhoprime2 <- function(u, family, standardize = FALSE)
{
  family.name <- match.arg(family$name, choices = FAMILY.NAMES)

  if(standardize)
    Mchi(u, cc = family$cc, psi = family.name, deriv = 2)
  else
    Mpsi(u, cc = family$cc, psi = family.name, deriv = 1)
}


