# The supported rho families. 

#FAMILY.NAMES <- c("bisquare", "ggw", "hampel", "huber", "lqq", "modified.optimal", "optimal", "welsh")
FAMILY.NAMES <- c("bisquare", "modified.optimal", "optimal")


#' Tukey bisquare rho object
#'
#' @param efficiency the desired efficiency of the corresponding regression 
#' estimator for Gaussian errors
#' @param breakdown.point the desired breakdown of the associated regression 
#' estimator
#'
#' @return A list with elements \code{name} (the string 'bisquare') and the corresponding tuning
#' constant \code{cc}
#'
#' @rdname rho
#' @author Kjell Konis
#'
#' @export
bisquare <- function(efficiency, breakdown.point)
{
  if(missing(breakdown.point))
    cc <- c("c" = findTuningConstFromGaussianEfficiency(efficiency, "bisquare"))
  else
    cc <- c("c" = 1.5477)

  list(name = "bisquare", cc = cc)
}

#' Optimal rho function object
#'
#' @param e the desired efficiency of the corresponding regression 
#' estimator for Gaussian errors
#'
#' @return A list with elements \code{name} (the string 'optimal') and the corresponding tuning
#' constant \code{cc}
#'
#' @rdname rho
#' @author Kjell Konis
#'
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


#' Modified optimal rho function object
#'
#' @param e the desired efficiency of the corresponding regression 
#' estimator for Gaussian errors
#'
#' @return A list with elements \code{name} (the string 'modified.optimal') and the corresponding tuning
#' constant \code{cc}
#'
#' @rdname rho
#' @author Kjell Konis
#'
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
#' @param family tuning parameters as computed by 
#' a rho family function. Currently the following are implemented
#' \code{\link{bisquare}}, \code{\link{optimal}}, and \code{\link{modified.optimal}}. 
#'
#' @return The value of \code{rho} at \code{u}
#'
#' @rdname rho
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @export
rho <- function(u, family, cc, standardize = TRUE)
{
  family.name <- match.arg(family, choices = FAMILY.NAMES)

  if(standardize)
    Mchi(u, cc = cc, psi = family.name, deriv = 0)
  else
    Mpsi(u, cc = cc, psi = family.name, deriv = -1)
}


#' The first derivative of Tukeys bisquare rho function
#'
#' @param u scalar or vector at which the derivative of rho is to be evaluated
#' @param family tuning parameters as computed by 
#' a rho family function. Currently the following are implemented
#' \code{\link{bisquare}}, \code{\link{optimal}}, and \code{\link{modified.optimal}}. 
#'
#' @return The value of the first derivative \code{rho} evaluated at \code{u}
#'
#' @rdname rhoprime
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @export
rhoprime <- function(u, family, cc, standardize = FALSE)
{
  family.name <- match.arg(family, choices = FAMILY.NAMES)

  if(standardize)
    Mchi(u, cc = cc, psi = family.name, deriv = 1)
  else
    Mpsi(u, cc = cc, psi = family.name, deriv = 0)
}


#' The second derivative of Tukey bisquare rho function
#'
#' @param u scalar or vector at which the second derivative of rho is to be evaluated
#' @param family tuning parameters as computed by 
#' a rho family function. Currently the following are implemented
#' \code{\link{bisquare}}, \code{\link{optimal}}, and \code{\link{modified.optimal}}. 
#'
#' @return The value of the second derivative of \code{rho} evaluated at \code{u}
#'
#' @rdname rhoprime2
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @export
rhoprime2 <- function(u, family, cc, standardize = FALSE)
{
  family.name <- match.arg(family, choices = FAMILY.NAMES)

  if(standardize)
    Mchi(u, cc = cc, psi = family.name, deriv = 2)
  else
    Mpsi(u, cc = cc, psi = family.name, deriv = 1)
}


