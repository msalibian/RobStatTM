# The supported rho families.

#FAMILY.NAMES <- c("bisquare", "ggw", "hampel", "huber", "lqq", "modopt", "optimal", "welsh")
FAMILY.NAMES <- c("bisquare", "modopt", "optimal")


#' Tuning parameter the rho loss functions
#'
#' This function computes the tuning constant that yields an MM-regression
#' estimator with a desired asymptotic efficiency when computed with a
#' rho function in the corresponding family. The output of this
#' function can be passed to the functions \link{lmrobdet.control},
#' \link{mscale} and \link{rho}.
#'
#' @param e the desired efficiency of the corresponding regression
#' estimator for Gaussian errors
#'
#' @return A length-1 vector with the corresponding tuning constant.
#'
#' @rdname bisquare
#' @author Kjell Konis
#'
#' @examples
#' # Tuning parameters for an 85%-efficient M-estimator at a Gaussian model
#' bisquare(.85)
#'
#' @export
bisquare <- function(e) #, breakdown.point)
{
  findTuningConstFromGaussianEfficiency(e, "bisquare")
}

#' Tuning parameter for a rho function in the (asymptotic bias-) optimal family
#'
#' This function computes the tuning constant that yields an MM-regression
#' estimator with a desired asymptotic efficiency when computed with a
#' rho function in the corresponding family. The output of this
#' function can be passed to the functions \link{lmrobdet.control},
#' \link{mscale} and \link{rho}.
#'
#' @param e the desired efficiency of the corresponding regression
#' estimator for Gaussian errors
#'
#' @return A vector with named elements containing the corresponding tuning
#' parameters.
#'
#' @author Kjell Konis
#'
#' @examples
#' # Tuning parameters for an 85%-efficient M-estimator at a Gaussian model
#' optimal(.85)
#'
#' @export
optimal <- function(e)
{
  if( e > .9999) {
    e <- .9999
    warning("Current implementation of \'optimal\' or \'modopt\' only allows efficiencies up to 99.99%. Efficiency set to 99.99% for this call.")
  }
  a <- findTuningConstFromGaussianEfficiency(e, "optimal")
  cc <- c(a, psiSupportFromTuningConst(a, "optimal"), 1.0)
  cc[5] <- Psi_optimal(cc[2], cc[1])
  cc[6] <- Psi_optimal(cc[3], cc[1]) - cc[5]
  names(cc) <- c("a", "lower", "upper", "c", "Psi(lower)", "rho(Inf)")
  cc
}


#' Tuning parameter for a rho function in the modified (asymptotic bias-) optimal family
#'
#' This function computes the tuning constant that yields an MM-regression
#' estimator with a desired asymptotic efficiency when computed with a
#' rho function in the corresponding family. The output of this
#' function can be passed to the functions \link{lmrobdet.control},
#' \link{mscale} and \link{rho}.
#'
#' @param e the desired efficiency of the corresponding regression
#' estimator for Gaussian errors
#'
#' @return A vector with named elements containing the corresponding tuning
#' parameters.
#'
#' @author Kjell Konis
#'
#' @examples
#' # Tuning parameters for an 85%-efficient M-estimator at a Gaussian model
#' modopt(.85)
#'
#' @export
modopt <- function(e)
{
  if( e > .9999) {
    e <- .9999
    warning("Current implementation of \'optimal\' or \'modopt\' only allows efficiencies up to 99.99%. Efficiency set to 99.99% for this call.")
  }
  a <- findTuningConstFromGaussianEfficiency(e, "modopt")
  cc <- c(a, DNORM1 / (DNORM1 - a), psiSupportFromTuningConst(a, "modopt")[2], 1.0)
  cc[5] <- Psi_optimal(1.0, cc[1])
  cc[6] <- (0.5 + cc[2] * (Psi_optimal(cc[3], cc[1]) - cc[5]))
  names(cc) <- c("a", "normConst", "upper", "c", "Psi(1)", "rho(Inf)")
  cc
}



#' Rho functions
#'
#' This function returns the value of the "rho" loss function used
#' to compute either an M-scale estimator or a robust regression
#' estimator. It currently can be used to compute the bisquare, optimal
#' and modified optimal loss functions.
#'
#' @param u point or vector at which rho is to be evaluated
#' @param family family string specifying the name of the family of loss function to be used (current valid
#' options are "bisquare", "optimal" and "modopt").
#' @param cc tuning parameters to be computed according to efficiency and / or breakdown
#' considerations. See \link{lmrobdet.control}, \link{bisquare}, \link{modopt}
#' and \link{optimal}.
#' @param standardize logical value determining whether the rho function is to be
#' standardized so that its maximum value is 1. See \link{Mpsi}.
#'
#' @return The value(s) of \code{rho} at \code{u}
#'
#' @rdname rho
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @examples
#' # Evaluate rho tuned for 85% efficiency
#' rho(u=1.1, family='bisquare', cc=bisquare(.85))
#' # Evaluate rho tuned for 50% breakdown
#' rho(u=1.1, family='optimal', cc=lmrobdet.control(bb=.5, family='optimal')$tuning.chi)
#'
#' @export
rho <- function(u, family=" bisquare", cc, standardize = TRUE)
{
  family.name <- match.arg(family, choices = FAMILY.NAMES)

  if(standardize)
    Mchi(u, cc = cc, psi = family.name, deriv = 0)
  else
    Mpsi(u, cc = cc, psi = family.name, deriv = -1)
}


#' The first derivative of the rho function
#'
#' @param u point or vector at which rho is to be evaluated
#' @param family family string specifying the name of the family of loss function to be used (current valid
#' options are "bisquare", "optimal" and "modopt").
#' @param cc tuning parameters to be computed according to efficiency and / or breakdown
#' considerations. See \link{lmrobdet.control}, \link{bisquare}, \link{modopt}
#' and \link{optimal}.
#' @param standardize logical value determining whether the rho function is to be
#' standardized so that its maximum value is 1. See \link{Mpsi}.
#'
#' @return The value of the first derivative \code{rho} evaluated at \code{u}
#'
#' @rdname rhoprime
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @examples
#' # Evaluate the derivative of a rho function tuned for 85% efficiency
#' rhoprime(u=1.1, family='bisquare', cc=bisquare(.85))
#' # Evaluate the derivative of a rho function tuned for 50% breakdown
#' rhoprime(u=1.1, family='optimal', cc=lmrobdet.control(bb=.5, family='optimal')$tuning.chi)
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


#' The second derivative of the rho function
#'
#' @param u point or vector at which rho is to be evaluated
#' @param family family string specifying the name of the family of loss function to be used (current valid
#' options are "bisquare", "optimal" and "modopt").
#' @param cc tuning parameters to be computed according to efficiency and / or breakdown
#' considerations. See \link{lmrobdet.control}, \link{bisquare}, \link{modopt}
#' and \link{optimal}.
#' @param standardize logical value determining whether the rho function is to be
#' standardized so that its maximum value is 1. See \link{Mpsi}.
#'
#' @return The value of the second derivative of \code{rho} evaluated at \code{u}
#'
#' @rdname rhoprime2
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @examples
#' # Evaluate the 2nd derivative of a rho function tuned for 85% efficiency
#' rhoprime2(u=1.1, family='bisquare', cc=bisquare(.85))
#' # Evaluate the 2nd derivative of a rho function tuned for 50% breakdown
#' rhoprime2(u=1.1, family='optimal', cc=lmrobdet.control(bb=.5, family='optimal')$tuning.chi)
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


