\name{chgDefaults-methods}
\title{Change Defaults (Parameters) of "Psi Function" Objects}
\docType{methods}
\alias{chgDefaults}% the generic
\alias{chgDefaults-methods}
\alias{chgDefaults,ANY-method}
\alias{chgDefaults,psi_func-method}
\description{
  To modify an object of class \code{\linkS4class{psi_func}},
  i.e. typically change the tuning parameters, the generic function
  \code{chgDefaults()} is called and works via the corresponding method.
}
\section{Methods}{
  \describe{
    \item{\code{object = "psi_func"}}{The method is used to change the default
      values for the tuning parameters, and returns an object of class
      \code{\linkS4class{psi_func}}, a copy of input \code{object}
      with the slot \code{tDefs} possibly changed;.
    }
  }
}
\seealso{
  \code{\link{psiFunc}}
}
\examples{
## Hampel's psi and rho:
H.38 <- chgDefaults(hampelPsi, k = c(1.5, 3.5, 8))
H.38
plot(H.38)
## for more see  ?psiFunc
}
\keyword{methods}

