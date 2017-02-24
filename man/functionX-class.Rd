\name{functionX-class}
\docType{class}
\alias{functionX-class}
\title{Class "functionX" of Psi-like Vectorized Functions}
\description{
  The class \code{"functionX"} of vectorized functions of one argument
  \code{x} and typically further tuning parameters.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("functionX", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Directly extends class \code{"function"}.}
  }
}
\section{Extends}{
  Class \code{"function"}, from data part.
  Class \code{"OptionalFunction"}, by class \code{"function"}.
  Class \code{"PossibleMethod"}, by class \code{"function"}.
}
\section{Methods}{
  No methods defined with class "functionX" in the signature.
}
\author{Martin Maechler}
\seealso{
  \code{\link{psiFunc}()}, and class descriptions of
  \code{\linkS4class{functionXal}} for \emph{functionals} of
  \code{"functionX"}, and \code{\linkS4class{psi_func}} which has
  several \code{functionX} slots.
}
% \examples{

% }
\keyword{classes}
