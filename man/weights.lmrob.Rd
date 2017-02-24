\name{weights.lmrob}
\title{Extract Robustness and Model Weights}
\alias{weights.lmrob}
\alias{weights.glmrob}
\description{
  \code{weights()} extracts robustness weights or fitting
  (or prior) weights from a \code{lmrob} or \code{glmrob} object.
}
\usage{
\method{weights}{lmrob}(object, type = c("prior", "robustness"), ...)
}
\arguments{
  \item{object}{
    an object of class \code{"lmrob"} or \code{"glmrob"}, typically the
    result of a call to \code{\link{lmrob}}, or \code{\link{glmrob}},
    respectively.}
  \item{type}{the type of weights to be returned.  Either
    \code{"prior"} (default), or  \code{"robustness"}.}
  \item{\dots}{not used currently.}
}
\details{
  The \dQuote{prior weights} correspond to the weights specified using
  the \dQuote{weights} argument when calling \code{lmrob}. The
  \dQuote{robustness weights} are the weights assigned by the
  M-estimator of regression, \eqn{\psi(r_i/S) / (r_i/S)}. The robust
  coefficient estimate then numericarlly corresponds to a weighted least
  squares fit using the product of both types of weights as weights.
}
\value{
  Weights extracted from the object \code{object}.
}
\author{Manuel Koller and Martin Maechler.}

\seealso{
  \code{\link{lmrob}}, \code{\link{glmrob}} and \code{\link{weights}}
}
