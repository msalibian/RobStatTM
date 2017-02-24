\name{summary.nlrob}
\alias{summary.nlrob}

\title{Summarizing Robust Fits of Nonlinear Regression Models }
\description{
  \code{summary} method for objects of class \code{"nlrob"}, i.e.,
  \code{\link{nlrob}()} results.  Currently it only works for
  \code{nlrob(*, method="M")}.
}
\usage{
\method{summary}{nlrob}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
}

\arguments{
  \item{object}{an object of class \code{"nlrob"}, usually, a result of
    a call to \code{\link{nlrob}}.}
  \item{correlation}{logical variable indicating whether
    to compute the correlation matrix of the estimated coefficients.}
  \item{symbolic.cor}{logical indicating whether
    to use symbols to display the above correlation matrix.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  The function \code{\link{summary.nlrob}} computes and returns an
  object of class \code{"summary.nlrob"} of summary statistics of the
  robustly fitted linear model given in \code{object}.
  There is a print method, \code{print.summary.lmrob()}, which nicely
  formats the output.

  The result keeps a large part of \code{object}'s components such as
  \code{residuals}, \code{cov} or \code{w}, and additionally contains
  \item{coefficients}{the matrix of coefficients, standard errors and
    p-values.}
  \item{correlation}{if the \code{correlation} argument was true, the
    correlation matrix of the parameters.}
%% maybe add more
}

\author{Andreas Ruckstuhl}

\seealso{\code{\link{nlrob}()}, also for examples.
}
\keyword{regression}
\keyword{nonlinear}
\keyword{robust}
