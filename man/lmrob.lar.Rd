\name{lmrob.lar}
\alias{lmrob.lar}
\title{Least Absolute Residuals / L1 Regression}
\description{
  To compute least absolute residuals (LAR) or  \dQuote{L1} regression,
  \code{lmrob.lar} implements the routine L1 in Barrodale and Roberts (1974),
  which is based on the simplex method of linear programming.  It is a
  copy of \code{lmRob.lar} (in early 2012) from the \pkg{robust} package.
}
\usage{
lmrob.lar(x, y, control, mf)
}
\arguments{
  \item{x}{numeric matrix for the predictors.}
  \item{y}{numeric vector for the response.}
  \item{control}{\code{\link{list}} as returned by
    \code{\link{lmrob.control}()} .}
  \item{mf}{dummy parameter.}
}
\details{
  This method is used for computing the M-S estimate and typically not
  to be used on its own.

  A description of the Fortran subroutines used can be found in Marazzi
  (1993).  In the book, the main method is named \code{RILARS}.
}
\value{
  A list that includes the following components:
  \item{coef }{The L1-estimate of the coefficient vector}
  \item{scale }{The residual scale estimate (mad)}
  \item{resid }{The residuals}
  \item{iter }{The number of iterations required by the simplex
    algorithm}
  \item{status }{Return status (0: optimal, but non unique solution, 1:
    optimal unique solution)}
  \item{converged }{Convergence status (always \code{TRUE}), needed for
    \code{\link{lmrob.fit}}.}
}
\references{
  Marazzi, A. (1993).
  \emph{Algorithms, routines, and S functions for robust statistics}.
  Wadsworth & Brooks/Cole, Pacific Grove, CA.
}
\author{
  Manuel Koller
}
\seealso{
  \code{\link[quantreg]{rq}} from package \pkg{quantreg}.
}
\examples{
data(stackloss)
X <- model.matrix(stack.loss ~ . , data = stackloss)
y <- stack.loss
lmrob.lar(X, y)
}
\keyword{ L1 }
\keyword{ regression }
