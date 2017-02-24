\name{lmrob.S}
\alias{lmrob.S}
\title{ S-regression estimators }
\description{
  Computes an S-estimator for linear regression,
  using the \dQuote{fast S} algorithm.% of Matias Salibian & Victor Yohai ...
}
\usage{
lmrob.S(x, y, control, trace.lev = control$trace.lev, mf = NULL)
}
\arguments{
  \item{x}{ design matrix }
  \item{y}{ response vector }
  \item{control}{ list as returned by \code{\link{lmrob.control}} }
  \item{trace.lev}{integer indicating if the progress of the algorithm
    should be traced (increasingly); default \code{trace.lev = 0} does
    no tracing.}
  \item{mf}{(optional) a model frame as returned by \code{\link{model.frame}}, 
    used only to compute outlier statistics, see \code{\link{outlierStats}}.}
}
\details{
  This function is used by \code{\link{lmrob.fit}} and not
  intended to be used on its own (because an S-estimator has too low
  efficiency \sQuote{on its own}).

  By default, the subsampling algorithm uses a customized LU
  decomposition which ensures a non singular subsample (if this is at
  all possible). This makes the Fast-S algorithm also feasible for
  categorical and mixed continuous-categorical data.

  One can revert to the old subsampling scheme by setting the parameter
  \code{subsampling} in \code{control} to \code{"simple"}.
}
\value{
  A list with components
  \item{coefficients}{numeric vector (length \eqn{p}) of S-regression coefficient estimates.}
  \item{scale}{the S-scale residual estimate}% 'residual estimate' ?? % resid. VAR !?
  % \item{cov}{covariance matrix (\eqn{p \times p}{p x p}) of the
  %   coefficient estimates.}
  \item{fitted.values}{numeric vector (length \eqn{n}) of the fitted
    values.}
  \item{residuals}{numeric vector (length \eqn{n}) of the residuals.}
  \item{rweights}{numeric vector (length \eqn{n}) of the robustness weights.}
  \item{k.iter}{(maximal) number of refinement iterations used.}
  \item{converged}{logical indicating if \bold{all} refinement
    iterations had converged.}
  \item{control}{the same list as the \code{control} argument.}
}
\seealso{\code{\link{lmrob}}, also for references.
}
\author{ Matias Salibian-Barrera and Manuel Koller (and Martin Maechler
  for minor details)
}
\examples{
set.seed(33)
x1 <- sort(rnorm(30)); x2 <- sort(rnorm(30)); x3 <- sort(rnorm(30))
X. <- cbind(x1, x2, x3)
y <-  10 + X. \%*\% (10*(2:4)) + rnorm(30)/10
y[1] <- 500   # a moderate outlier
X.[2,1] <- 20 # an X outlier
X1  <- cbind(1, X.)

(m.lm <- lm(y ~ X.))
set.seed(12)
m.lmS <- lmrob.S(x=X1, y=y,
                 control = lmrob.control(nRes = 20), trace.lev=1)
m.lmS[c("coefficients","scale")]
all.equal(unname(m.lmS$coef), 10 * (1:4), tolerance = 0.005)
stopifnot(all.equal(unname(m.lmS$coef), 10 * (1:4), tolerance = 0.005),
          all.equal(m.lmS$scale, 1/10, tolerance = 0.09))
}
\keyword{robust}
\keyword{regression}
