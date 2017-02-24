\name{lmrob.fit}
\alias{lmrob.fit}
\alias{lmrob.fit.MM}
\title{ MM-type estimator for regression }
\description{
  Compute MM-type estimators of regression:  An S-estimator is
  used as starting value, and an M-estimator with fixed scale and
  redescending psi-function is used from there. Optionally a D-step
  (Design Adaptive Scale estimate) as well as a second M-step is
  calculated.
}
\usage{
lmrob.fit(x, y, control, init = NULL, mf = NULL)
}
\arguments{
  \item{x}{design matrix (\eqn{n \times p}{n x p}) typically including a
    column of \code{1}s for the intercept.}
  \item{y}{numeric response vector (of length \eqn{n}).}
  \item{control}{a list of control parameters as returned
    by \code{\link{lmrob.control}}, used for both the initial S-estimate
    and the subsequent M- and D-estimates.}
  \item{init}{optional \code{\link{list}} of initial estimates.  See
    \emph{Details}.}
  \item{mf}{(optional) a model frame as returned by \code{\link{model.frame}},
    used only to compute outlier statistics, see \code{\link{outlierStats}}.}
}
\details{This function is the basic fitting function for MM-type estimation,
  called by \code{\link{lmrob}} and typically not to be used on its own.

  If given, \code{init} must be a list of initial estimates containing
  at least the initial coefficients and scale as \code{coefficients} and
  \code{scale}.  Otherwise it calls \code{\link{lmrob.S}(..)} and uses it
  as initial estimator.
}
\value{
  A list with components
  \item{fitted.values}{\eqn{X \beta}{X beta}, i.e., \code{X \%*\% coefficients}.}
  \item{residuals}{the raw residuals, \code{y - fitted.values}}
  \item{rweights}{robustness weights derived from the final M-estimator
    residuals (even when not converged).}
  \item{rank}{}
  \item{degree.freedom}{\code{n - rank}}% more!
  \item{coefficients}{estimated regression coefficient vector}
  \item{scale}{the robustly estimated error standard deviation}% = final.MM$scale,
  \item{cov}{variance-covariance matrix of \code{coefficients}, if the
    RWLS iterations have converged (and \code{control$cov} is not \code{"none"}).}
  \item{control}{}% = control,
  \item{iter}{}% = final.MM$iter, <<<<<< also 'init.S' !
  \item{converged}{logical indicating if the RWLS iterations have converged.}
  \item{init.S}{the whole initial S-estimator result, including its own
    \code{converged} flag, see \code{\link{lmrob.S}} (only for MM-estimates).}
  \item{init}{A similar list that contains the results of intermediate
    estimates (not for MM-estimates).}
}
\author{ Matias Salibian-Barrera, Martin Maechler and Manuel Koller}
\seealso{
  \code{\link{lmrob}},
  \code{\link{lmrob..M..fit}},
  \code{\link{lmrob..D..fit}},
  \code{\link{lmrob.S}}
}
\keyword{robust}
\keyword{regression}
