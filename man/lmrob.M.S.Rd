\name{lmrob.M.S}
\alias{lmrob.M.S}
\title{ M-S regression estimators }
\description{
  Computes an M-S-estimator for linear regression using the
  \dQuote{M-S} algorithm.
}
\usage{
lmrob.M.S(x, y, control, mf, split)
}
\arguments{
  \item{x}{numeric matrix (a \code{\link{model.matrix}}) of the
    predictors.}
  \item{y}{numeric vector for the response }
  \item{control}{ list as returned by \code{\link{lmrob.control}}.}
  \item{mf}{a model frame as returned by \code{\link{model.frame}}.}
  \item{split}{(optional) list as returned by \code{\link{splitFrame}}.}
}
\details{
  This function is used by \code{\link{lmrob}} and not intended to be
  used on its own (because an M-S-estimator has too low efficiency
  \sQuote{on its own}).

  An M-S estimator is a combination of an S-estimator for the
  continuous variables and an L1-estimator for the categorical
  variables.

  The S-estimator is estimated using a subsampling algorithm.  If the
  model includes interactions between categorical (\code{\link{factor}})
  and continuous variables, the subsampling algorithm might fail.  In
  this case, one can choose to assign the interaction to the categorical
  side of variables rather than to the continuous side.  This can be
  accomplished via the control argument \code{split.type} or by
  specifying \code{split}, see \code{\link{splitFrame}}.

  Note that the return status \code{converged} does not refer to the
  actual convergence status.  The algorithm used does not guarantee
  convergence and thus true convergence is almost never reached. This
  is, however, not a problem if the estimate is only used as initial
  estimate part of an MM or SMDM estimate.

  The algorithm sometimes produces the warning message \dQuote{Skipping
    design matrix equilibration (dgeequ): row ?? is exactly zero.}. This
    is just an artifact of the algorithm and can be ignored safely.
}
\value{
  A list with components
  \item{coefficients}{numeric vector (length \eqn{p}) of M-S-regression
    coefficient estimates.}
  \item{scale}{the M-S-scale residual estimate}
  \item{residuals}{numeric vector (legnth \eqn{n}) of the residuals.}
  \item{rweights}{numeric vector (length \eqn{n}) of the robustness weights.}
  \item{control}{the same list as the \code{control} argument.}
  \item{converged}{Convergence status (always \code{TRUE}), needed for
    \code{\link{lmrob.fit}}.}
}
\references{
  Maronna, R. A., and Yohai, V. J. (2000).
  Robust regression with both continuous and categorical predictors.
  \emph{Journal of Statistical Planning and Inference} \bold{89}, 197--214.
}
\author{
  Manuel Koller
}

\seealso{
  \code{\link{lmrob}}; for a description of the available split types, see
  \code{\link{splitFrame}}.

  \code{\link[robust]{lmRob}} in package \pkg{robust} uses a version of
  the M-S algorithm automatically when the formula contains factors.
  Our version however follows Maronna and Yohai (2000) more closely.
}
\examples{
data(education)
education <- within(education, Region <- factor(Region))
flm <- lm(Y ~ Region + X1 + X2 + X3, education)
x <- model.matrix(flm)
y <- education$Y # == model.response(model.frame(flm))
set.seed(17)
f.MS <- lmrob.M.S(x, y, control = lmrob.control(),
                  mf = model.frame(flm))

## The typical use of the "M-S" estimator -- as initial estimate :
fmMS <- lmrob(Y ~ Region + X1 + X2 + X3, education,
              init = "M-S")
}

\keyword{ M-S }
\keyword{ robust }
\keyword{ regression }
