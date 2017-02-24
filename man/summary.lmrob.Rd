\name{summary.lmrob}
\title{Summary Method for "lmrob" Objects}
%
\alias{summary.lmrob}
\alias{hatvalues.lmrob}
\alias{.lmrob.hat}
\alias{vcov.lmrob}
\alias{print.summary.lmrob}
\alias{model.matrix.lmrob}
%
\description{
  Summary method for \R object of class \code{"lmrob"} and
  \code{\link{print}} method for the summary object.

  Further, methods \code{\link{fitted}()}, \code{\link{residuals}()}
  work (via the default methods), and
  \code{\link{predict}()} (see \code{\link{predict.lmrob}},
  \code{\link{vcov}()}, \code{\link{weights}()} (see
  \code{\link{weights.lmrob}}), \code{\link{model.matrix}()},
  \code{\link{confint}()}, \code{\link{dummy.coef}()},
  \code{\link{hatvalues}()}, etc.,
  have explicitly defined \code{lmrob} methods.  \code{.lmrob.hat()} is
  the lower level \dQuote{work horse} of the \code{hatvalues()} method.
}
\usage{
\method{summary}{lmrob}(object, correlation = FALSE,
        symbolic.cor = FALSE, \dots)
\method{print}{summary.lmrob}(x, digits = max(3, getOption("digits") - 3),
      symbolic.cor= x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"),
      showAlgo = TRUE, \dots)

\method{vcov}{lmrob}(object, cov = object$control$cov, \dots)
\method{model.matrix}{lmrob}(object, \dots)

% not yet
% .lmrob.hat(x, w = rep(1, NROW(x)), wqr = qr(sqrt(w) * x))
}
\arguments{
  \item{object}{an \R object of class \code{lmrob}, typically created by
    \code{\link{lmrob}}.}
  \item{correlation}{logical variable indicating whether
    to compute the correlation matrix of the estimated coefficients.}
  \item{symbolic.cor}{logical indicating whether
    to use symbols to display the above correlation matrix.}

  \item{x}{an \R object of class \code{summary.lmrob}, typically
    resulting from \code{summary(\link{lmrob}(..),..)}.}
  \item{digits}{number of digits for printing, see \code{digits} in
    \code{\link{options}}.}
  \item{signif.stars}{logical variable indicating
    whether to use stars to display different levels of
    significance in the individual t-tests.}
  \item{showAlgo}{optional \code{\link{logical}} indicating if the
    algorithmic parameters (as mostly inside the \code{control} part)
    should be shown.}
  \item{cov}{covariance estimation function to use, a
    \code{\link{function}} or \link{character} string naming the
    function; \pkg{robustbase} currently provides \code{".vcov.w"} and
    \code{".vcov.avar1"}, see \emph{Details} of \code{\link{lmrob}}.
    Particularly useful when \code{object} is the result of
    \code{lmrob(.., cov = "none")}, where \preformatted{  object$cov <- vcov(object, cov = ".vcov.w")}
    allows to \emph{update} the fitted object.}
  \item{\dots}{potentially more arguments passed to methods.}
}
\value{
  \code{summary(object)} returns an object of S3 class
  \code{"summary.lmrob"}, basically a \code{\link{list}} with components
  "call", "terms", "residuals", "scale", "rweights", "converged",
  "iter", "control" all copied from \code{object}, and further
  components, partly for compatibility with \code{\link{summary.lm}},
  \item{coefficients}{a \code{\link{matrix}} with columns \code{"Estimate"},
    \code{"Std. Error"}, \code{"t value"}, and \code{"PR(>|t|)"}, where
    "Estimate" is identical to \code{\link{coef}(object)}. Note that
    \code{\link{coef}(<summary.obj>)} is slightly preferred to access
    this matrix.}
  \item{df}{degrees of freedom, in an \code{\link{lm}} compatible way.}
  \item{sigma}{identical to \code{\link{sigma}(object)}.}
  \item{aliased}{..}%FIXME
  \item{cov}{derived from \code{object$cov}.}% FIXME: say more
  \item{r.squared}{robust \dQuote{R squared} or \eqn{R^2}, a coefficient
    of determination:  This is the consistency corrected robust
    coefficient of determination by Renaud and Victoria-Feser (2010).}
  \item{adj.r.squared}{an adjusted R squared, see \code{r.squared}.}
}
\references{
  Renaud, O. and Victoria-Feser, M.-P. (2010).
  A robust coefficient of determination for regression,
  \emph{Journal of Statistical Planning and Inference} \bold{140}, 1852-1862.
}
\seealso{\code{\link{lmrob}}, \code{\link{predict.lmrob}},
  \code{\link{weights.lmrob}}, \code{\link{summary.lm}},
  \code{\link{print}}, \code{\link{summary}}.
}
\examples{
mod1 <- lmrob(stack.loss ~ ., data = stackloss)
sa <- summary(mod1)  # calls summary.lmrob(....)
sa                   # dispatches to call print.summary.lmrob(....)

## correlation between estimated coefficients:
cov2cor(vcov(mod1))

cbind(fit = fitted(mod1), resid = residuals(mod1),
      wgts= weights(mod1, type="robustness"),
      predict(mod1, interval="prediction"))

data(heart)
sm2 <- summary( m2 <- lmrob(clength ~ ., data = heart) )
sm2
}
\keyword{robust}
\keyword{regression}
