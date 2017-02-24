\name{summary.lts}
\alias{summary.lts}
\alias{print.summary.lts}
\title{Summary Method for LTS objects}
\usage{
\method{summary}{lts}(object, correlation = FALSE, \dots)

\method{print}{summary.lts}(x, digits = max(3, getOption("digits") - 3),
     signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"lts"}, usually, a result of a call to \code{\link{ltsReg}}.}
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed.}
  \item{x}{an object of class \code{"summary.lts"}, usually, a result of a call to \code{summary.lts}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical indicating if \dQuote{significance stars}
    should be printer, see \code{\link{printCoefmat}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  \code{summary} method for class \code{"lts"}.
}
\details{
  These functions compute and print summary statistics for weighted
  least square estimates with weights based on LTS estimates.  Therefore
  the statistics are similar to those for LS but all terms are
  multiplied by the corresponding weight.

  Correlations are printed to two decimal places: to see the actual correlations
  print \code{summary(object)$correlation} directly.
}
\value{
  The function \code{summary.lts} computes and returns a list of summary
  statistics of the fitted linear model given in \code{object}, using
  the components of this object (list elements).
  \item{residuals}{the residuals - a vector like the response \code{y}
    containing the residuals from the weighted least squares regression.}
  \item{coefficients}{a \eqn{p \times 4}{p x 4} matrix with columns for
    the estimated coefficient, its standard error, t-statistic and
    corresponding (two-sided) p-value. }
  \item{sigma}{the estimated scale of the reweighted residuals
    \deqn{\hat\sigma^2 = \frac{1}{n-p}\sum_i{R_i^2},}{%
      sigma^2 = 1/(n-p) Sum(R[i]^2),}
    where \eqn{R_i}{R[i]} is the \eqn{i}-th residual, \code{residuals[i]}.}
  \item{df}{degrees of freedom, a 3-vector \eqn{(p, n-p, p*)}, the last
    being the number of non-aliased coefficients.}
  \item{fstatistic}{(for models including non-intercept terms)
    a 3-vector with the value of the F-statistic with
    its numerator and denominator degrees of freedom.}
  \item{r.squared}{\eqn{R^2}, the \dQuote{fraction of variance explained by
    the model},
    \deqn{R^2 = 1 - \frac{\sum_i{R_i^2}}{\sum_i(y_i- y^*)^2},}{%
      R^2 = 1 - Sum(R[i]^2) / Sum((y[i]- y*)^2),}
    where \eqn{y^*}{y*} is the mean of \eqn{y_i}{y[i]} if there is an
    intercept and zero otherwise.}
  \item{adj.r.squared}{the above \eqn{R^2} statistic
    \dQuote{\emph{adjusted}}, penalizing for higher \eqn{p}.}
  \item{cov.unscaled}{a \eqn{p \times p}{p x p} matrix of (unscaled)
    covariances of the \eqn{\hat\beta_j}{coef[j]}, \eqn{j=1, \dots, p}.}
  \item{correlation}{the correlation matrix corresponding to the above
    \code{cov.unscaled}, if \code{correlation = TRUE} is specified.}
}
\seealso{
  \code{\link{ltsReg}}; the generic \code{\link{summary}}.
}
\examples{
data(Animals2)
ltsA <- ltsReg(log(brain) ~ log(body), data = Animals2)
(slts <- summary(ltsA))
## non-default options for printing the summary:
print(slts, digits = 5, signif.stars = FALSE)
}
\keyword{regression}
\keyword{robust}
