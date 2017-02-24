\name{outlierStats}
\alias{outlierStats}
\title{Robust Regression Outlier Statistics}
\description{
  Simple statistics about observations with robustness
  weight of almost zero for models that include factor
  terms. The number of rejected observations and the mean
  robustness weights are computed for each level of each
  factor included in the model.
}
\usage{
outlierStats(object, x = object$x, control = object$control,
             epsw = control$eps.outlier, epsx = control$eps.x,
             warn.limit.reject = control$warn.limit.reject,
             warn.limit.meanrw = control$warn.limit.meanrw)
}
\arguments{
  \item{object}{object of class \code{"lmrob"}, typically the result of a call
    to \code{\link{lmrob}}.}
  \item{x}{design matrix}
  \item{control}{list as returned by \code{\link{lmrob.control}}.}
  \item{epsw}{limit on the robustness weight below which an observation is considered
    to be an outlier.  Either a \code{numeric(1)} or a
    \code{\link{function}} that takes the number of observations as an argument.}
  \item{epsx}{limit on the absolute value of the elements of the design matrix
    below which an element is considered zero.  Either a numeric(1) or a
    function that takes the maximum absolute value in the design matrix
    as an argument.}
  \item{warn.limit.reject}{limit of ratio
    \eqn{\#\mbox{rejected} / \#\mbox{obs in level}}{# rejected / # obs in level}
    above (\eqn{\geq}{>=}) which a warning is produced.  Set to
    \code{NULL} to disable warning.}
  \item{warn.limit.meanrw}{limit of the mean robustness per factor level
    below which (\eqn{\leq}{<=}) a warning is produced.  Set to
    \code{NULL} to disable warning.}
}
\details{
  For models that include factors, the fast S-algorithm used by
  \code{\link{lmrob}} can produce \dQuote{bad} fits for some of the
  factor levels, especially if there are many levels with only a
  few observations. Such a \dQuote{bad} fit is characterized as a
  fit where most of the observations in a level of a factor are
  rejected, i.e., are assigned robustness weights of zero or nearly
  zero. We call such a fit a \dQuote{local exact fit}.

  If a local exact fit is detected, then we recommend to increase some
  of the control parameters of the \dQuote{fast S}-algorithm. As a first
  aid solution in such cases, one can use \code{setting="KS2014"}, see also
  \code{\link{lmrob.control}}.

  This function is called internally by \code{\link{lmrob}} to issue a
  warning if a local exact fit is detected. The output is available as
  \code{ostats} in objects of class \code{"lmrob"} (only if the statistic
  is computed).
}
\value{
  A data frames for each column with any zero elementes as well as an
  overall statistic. The data frame consist of the names of the
  coefficients in question, the number of non-zero observation in that
  level (\code{N.nonzero}), the number of rejected observations
  (\code{N.rejected}), the ratio of rejected observations to the
  number of observations in that level (\code{Ratio}) and the mean
  robustness weight of all the observations in the corresponding level
  (\code{Mean.RobWeight}).
}
\references{
  Koller, M. and Stahel, W.A. (2014)
  Nonsingular subsampling for regression S~estimators with categorical
  predictors. Under review.
}
\author{Manuel Koller}
\seealso{
  \code{\link{lmrob.control}} for the default values of the control
  parameters;  \code{\link{summarizeRobWeights}}.
}
\examples{
## artificial data example
data <- expand.grid(grp1 = letters[1:5], grp2 = letters[1:5], rep=1:3)
set.seed(101)
data$y <- c(rt(nrow(data), 1))
## compute outlier statistics for all the estimators
control <- lmrob.control(method = "SMDM",
                         compute.outlier.stats = c("S", "MM", "SMD", "SMDM"))
## warning is only issued for some seeds
set.seed(2)
fit1 <- lmrob(y ~ grp1*grp2, data, control = control)
## do as suggested:
fit2 <- lmrob(y ~ grp1*grp2, data, setting = "KS2014")

## the plot function should work for such models as well
plot(fit1)

\dontrun{
  ## access statistics:
  fit1$ostats ## SMDM
  fit1$init$ostats ## SMD
  fit1$init$init$ostats ## SM
  fit1$init$init$init.S$ostats ## S
}%dont
}
\keyword{robust}
\keyword{regression}
