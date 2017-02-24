\name{BYlogreg}
\alias{BYlogreg}
\title{Bianco-Yohai Estimator for Robust Logistic Regression}
\encoding{utf8}
\description{
  Computation of the estimator of Bianco and Yohai (1996) in logistic regression.
  Now provides both the \emph{weighted} and regular (unweighted) BY-estimator.

  By default, an intercept term is included and p parameters are estimated.
  For more details, see the reference.

  Note: This function is for \dQuote{back-compatibility} with the
  \code{BYlogreg()} code web-published at
  \url{http://www.econ.kuleuven.be/public/NDBAE06/programs/roblog/};
  instead, the recommended interface is \code{\link{glmrob}(*, method =
  "BY")} or \code{... method = "WBY" ..}, see \code{\link{glmrob}}.
}
\usage{
BYlogreg(x0, y, initwml = TRUE, addIntercept = TRUE,
         const = 0.5, kmax = 1000, maxhalf = 10, sigma.min = 1e-4,
         trace.lev = 0)
}
\arguments{
  \item{x0}{a numeric \eqn{n \times (p-1)}{n * (p-1)} matrix containing
    the explanatory variables.}
  \item{y}{numeric \eqn{n}-vector of binomial (0 - 1) responses.}
  \item{initwml}{logical for selecting one of the two possible methods
    for computing the initial value of the optimization process.

    If \code{initwml} is true (default), a weighted ML estimator is
    computed with weights derived from the MCD estimator
    computed on the explanatory variables.

    If \code{initwml} is false, a classical ML fit is perfomed.  When
    the explanatory variables contain binary observations, it is
    recommended to set initwml to FALSE or to modify the code of the
    algorithm to compute the weights only on the continuous variables.
  }
  \item{addIntercept}{logical indicating that a column of \code{1} must be
    added the \eqn{x} matrix.}
  \item{const}{tuning constant used in the computation of the estimator
    (default=0.5).}
  \item{kmax}{maximum number of iterations before convergence (default=1000).}
  \item{maxhalf}{max number of step-halving (default=10).}
  \item{sigma.min}{smallest value of the scale parameter before
    implosion (and hence non-convergence) is assumed.}
  \item{trace.lev}{logical (or integer) indicating if intermediate results
    should be printed; defaults to \code{0} (the same as \code{FALSE}).}
}
%% \details{
%%   If necessary, more details than the description above
%% }
\value{
  a list with components
  \item{convergence}{logical indicating if convergence was achieved}
  \item{objective}{the value of the objective function at the minimum}
  \item{coefficients}{vector of parameter estimates}
  \item{vcov}{variance-covariance matrix of the coefficients (if convergence is TRUE).}
  \item{sterror}{standard errors, i.e., simply \code{sqrt(diag(.$vcov))}, if convergence.}
}
\references{
  Croux, C., and Haesbroeck, G. (2003)
  Implementing the Bianco and Yohai estimator for Logistic Regression,
  \emph{Computational Statistics and Data Analysis} \bold{44}, 273--295.

  Ana M. Bianco and Víctor J. Yohai (1996)
  Robust estimation in the logistic regression model.
  In Helmut Rieder, \emph{Robust Statistics, Data Analysis, and
    Computer Intensive Methods}, Lecture Notes in Statistics \bold{109},
  pages 17--34.
}
\author{
  Originally, Christophe Croux and Gentiane Haesbroeck, with
  thanks to Kristel Joossens and Valentin Todorov for improvements.

  Speedup, tweaks, more \dQuote{control} arguments: Martin Maechler.
}
\seealso{
  The more typical way to compute BY-estimates (via
  \code{\link{formula}} and methods):
  \code{\link{glmrob}(*, method = "WBY")} and \code{.. method = "BY"}.
}
\examples{
set.seed(17)
x0 <- matrix(rnorm(100,1))
y  <- rbinom(100, size=1, prob= 0.5) # ~= as.numeric(runif(100) > 0.5)
BY <- BYlogreg(x0,y)
BY <- BYlogreg(x0,y, trace.lev=TRUE)

## The "Vaso Constriction"  aka "skin" data:
data(vaso)
vX <- model.matrix( ~ log(Volume) + log(Rate), data=vaso)
vY <- vaso[,"Y"]
head(cbind(vX, vY))# 'X' does include the intercept

vWBY <- BYlogreg(x0 = vX, y = vY, addIntercept=FALSE) # as 'vX' has it already
v.BY <- BYlogreg(x0 = vX, y = vY, addIntercept=FALSE, initwml=FALSE)
## they are relatively close:
stopifnot( all.equal(vWBY, v.BY, tolerance = 2e-4) )
}
\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
