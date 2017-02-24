\name{glmrob..control}
\title{Controlling Robust GLM Fitting by Different Methods}
\alias{glmrobMqle.control}
\alias{glmrobMT.control}
\alias{glmrobBY.control}
\description{
  These are auxiliary functions as user interface for \code{\link{glmrob}} fitting
  when the different methods, \code{"Mqle"}, \code{"BY"}, or
  \code{"MT"} are used.  Typically only used when calling \code{\link{glmrob}}.
}

\usage{
glmrobMqle.control(acc = 1e-04, test.acc = "coef", maxit = 50, tcc = 1.345)
glmrobBY.control  (maxit = 1000, const = 0.5, maxhalf = 10)
glmrobMT.control  (cw = 2.1, nsubm = 500, acc = 1e-06, maxit = 200)
}
\arguments{
  \item{acc}{positive convergence tolerance;
    the iterations converge when ???}
  \item{test.acc}{Only "coef" is currently implemented}
  \item{maxit}{integer giving the maximum number of iterations. }
  \item{tcc}{tuning constant c for Huber's psi-function}
  \item{const}{for "BY", the normalizing constant ..}% FIXME
  \item{maxhalf}{for "BY"; the number of halving steps when the gradient
    itself no longer improves.  We have seen examples when increasing
    \code{maxhalf} was of relevance.}
  \item{cw}{tuning constant c for Tukey's biweight psi-function}
  \item{nsubm}{the number of subsamples to take for finding an initial
    estimate for \code{method = "MT"}.}
}
%% \details{
%% }
\value{
  A \code{\link{list}} with the arguments as components.
}

\author{Andreas Ruckstuhl and Martin Maechler}

\seealso{\code{\link{glmrob}}}

\examples{
str(glmrobMqle.control())
str(glmrobBY.control())
str(glmrobMT.control())
}

\keyword{robust}
\keyword{regression}
\keyword{nonlinear}

