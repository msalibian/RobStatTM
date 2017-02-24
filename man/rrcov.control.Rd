\name{rrcov.control}
\alias{rrcov.control}
\title{Control Settings for covMcd and ltsReg}
\description{
  Auxiliary function for passing the estimation options as parameters to the
  estimation functions.
}
\usage{
rrcov.control(alpha = 1/2, method = c("covMcd", "covComed", "ltsReg"),
              nsamp = 500, nmini = 300, kmini = 5,
              seed = NULL, tolSolve = 1e-14,
              scalefn = "hrv2012", maxcsteps = 200,
              trace = FALSE,
              wgtFUN = "01.original", beta,
              use.correction = identical(wgtFUN, "01.original"),
              adjust = FALSE)
}
\arguments{
  \item{alpha}{This parameter controls the size of the subsets over
    which the determinant is minimized, i.e., \code{alpha*n} observations
    are used for computing the determinant.  Allowed values are between 0.5
    and 1 and the default is 0.5. }
  \item{method}{a string specifying the \dQuote{main} function for which
    \code{rrcov.control()} is used.  This currently only makes a
    difference to determine the default for \code{beta}.}
  \item{nsamp}{number of subsets used for initial estimates or \code{"best"}
    or \code{"exact"}. Default is \code{nsamp = 500}.
    If \code{nsamp="best"} exhaustive enumeration is done, as far as
    the number of trials do not exceed 5000. If \code{nsamp="exact"}
    exhaustive enumeration will be attempted however many samples
    are needed. In this case a warning message will be displayed
    saying that the computation can take a very long time. }
  \item{nmini, kmini}{for \code{\link{covMcd}}: For large \eqn{n}, the algorithm
    splits the data into maximally \eqn{kmini} subsets of targetted size
    \code{nmini}.  See \code{\link{covMcd}} for more details.}%--- ./covMcd.Rd
  \item{seed}{initial seed for R's random number generator; see
    \code{\link{.Random.seed}} and the description of the \code{seed}
    argument in \code{\link{lmrob.control}}.}
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in \code{\link{mahalanobis}}.}
  \item{scalefn}{(for deterministic \code{\link{covMcd}()}:) a character
    string or \code{\link{function}} for computing a robust scale
    estimate.  The current default \code{"hrv2012"} uses the recommendation
    of Hubert et al (2012); see \code{\link{covMcd}} for more.}
  \item{maxcsteps}{integer specifying the maximal number of
    concentration steps for the deterministic MCD.}
  \item{trace}{logical or integer indicating whether to print
    intermediate results.  Default is \code{trace = FALSE}.}
  \item{wgtFUN}{a character string or \code{\link{function}}, specifying
    how the weights for the reweighting step should be computed, see
    \code{\link{ltsReg}}, \code{\link{covMcd}} or
    \code{\link{covComed}}, respectively.  The default is specified by
    \code{"01.original"}, as the resulting weights are 0 or 1.  Alternative
    string specifications need to match \code{names(.wgtFUN.covComed)} -
    which currently is experimental.}
  \item{beta}{a quantile, experimentally used for some of the prespecified
    \code{wgtFUN}s, see e.g., \code{\link{.wgtFUN.covMcd}} and
    \code{\link{.wgtFUN.covComed}}.}
  \item{use.correction}{whether to use finite sample correction factors.
    Defaults to \code{TRUE}.}
  \item{adjust}{(for \code{\link{ltsReg}()}:) whether to perform
    intercept adjustment at each step.  Because this can be quite time
    consuming, the default is \code{adjust = FALSE}.}
}
\seealso{
  For details, see the documentation about \code{\link{ltsReg}} and
  \code{\link{covMcd}}, respectively.
}
\value{
  A list with components, as the parameters passed by the invocation
}
\author{Valentin Todorov}
\examples{
data(Animals, package = "MASS")
brain <- Animals[c(1:24, 26:25, 27:28),]
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])

ctrl <- rrcov.control(alpha=0.75, trace=TRUE)
covMcd(hbk.x,      control = ctrl)
covMcd(log(brain), control = ctrl)
}
\keyword{robust}
\keyword{multivariate}
