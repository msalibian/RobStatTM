\name{covComed}
\title{Co-Median Location and Scatter "Covariance" Estimator}
\alias{covComed}
%%TODO: these two are not yet imported and hence not yet documented below:
\alias{comedian}
\alias{COM}
\alias{.wgtFUN.covComed}
%
\concept{Co-median}
\concept{Comedian}
\description{
  Compute (versions of) the (multivariate) \dQuote{Comedian} covariance,
  i.e., multivariate location and scatter estimator
}
\usage{
covComed(X, n.iter = 2, reweight = FALSE, tolSolve = control$tolSolve,
         trace = control$trace, wgtFUN = control$wgtFUN,
         control = rrcov.control())
% comedian(...)
% COM(...)
}
\arguments{
  \item{X}{data matrix of dimension, say \eqn{n \times p}{n x p}.}
  \item{n.iter}{number of comedian() iterations.  Can be as low as zero.}
  \item{reweight}{logical indicating if the final distances and weights
    should be recomputed from the final \code{cov} and \code{center}.
    The default is currently \code{FALSE} because that was implicit in
    the first version of the \R code.}
  \item{tolSolve}{a numerical tolerance passed to \code{\link{solve}}.}
  \item{trace}{logical (or integer) indicating if intermediate results
    should be printed; defaults to \code{FALSE}; values \eqn{\ge 2}{>= 2}
    also produce print from the internal (Fortran) code.}
  \item{wgtFUN}{a character string or \code{\link{function}}, specifying
    how the weights for the reweighting step should be computed.  The default,
    \code{wgtFUN = "01.original"} corresponds to 0-1 weights as proposed
    originally.  Other predefined string options are available, though
    experimental, see the experimental \code{.wgtFUN.covComed} object.}
  \item{control}{a list with estimation options - this includes those
    above provided in the function specification, see
    \code{\link{rrcov.control}} for the defaults.  If \code{control} is
    supplied, the parameters from it will be used.  If parameters are
    passed also in the invocation statement, they will override the
    corresponding elements of the control object.}
}
\details{
  .. not yet ..
}
\value{
  an object of class \code{"covComed"} which is basically a list with components
 \item{comp1 }{Description of 'comp1'}
 \item{comp2 }{Description of 'comp2'}
 ... FIXME ...
}
\references{
  Falk, M. (1997)
  On mad and comedians.
  \emph{Annals of the Institute of Statistical Mathematics} \bold{49}, 615--644.

  Falk, M. (1998).
  A note on the comedian for elliptical distributions.
  \emph{Journal of Multivariate Analysis} \bold{67}, 306--317.
}
\author{
  Maria Anna di Palma (initial), Valentin Todorov and Martin Maechler
}
\seealso{
 \code{\link{covMcd}}, etc
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
(cc1 <- covComed(hbk.x))
(ccW <- covComed(hbk.x, reweight=TRUE))
cc0  <- covComed(hbk.x, n.iter=0)
cc0W <- covComed(hbk.x, n.iter=0, reweight=TRUE)

stopifnot(all.equal(unclass(cc0), # here, the 0-1 weights don't change:
                    cc0W[names(cc0)], tol=1e-12),
          which(cc1$weights == 0) == 1:14,
          which(ccW$weights == 0) == 1:14,
          which(cc0$weights == 0) == 1:14)
%% TODO : more

## Martin's smooth reweighting:

## List of experimental pre-specified wgtFUN() creators:
## Cutoffs may depend on  (n, p, control$beta) :
str(.wgtFUN.covComed)
}
