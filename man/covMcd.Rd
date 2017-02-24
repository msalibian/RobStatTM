\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\name{covMcd}
\title{Robust Location and Scatter Estimation via MCD}
\alias{covMcd}
\alias{print.mcd}
\alias{.MCDcons}
\alias{.MCDcnp2}
\alias{.MCDcnp2.rew}
\alias{.MCDsingularityMsg}
\alias{.wgtFUN.covMcd}
%
\concept{High breakdown point}
\description{
  Compute the Minimum Covariance Determinant (MCD) estimator,
  a robust multivariate location and scale estimate with a high
  breakdown point, via the \sQuote{Fast MCD} or \sQuote{Deterministic
    MCD} (\dQuote{DetMcd}) algorithm.
}
\usage{
covMcd(x, cor = FALSE, raw.only = FALSE,
       alpha =, nsamp =, nmini =, kmini =,
       scalefn =, maxcsteps =,
       initHsets = NULL, save.hsets = FALSE, names = TRUE, % full.h = save.hsets,
       seed =, tolSolve =, trace =,
       use.correction =, wgtFUN =, control = rrcov.control())
}
\arguments{
  \item{x}{a matrix or data frame.}
  \item{cor}{should the returned result include a correlation matrix?
    Default is \code{cor = FALSE}.}
  \item{raw.only}{should only the \dQuote{raw} estimate be returned,
    i.e., no (re)weighting step be performed; default is false.}
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the determinant is minimized; roughly \code{alpha*n},
    (see \sQuote{Details} below)
    observations are used for computing the determinant.  Allowed values
    are between 0.5 and 1 and the default is 0.5.}
  \item{nsamp}{number of subsets used for initial estimates or \code{"best"},
    \code{"exact"}, or \code{"deterministic"}.  Default is \code{nsamp = 500}.
    For \code{nsamp = "best"} exhaustive enumeration is done, as long as
    the number of trials does not exceed 100'000 (\code{= nLarge}).
    For \code{"exact"}, exhaustive enumeration will be attempted however
    many samples are needed.  In this case a warning message may be
    displayed saying that the computation can take a very long time.

    For \code{"deterministic"}, the \emph{deterministic} MCD is computed; as
    proposed by Hubert et al. (2012) it starts from the \eqn{h} most
    central observations of \emph{six} (deterministic) estimators.
  }
  \item{nmini, kmini}{for \eqn{n \ge 2 \times n_0}{n >= 2 n_0},
    \eqn{n_0 := \code{nmini}}, the algorithm splits the data into
    maximally \code{kmini} (by default 5) subsets, of size
    approximately, but at least \code{nmini}.  When \code{nmini*kmini < n},
    the initial search uses only a \emph{subsample} of size \code{nmini*kmini}.
    %% FIXME? -- more accurately ==> ../src/rffastmcd.f
    The original algorithm had \code{nmini = 300} and \code{kmini = 5}
    hard coded.}
  \item{scalefn}{for the deterministic MCD: \code{\link{function}} to
    compute a robust scale estimate or character string specifying a
    rule determining such a function.  The default, currently
    \code{"hrv2012"}, uses the recommendation of Hubert, Rousseeuw and
    Verdonck (2012) who recommend \code{\link{Qn}}
    for \eqn{n < 1000} and \code{\link{scaleTau2}} for larger n.  Alternatively,
    \code{scalefn = "v2014"}, uses that rule with cutoff \eqn{n = 5000}.}
  \item{maxcsteps}{maximal number of concentration steps in the
    deterministic MCD; should not be reached.}
  \item{initHsets}{NULL or a \eqn{K x h} integer matrix of initial
    subsets of observations of size \eqn{h} (specified by the indices in
    \code{1:n}).}
  \item{save.hsets}{(for deterministic MCD) logical indicating if the
    initial subsets should be returned as \code{initHsets}.}
  \item{names}{logical; if true (as by default), several parts of the
    result have a \code{\link{names}} or \code{\link{dimnames}}
    respectively, derived from data matrix \code{x}.}
  \item{seed}{initial seed for random generator, like
    \code{\link{.Random.seed}}, see \code{\link{rrcov.control}}.}
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in \code{\link{mahalanobis}}.}
  \item{trace}{logical (or integer) indicating if intermediate results
    should be printed; defaults to \code{FALSE}; values \eqn{\ge 2}{>= 2}
    also produce print from the internal (Fortran) code.}
  \item{use.correction}{ whether to use finite sample correction
    factors; defaults to \code{TRUE}.}
  \item{wgtFUN}{a character string or \code{\link{function}}, specifying
    how the weights for the reweighting step should be computed.  Up to
    April 2013, the only option has been the original proposal in (1999),
    now specified by \code{wgtFUN = "01.original"} (or via
    \code{control}).  Since \pkg{robustbase} version 0.92-3, Dec.2014,
    other predefined string options are available, though experimental,
    see the experimental \code{.wgtFUN.covMcd} object.}
  \item{control}{a list with estimation options - this includes those
    above provided in the function specification, see
    \code{\link{rrcov.control}} for the defaults.  If \code{control} is
    supplied, the parameters from it will be used.  If parameters are
    passed also in the invocation statement, they will override the
    corresponding elements of the control object.}
}
\details{
  The minimum covariance determinant estimator of location and scatter
  implemented in \code{covMcd()} is similar to \R function
  \code{\link[MASS]{cov.mcd}()} in \pkg{MASS}.  The MCD method looks for
  the \eqn{h (> n/2)} (\eqn{h = h(\alpha,n,p) =}
  \code{\link{h.alpha.n}(alpha,n,p)}) observations (out of \eqn{n})
  whose classical covariance matrix has the lowest possible determinant.

  The raw MCD estimate of location is then the average of these \eqn{h} points,
  whereas the raw MCD estimate of scatter is their covariance matrix,
  multiplied by a consistency factor (\code{.MCDcons(p, h/n)}) and (if
  \code{use.correction} is true) a finite sample correction factor
  (\code{.MCDcnp2(p, n, alpha)}), to make it consistent at the
  normal model and unbiased at small samples.  Both rescaling factors
  (consistency and finite sample) are returned in the length-2 vector
  \code{raw.cnp2}.

  The implementation of \code{covMcd} uses the Fast MCD algorithm of
  Rousseeuw and Van Driessen (1999) to approximate the minimum
  covariance determinant estimator.

  Based on these raw MCD estimates, (unless argument \code{raw.only} is
  true), a reweighting step is performed, i.e., \code{V <- \link{cov.wt}(x,w)},
  where \code{w} are weights determined by \dQuote{outlyingness} with
  respect to the scaled raw MCD.  Again, a consistency factor and
  (if \code{use.correction} is true) a finite sample correction factor
  (\code{.MCDcnp2.rew(p, n, alpha)}) are applied.
  The reweighted covariance is typically considerably more efficient
  than the raw one, see Pison et al. (2002).

  The two rescaling factors for the reweighted estimates are returned in
  \code{cnp2}.  Details for the computation of the finite sample
  correction factors can be found in Pison et al. (2002).
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}, based on
  work written for S-plus by Peter Rousseeuw and Katrien van Driessen
  from University of Antwerp.% no E-mails for spam-protection

  Visibility of (formerly internal) tuning parameters, notably
  \code{wgtFUN()}: Martin Maechler
}
\value{
  An object of class \code{"mcd"} which is basically a
  \code{\link{list}} with components
  \item{center}{the final estimate of location.}
  \item{cov}{the final estimate of scatter.}
  \item{cor}{the (final) estimate of the correlation matrix (only if
    \code{cor = TRUE}).}
  \item{crit}{the value of the criterion, i.e., the logarithm of the
    determinant.  Previous to Nov.2014, it contained the determinant
    itself which can under- or overflow relatively easily.}
  \item{best}{the best subset found and used for computing the raw
    estimates, with \code{\link{length}(best) == quan =
      \link{h.alpha.n}(alpha,n,p)}.}
  \item{mah}{mahalanobis distances of the observations using the final
    estimate of the location and scatter.}
  \item{mcd.wt}{weights of the observations using the final estimate of
    the location and scatter.}
  \item{cnp2}{a vector of length two containing the consistency
    correction factor and the finite sample correction factor of
    the final estimate of the covariance matrix.}
  \item{raw.center}{the raw (not reweighted) estimate of location.}
  \item{raw.cov}{the raw (not reweighted) estimate of scatter.}
  \item{raw.mah}{mahalanobis distances of the observations based on the
    raw estimate of the location and scatter.}
  \item{raw.weights}{weights of the observations based on the raw
    estimate of the location and scatter.}
  \item{raw.cnp2}{a vector of length two containing the consistency
    correction factor and the finite sample correction factor of
    the raw estimate of the covariance matrix.}
  \item{X}{the input data as numeric matrix, without \code{\link{NA}}s.}
  \item{n.obs}{total number of observations.}
  \item{alpha}{the size of the subsets over which the determinant is
    minimized (the default is \eqn{(n+p+1)/2}).}
  \item{quan}{the number of observations, \eqn{h}, on which the MCD is
    based.  If \code{quan} equals \code{n.obs}, the MCD is the classical
    covariance matrix.}
  \item{method}{character string naming the method (Minimum Covariance
    Determinant), starting with \code{"Deterministic"} when
    \code{nsamp="deterministic"}.}
  \item{iBest}{(for the deterministic MCD) contains indices from 1:6
    denoting which of the (six) initial subsets lead to the best set found.}
  \item{n.csteps}{(for the deterministic MCD) for each of the initial
    subsets, the number of C-steps executed till convergence.}
  \item{call}{the call used (see \code{\link{match.call}}).}
}
\references{
  Rousseeuw, P. J. and Leroy, A. M. (1987)
  \emph{Robust Regression and Outlier Detection.} Wiley.

  Rousseeuw, P. J. and van Driessen, K. (1999)
  A fast algorithm for the minimum covariance determinant estimator.
  \emph{Technometrics} \bold{41}, 212--223.

  Pison, G., Van Aelst, S., and Willems, G. (2002)
  Small Sample Corrections for LTS and MCD,
  \emph{Metrika} \bold{55}, 111--123.% ~/save/papers/robust-diverse/Pison_VanAelst_Willems.pdf

  Hubert, M., Rousseeuw, P. J. and Verdonck, T. (2012)
  A deterministic algorithm for robust location and scatter.
  Journal of Computational and Graphical Statistics \bold{21}, 618--637.
}
\seealso{
  \code{\link[MASS]{cov.mcd}} from package \CRANpkg{MASS};
  \code{\link{covOGK}} as cheaper alternative for larger dimensions.

  \code{\link[robustX]{BACON}} and \code{\link[robustX]{covNNC}},
  from package \CRANpkg{robustX};
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
set.seed(17)
(cH <- covMcd(hbk.x))
cH0 <- covMcd(hbk.x, nsamp = "deterministic")
with(cH0, stopifnot(quan == 39,
     iBest == c(1:4,6), # 5 out of 6 gave the same
     identical(raw.weights, mcd.wt),
     identical(which(mcd.wt == 0), 1:14), all.equal(crit, -1.045500594135)))

## the following three statements are equivalent
c1 <- covMcd(hbk.x, alpha = 0.75)
c2 <- covMcd(hbk.x, control = rrcov.control(alpha = 0.75))
## direct specification overrides control one:
c3 <- covMcd(hbk.x, alpha = 0.75,
             control = rrcov.control(alpha=0.95))
c1

## Martin's smooth reweighting:

## List of experimental pre-specified wgtFUN() creators:
## Cutoffs may depend on  (n, p, control$beta) :
str(.wgtFUN.covMcd)

cMM <- covMcd(hbk.x, wgtFUN = "sm1.adaptive")

ina <- which(names(cH) == "call")
all.equal(cMM[-ina], cH[-ina]) # *some* differences, not huge (same 'best'):
stopifnot(all.equal(cMM[-ina], cH[-ina], tol = 0.2))
}
\keyword{robust}
\keyword{multivariate}
