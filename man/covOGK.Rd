\name{covOGK}
\alias{covOGK}
\alias{covGK}
\alias{s_mad}
\alias{s_IQR}
\alias{hard.rejection}
%
\title{Orthogonalized Gnanadesikan-Kettenring (OGK) Covariance Matrix Estimation}
\description{
  Computes the orthogonalized pairwise covariance matrix estimate described in
  in Maronna and Zamar (2002).  The pairwise proposal goes back to
  Gnanadesikan and Kettenring (1972).
}
\usage{
covOGK(X, n.iter = 2, sigmamu, rcov = covGK, weight.fn = hard.rejection,
       keep.data = FALSE, \dots)

covGK (x, y, scalefn = scaleTau2, \dots)
s_mad(x, mu.too = FALSE, na.rm = FALSE)
s_IQR(x, mu.too = FALSE, na.rm = FALSE)
}
\arguments{
  \item{X}{data in something that can be coerced into a numeric matrix.}
  \item{n.iter}{number of orthogonalization iterations.  Usually 1 or 2;
    values greater than 2 are unlikely to have any significant effect on
    the estimate (other than increasing the computing time).}
  \item{sigmamu, scalefn}{a function that computes univariate robust
    location and scale estimates.  By default it should return a single
    numeric value containing the robust scale (standard deviation)
    estimate.  When \code{mu.too} is true, \code{sigmamu()} should
    return a numeric vector of length 2 containing robust location and
    scale estimates.  See \code{\link{scaleTau2}}, \code{\link{s_Qn}},
    \code{\link{s_Sn}}, \code{s_mad} or \code{s_IQR} for examples to be
    used as \code{sigmamu} argument.}
  \item{rcov}{function that computes a robust covariance estimate
    between two vectors.  The default, Gnanadesikan-Kettenring's
    \code{covGK}, is simply \eqn{(s^2(X+Y) - s^2(X-Y))/4} where
    \eqn{s()} is the scale estimate \code{sigmamu()}.}
  \item{weight.fn}{a function of the robust distances and the number of
    variables \eqn{p} to compute the weights used in the reweighting step.}
  \item{keep.data}{logical indicating if the (untransformed) data matrix
    \code{X} should be kept as part of the result.}
  \item{\dots}{additional arguments; for \code{covOGK} to be passed to
    \code{sigmamu()} and \code{weight.fn()}; for \code{covGK} passed to \code{scalefn}.}
%%% covGK():
  \item{x,y}{numeric vectors of the same length, the covariance of which
    is sought in \code{covGK} (or the scale, in \code{s_mad} or
    \code{s_IQR}).}
%%% s_mad(), s_IQR():
  \item{mu.too}{logical indicating if both location and scale should be
    returned or just the scale (when \code{mu.too=FALSE} as by default).}
  \item{na.rm}{if \code{TRUE} then \code{\link{NA}} values are stripped
    from \code{x} before computation takes place.}
}
\details{
  Typical default values for the \emph{function} arguments
  \code{sigmamu}, \code{rcov}, and \code{weight.fn}, are
  available as well, see the \emph{Examples} below,
  \bold{but} their names and calling sequences are
  still subject to discussion and may be changed in the future.

  The current default, \code{weight.fn = hard.rejection} corresponds to
  the proposition in the litterature, but Martin Maechler strongly
  believes that the hard threshold currently in use is too arbitrary,
  and further that \emph{soft} thresholding should be used instead, anyway.
}
\value{
  \code{covOGK()} currently returns a list with components
  \item{center}{robust location: numeric vector of length \eqn{p}.}
  \item{cov}{robust covariance matrix estimate: \eqn{p\times p}{p x p}
    matrix.}
  \item{wcenter, wcov}{re-\bold{w}eighted versions of \code{center} and
    \code{cov}.}
  \item{weights}{the robustness weights used.}
  \item{distances}{the mahalanobis distances computed using
    \code{center} and \code{cov}.}
  \dots\dots
  \cr
  \bold{but note that this might be radically changed to returning an
    S4 classed object!}

  \code{covGK()} is a trivial 1-line function returning the covariance
  estimate
  \deqn{\hat c(x,y) = \left(\hat \sigma(x+y)^2 - \hat \sigma(x-y)^2 \right)/4,%
        }{  c^(x,y) = [s^(x+y)^2 - s^(x-y)^2]/4,}%
  where \eqn{\hat \sigma(u)}{s^(u)} is the scale estimate of \eqn{u}
  specified by \code{scalefn}.

  \code{s_mad()}, and \code{s_IQR()} return the
  scale estimates \code{\link[stats]{mad}} or \code{\link[stats]{IQR}}
  respectively, where the \code{s_*} functions return a length-2 vector
  (mu, sig) when \code{mu.too = TRUE}, see also \code{\link{scaleTau2}}.
}
\references{
  Maronna, R.A. and Zamar, R.H. (2002)
  Robust estimates of location and dispersion of high-dimensional datasets;
  \emph{Technometrics} \bold{44}(4), 307--317.

  Gnanadesikan, R. and John R. Kettenring (1972)
  Robust estimates, residuals, and outlier detection with multiresponse data.
  \emph{Biometrics} \bold{28}, 81--124.
}
\author{Kjell Konis \email{konis@stats.ox.ac.uk}, with modifications by
    Martin Maechler.}
\seealso{\code{\link{scaleTau2}},
  \code{\link{covMcd}}, \code{\link[MASS]{cov.rob}}.
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])

cO1 <- covOGK(hbk.x, sigmamu = scaleTau2)
cO2 <- covOGK(hbk.x, sigmamu = s_Qn)
cO3 <- covOGK(hbk.x, sigmamu = s_Sn)
cO4 <- covOGK(hbk.x, sigmamu = s_mad)
cO5 <- covOGK(hbk.x, sigmamu = s_IQR)

%% FIXME: Add time comparison, here or in "vignette", "demo", "...

data(toxicity)
cO1tox <- covOGK(toxicity, sigmamu = scaleTau2)
cO2tox <- covOGK(toxicity, sigmamu = s_Qn)

## nice formatting of correlation matrices:
as.dist(round(cov2cor(cO1tox$cov), 2))
as.dist(round(cov2cor(cO2tox$cov), 2))

## "graphical"
symnum(cov2cor(cO1tox$cov))
symnum(cov2cor(cO2tox$cov), legend=FALSE)
}
\keyword{robust}
\keyword{multivariate}
