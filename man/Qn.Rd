\name{Qn}
\alias{Qn}
\alias{Qn.old}
\alias{s_Qn}
%
\title{Robust Location-Free Scale Estimate More Efficient than MAD}
\description{
  Compute the robust scale estimator \eqn{Q_n}{Qn}, an efficient
  alternative to the MAD.

  See the references for more.
}
\usage{
Qn(x, constant = 2.21914, finite.corr = missing(constant))

s_Qn(x, mu.too = FALSE, \dots)
}
\arguments{
  \item{x}{numeric vector of observations.}
  \item{constant}{number by which the result is multiplied; the default
    achieves consistency for normally distributed data.  Note that until
    Nov. 2010, \dQuote{thanks} to a typo in the very first papers, a slightly
    wrong default constant, 2.2219, was used instead of the correct one
    which is equal to \code{1 / (sqrt(2) * qnorm(5/8))} (as mentioned
    already on p.1277, after (3.7) in Rousseeuw and Croux (1993)).

    If you need the old slightly off version for historical
    reproducibility, you can use \code{Qn.old()}.

    Note that the relative difference is only about 1 in 1000, and that
    the correction should not affect the finite sample corrections for
    \eqn{n \le 9}{n <= 9}.
  }
  \item{finite.corr}{logical indicating if the finite sample bias
    correction factor should be applied.  Defaults to \code{TRUE} unless
    \code{constant} is specified.}
  \item{mu.too}{logical indicating if the \code{\link[stats]{median}(x)} should
    also be returned for \code{s_Qn()}.}
  \item{\dots}{potentially further arguments for \code{s_Qn()} passed to
    \code{Qn()}.}
}
\value{
  \code{Qn()} returns a number, the \eqn{Q_n}{Qn} robust scale
  estimator, scaled to be consistent for \eqn{\sigma^2} and
  i.i.d. Gaussian observatsions, optionally bias corrected for finite
  samples.

  \code{s_Qn(x, mu.too=TRUE)} returns a length-2 vector with location
  (\eqn{\mu}) and scale; this is typically only useful for
  \code{\link{covOGK}(*, sigmamu = s_Qn)}.
}
\details{
  As the (default, consistency) constant needed to be corrected,
  the finite sample correction has been based on a much more extensive
  simulation, and on a 3rd or 4th degree polynomial model in \eqn{1/n}
  for odd or even n, respectively.
}
\references{
  Rousseeuw, P.J. and Croux, C. (1993)
  Alternatives to the Median Absolute Deviation,
  \emph{Journal of the American Statistical Association} \bold{88}, 1273--1283.
  % MM: ~/save/papers/Rousseeuw/93/R+Croux_MAD_Sn_Qn.pdf

  Christophe Croux and Peter J. Rousseeuw (1992)
  Time-Efficient Algorithms for Two Highly Robust Estimators of Scale,
  \emph{Computational Statistics, Vol. 1}, ed. Dodge and Whittaker,
  Physica-Verlag Heidelberg, 411--428; available via Springer Link.
  % MM: ~/save/papers/robust-diverse/Croux-Rousseeuw-Timeff_Scale_1992.pdf
  %% no longer \url{http://win-www.uia.ac.be/u/statis/abstract/Timeff92.htm}.

  About the typo in the \code{constant}:\cr
  Christophe Croux (2010)
  Private e-mail, Fri Jul 16, w/ Subject
  \emph{Re: Slight inaccuracy of Qn implementation \dots\dots}.
}
\seealso{\code{\link[stats]{mad}} for the \sQuote{most robust} but much less efficient
  scale estimator; \code{\link{Sn}} for a similar faster but less
  efficient alternative.  Finally, \code{\link{scaleTau2}} which some
  consider \dQuote{uniformly} better than Qn or competitors.
}
\author{Original Fortran code:
  Christophe Croux and Peter Rousseeuw \email{rousse@wins.uia.ac.be}.
  \cr
  Port to C and R: Martin Maechler, \email{maechler@R-project.org}
}
\examples{
set.seed(153)
x <- sort(c(rnorm(80), rt(20, df = 1)))
s_Qn(x, mu.too = TRUE)
Qn(x, finite.corr = FALSE)
}
\keyword{robust}
\keyword{univar}
