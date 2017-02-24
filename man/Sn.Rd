\name{Sn}
\alias{Sn}
\alias{s_Sn}
%
\title{Robust Location-Free Scale Estimate More Efficient than MAD}
\description{
  Compute the robust scale estimator \eqn{S_n}{Sn}, an efficient
  alternative to the MAD.
}
\usage{
Sn(x, constant = 1.1926, finite.corr = missing(constant))

s_Sn(x, mu.too = FALSE, \dots)
}
\arguments{
  \item{x}{numeric vector of observations.}
  \item{constant}{number by which the result is multiplied; the default
    achieves consisteny for normally distributed data.}
  \item{finite.corr}{logical indicating if the finite sample bias
    correction factor should be applied.  Default to \code{TRUE} unless
    \code{constant} is specified.}
  \item{mu.too}{logical indicating if the \code{\link[stats]{median}(x)} should
    also be returned for \code{s_Sn()}.}
  \item{\dots}{potentially further arguments for \code{s_Sn()} passed to
    \code{Sn()}.}
}
\value{
  \code{Sn()} returns a number, the \eqn{S_n}{Sn} robust scale estimator, scaled to be
  consistent for \eqn{\sigma^2} and i.i.d. Gaussian observatsions,
  optionally bias corrected for finite samples.

  \code{s_Sn(x, mu.too=TRUE)} returns a length-2 vector with location
  (\eqn{\mu}) and scale; this is typically only useful for
  \code{\link{covOGK}(*, sigmamu = s_Sn)}.

}
\details{
  ............  FIXME ........
}
\references{
  Rousseeuw, P.J. and Croux, C. (1993)
  Alternatives to the Median Absolute Deviation,
  \emph{Journal of the American Statistical Association} \bold{88}, 1273--1283.
}
\seealso{\code{\link[stats]{mad}} for the \sQuote{most robust} but much
  less efficient scale estimator;
  \code{\link{Qn}} for a similar more efficient but slower alternative;
  \code{\link{scaleTau2}}.
}
\author{Original Fortran code:
  Christophe Croux and Peter Rousseeuw \email{rousse@wins.uia.ac.be}.
  \cr
  Port to C and R: Martin Maechler, \email{maechler@R-project.org}
}
\examples{
x <- c(1:10, 100+1:9)# 9 outliers out of 19
Sn(x)
Sn(x, c=1)# 9
Sn(x[1:18], c=1)# 9
set.seed(153)
x <- sort(c(rnorm(80), rt(20, df = 1)))
s_Sn(x, mu.too=TRUE)
}
\keyword{robust}
\keyword{univar}
