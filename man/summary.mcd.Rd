\name{summary.mcd}
\alias{summary.mcd}
\alias{print.summary.mcd}
\title{Summary Method for MCD objects}
\usage{
\method{summary}{mcd}(object, \dots)
\method{print}{summary.mcd}(x, digits = max(3, getOption("digits") - 3),
     print.gap = 2, \dots)
}
\arguments{
  \item{object,x}{an object of class \code{"mcd"} (or \code{"summary.mcd"});
    usually, a result of a call to \code{\link{covMcd}}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{print.gap}{number of horizontal spaces between numbers; see also
    \code{\link{print.default}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  \code{\link{summary}} method for class \code{"mcd"}.
}
\details{
  \code{summary.mcd()}, the S3 method, simply returns an (S3) object of
  \code{\link{class} "summary.mcd"} for which there's a
  \code{\link{print}} method:

  \code{print.summary.mcd} prints summary statistics for the weighted covariance
  matrix and location estimates with weights based on MCD estimates.
  While the function \code{\link{print.mcd}} prints only the robust estimates
  of the location and the covariance matrix, \code{print.summary.mcd} will
  print also the correlation matrix (if requested in the call to
  \code{covMcd} with \code{cor=TRUE}), the eigenvalues of the covariance
  or the correlation matrix and the robust (\dQuote{Mahalanobis}) distances.
}
\value{
  \code{summary.mcd} returns an \code{summary.mcd} object, whereas the
  \code{print} methods returns its first argument via
  \code{\link{invisible}}, as all \code{print} methods do.
}
\seealso{
  \code{\link{covMcd}}, \code{\link{summary}}
}
\examples{
data(Animals, package = "MASS")
brain <- Animals[c(1:24, 26:25, 27:28),]
lbrain <- log(brain)
summary(cLB <- covMcd(lbrain))
}
\keyword{multivariate}
\keyword{robust}
