\name{rankMM}
\title{Simple Matrix Rank}
\alias{rankMM}
\description{
  Compute the rank of a matrix \code{A} in simple way, based on the SVD,
  \code{\link{svd}()}, and \dQuote{the same as Matlab}.
}
\usage{
rankMM(A, tol = NULL, sv = svd(A, 0, 0)$d)
}
\arguments{
  \item{A}{a numerical matrix, maybe non-square.  When \code{sv} is
    specified, only \code{dim(A)} is made use of.}
  \item{tol}{numerical tolerance (compared to singular values).  By
    default, when \code{NULL}, the tolerance is determined from the
    maximal value of \code{sv} and the computer epsilon.}
  \item{sv}{vector of \emph{non-increasing} singular values of \code{A}, (to be
    passed if already known).}
}
\seealso{
  There are more sophisticated proposals for computing the rank of a
  matrix; for a couple of those, see \code{\link[Matrix]{rankMatrix}} in the
  \pkg{Matrix} package.
}
\value{
  an integer from the set \code{0:min(dim(A))}.
}
\author{
  Martin Maechler, Date: 7 Apr 2007
}
\examples{
rankMM # - note the simple function definition

hilbert <- function(n) { i <- seq_len(n); 1/outer(i - 1L, i, "+") }
hilbert(4)
H12 <- hilbert(12)
rankMM(H12)        # 11 - numerically more realistic
rankMM(H12, tol=0) # -> 12
## explanation :
round(log10(svd(H12, 0,0)$d), 1)
}
\keyword{algebra}
\keyword{array}

