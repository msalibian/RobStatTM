\name{colMedians}
\title{Fast Row or Column-wise Medians of a Matrix}
\alias{colMedians}
\alias{rowMedians}
\description{
  Calculates the median for each row (column) of a matrix \code{x}.
  This is the same as but more efficient than \code{apply(x, MM, median)}
  for MM=2 or MM=1, respectively.
}
\usage{
colMedians(x, na.rm = FALSE, hasNA = TRUE, keep.names=TRUE)
rowMedians(x, na.rm = FALSE, hasNA = TRUE, keep.names=TRUE)
}
\arguments{
 \item{x}{a \code{\link{numeric}} \eqn{n \times p}{n x p} \code{\link{matrix}}.}
 \item{na.rm}{if \code{\link{TRUE}}, \code{\link{NA}}s are excluded
   first, otherwise not.}
 \item{hasNA}{logical indicating if \code{x} may contain \code{\link{NA}}s.
   If set to \code{FALSE}, no internal NA handling is performed which
   typically is faster.}
\item{keep.names}{logical indicating if row or column names of \code{x}
  should become \code{\link{names}} of the result - as is the case for
  \code{\link{apply}(x, MM, median)}.}
}
\value{
  a \code{\link{numeric}} vector of length \eqn{n} or \eqn{p}, respectively.
}
\section{Missing values}{
  Missing values are excluded before calculating the medians
  \emph{unless} \code{hasNA} is false.  Note that \code{na.rm} has no
  effect and is automatically false when \code{hasNA} is false, i.e.,
  internally, before computations start, the following is executed:
\preformatted{if (!hasNA)        ## If there are no NAs, don't try to remove them
     narm <- FALSE}
}
\details{
  The implementation of \code{rowMedians()} and \code{colMedians()}
  is optimized for both speed and memory.
  To avoid coercing to \code{\link{double}}s (and hence memory allocation), there
  is a special implementation for \code{\link{integer}} matrices.
  That is, if \code{x} is an \code{\link{integer}} \code{\link{matrix}}, then
  \code{rowMedians(as.double(x))} (\code{rowMedians(as.double(x))})
  would require three times the memory of \code{rowMedians(x)}
  (\code{colMedians(x)}), but all this is avoided.
}
\author{Henrik Bengtsson, Harris Jaffee, Martin Maechler}

\seealso{
  See \code{\link{wgt.himedian}()} for a weighted hi-median, and
  \code{\link[matrixStats]{colWeightedMedians}()} etc from package
  \pkg{matrixStats} for \emph{weighted} medians.\cr
  For mean estimates, see \code{rowMeans()} in \code{\link{colSums}}().
}

\examples{
set.seed(1); n <- 234; p <- 543 # n*p = 127'062
x <- matrix(rnorm(n*p), n, p)
x[sample(seq_along(x), size= n*p / 256)] <- NA
R1 <- system.time(r1 <- rowMedians(x, na.rm=TRUE))
C1 <- system.time(y1 <- colMedians(x, na.rm=TRUE))
R2 <- system.time(r2 <- apply(x, 1, median, na.rm=TRUE))
C2 <- system.time(y2 <- apply(x, 2, median, na.rm=TRUE))
R2 / R1 # speedup factor: ~= 4   {platform dependent}
C2 / C1 # speedup factor: ~= 5.8 {platform dependent}
stopifnot(all.equal(y1, y2, tol=1e-15),
          all.equal(r1, r2, tol=1e-15))

(m <- cbind(x1=3, x2=c(4:1, 3:4,4)))
stopifnot(colMedians(m) == 3,
          all.equal(colMeans(m), colMedians(m)),# <- including names !
          all.equal(rowMeans(m), rowMedians(m)))
}
\keyword{array}
\keyword{robust}
\keyword{univar}
