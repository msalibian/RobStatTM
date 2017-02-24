\name{classPC}
\title{Compute Classical Principal Components via SVD or Eigen}
\alias{classPC}
\concept{PCA}
\description{
  Compute classical principal components (PC) via SVD (\code{\link{svd}}
  or eigenvalue decomposition (\code{\link{eigen}}) with non-trivial
  rank determination.
}
\usage{
classPC(x, scale = FALSE, center = TRUE, signflip = TRUE,
        via.svd = n > p, scores = FALSE)
}
\arguments{
  \item{x}{a numeric \code{\link{matrix}}.}
  \item{scale}{logical indicating if the matrix should be scaled; it is
    mean centered in any case (via
    \code{\link{scale}(*, scale=scale)}c}
  \item{center}{logical or numeric vector for \dQuote{centering} the matrix.}
  \item{signflip}{logical indicating if the sign(.) of the loadings
    should be determined should flipped such that the absolutely largest
    value is always positive.}
  \item{via.svd}{logical indicating if the computation is
    via SVD or Eigen decomposition; the latter makes sense
    typically only for n <= p.}
  \item{scores}{logical indicating}
}
\author{
Valentin Todorov; efficiency tweaks by Martin Maechler
}
\value{
  a \code{\link{list}} with components
  \item{rank}{the (numerical) matrix rank of \code{x}; an integer
    number, say \eqn{k}, from \code{0:min(dim(x))}.  In the \eqn{n > p} case,
    it is \code{\link{rankMM}(x)}.}
  \item{eigenvalues}{the \eqn{k} eigenvalues, in the \eqn{n > p} case,
    proportional to the variances.}
  \item{loadings}{the loadings, a \eqn{p \times k}{p * k} matrix.}
  \item{scores}{if the \code{scores} argument was true, the \eqn{n \times
      k}{n * k} matrix of scores, where \eqn{k} is the \code{rank} above.}
  \item{center}{a numeric \eqn{p}-vector of means, unless the
    \code{center} argument was false.}
  \item{scale}{if the \code{scale} argument was not false, the
    \code{scale} used, a \eqn{p}-vector.}
}
%% \details{
%% }
%% \references{
%% }
\seealso{
  In spirit very similar to \R's standard \code{\link{prcomp}} and
  \code{\link{princomp}}, one of the main differences being how the
  \emph{rank} is determined via a non-trivial tolerance.
}
\examples{
set.seed(17)
x <- matrix(rnorm(120), 10, 12) # n < p {the unusual case}
pcx  <- classPC(x)
(k <- pcx$rank) # = 9  [after centering!]
pc2  <- classPC(x, scores=TRUE)
pcS  <- classPC(x, via.svd=TRUE)
all.equal(pcx, pcS, tol = 1e-8)
## TRUE: eigen() & svd() based PC are close here
pc0 <- classPC(x, center=FALSE, scale=TRUE)
pc0$rank # = 10  here *no* centering (as E[.] = 0)

## Loadings are orthnormal:
zapsmall( crossprod( pcx$loadings ) )

## PC Scores are roughly orthogonal:
S.S <- crossprod(pc2$scores)
print.table(signif(zapsmall(S.S), 3), zero.print=".")
stopifnot(all.equal(pcx$eigenvalues, diag(S.S)/k))

## the usual n > p case :
pc.x <- classPC(t(x))
pc.x$rank # = 10, full rank in the n > p case

cpc1 <- classPC(cbind(1:3)) # 1-D matrix
stopifnot(cpc1$rank == 1,
          all.equal(cpc1$eigenvalues, 1),
          all.equal(cpc1$loadings, 1))
\dontshow{
stopifnot(classPC(x, center=FALSE)$rank == min(dim(x)))
ii <- names(pcx); ii <- ii[ii != "scores"]
stopifnot(all.equal(pcx[ii], pc2[ii], tol=0),
	  all.equal(pcx, pcS, tol=1e-8),
	  length(pc.x$center) == 10, identical(pc0$center, FALSE),
          all.equal(crossprod(pcx $loadings), diag(9)),
          all.equal(crossprod(pc.x$loadings), diag(10)),
          all.equal(colSums(abs(pcx$loadings)),
                    c(2.69035673, 2.78449399, 3.00148438,
                      2.9016688,  2.49400759, 2.90477204,
                      3.01639807, 2.4217181, 2.64665957)),
	  length(pc0$scale) == 12)
}% dont..
}
\keyword{multivariate}

