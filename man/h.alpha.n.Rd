\name{h.alpha.n}
\alias{h.alpha.n}
\title{Compute h, the subsample size for MCD and LTS}
\description{
 Compute h(alpha) which is the size of the subsamples to be used
 for MCD and LTS.  Given \eqn{\alpha = alpha}{alpha}, \eqn{n} and
 \eqn{p}, \eqn{h} is an \emph{integer}, \eqn{h \approx \alpha n}{h ~=
   alpha*n}, where the exact formula also depends on \eqn{p}.

 For \eqn{\alpha = 1/2}, \code{h == floor(n+p+1)/2}; for the general
 case, it's simply
 \code{n2 <- (n+p+1) \%/\% 2; floor(2*n2 - n + 2*(n-n2)*alpha)}.
}
\usage{
h.alpha.n(alpha, n, p)
}
\arguments{
  \item{alpha}{fraction, numeric (vector) in [0.5, 1], see, e.g.,
    \code{\link{covMcd}}.}
  \item{n}{integer (valued vector), the sample size.}
  \item{p}{integer (valued vector), the dimension.}
}
\value{
  numeric vector of \eqn{h(\alpha, n,p)}; when any of the arguments of
  length greater than one, the usual \R arithmetic (recycling) rules are
  used.
}
\seealso{\code{\link{covMcd}} and \code{\link{ltsReg}} which are
  \emph{defined} by \eqn{h = h(\alpha,n,p)} and hence both use
  \code{h.alpha.n}.
}
\examples{
n <- c(10:20,50,100)
p <- 5
## show the simple "alpha = 1/2" case:
cbind(n=n, h= h.alpha.n(1/2, n, p), n2p = floor((n+p+1)/2))

## alpha = 3/4 is recommended by some authors :
n <- c(15, 20, 25, 30, 50, 100)
cbind(n=n, h= h.alpha.n(3/4, n, p = 6))
}
\keyword{arith}
