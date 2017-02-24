\name{wgt.himedian}
\alias{wgt.himedian}
\title{Weighted Hi-Median}
\description{
  Compute the weighted Hi-Median of \code{x}.
}
\usage{
wgt.himedian(x, weights = rep(1, n))
}
\arguments{
  \item{x}{numeric vector}
  \item{weights}{numeric vector of weights; of the same length as \code{x}.}
}
% \value{
%   a number
% }
%\author{Martin Maechler}
\note{
  this is rather a by-product of the code used in \code{\link{Sn}} and
  \code{\link{Qn}}.  We currently plan to replace it with more general
  weighted quantiles.
}
\seealso{\code{\link{median}}; %%.....
  also \code{\link[Hmisc]{wtd.quantile}} from package \pkg{Hmisc}.
}
\examples{
x <- c(1:6, 20)
median(x) ## 4
stopifnot(all.equal(4, wgt.himedian(x)),
          all.equal(6, wgt.himedian(x, c(rep(1,6), 5))))
}
\keyword{univar}
\keyword{robust}
