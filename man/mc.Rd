\name{mc}
\alias{mc}
\title{Medcouple, a Robust Measure of Skewness}
\description{
  Compute the \sQuote{medcouple}, a \emph{robust} concept and estimator
  of skewness.  The medcouple is defined as a scaled median difference
  of the left and right half of distribution, and hence \emph{not} based
  on the third moment as the classical skewness.
}
\usage{
mc(x, na.rm = FALSE, doReflect = (length(x) <= 100),
   eps1 = .Machine$double.eps, eps2 = .Machine$double.xmin,
   maxit = 100, trace.lev = 0, full.result = FALSE)
}
\arguments{
  \item{x}{a numeric vector}
  \item{na.rm}{logical indicating how missing values (\code{\link{NA}}s)
    should be dealt with.}
  \item{doReflect}{logical indicating if the internal MC should also be
    computed on the \emph{reflected} sample \code{-x}, with final result
    \code{(mc.(x) - mc.(-x))/2}.  This makes sense since the internal
    MC, \code{mc.()} computes the himedian() which can differ slightly from
    the median.}%% only whenever sum(x <= med) * sum(x >= med) is even
  \item{eps1,eps2}{tolerance in the algorithm; only change with care!}
  \item{maxit}{maximul number of iterations; typically a few should be
    sufficient.}
  \item{trace.lev}{integer specifying how much diagnostic output the
    algorithm (in C) should produce.  No output by default, most output
    for \code{trace.lev = 5}.}
  \item{full.result}{logical indicating if the full return values (from
    C) should be returned as a list via \code{attr(*, "mcComp")}.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  a number between -1 and 1, which is the medcouple, \eqn{MC(x)}.
  For \code{r <- mc(x, full.result = TRUE, ....)}, then
  \code{attr(r, "mcComp")} is a list with components
  \item{medc}{the medcouple  \eqn{mc.(x)}.}
  \item{medc2}{the medcouple \eqn{mc.(-x)} if \code{doReflect=TRUE}.}
  \item{eps}{tolerances used.}
  \item{iter,iter2}{number of iterations used.}
  \item{converged,converged2}{logical specifying \dQuote{convergence}.}
}
\references{
  Guy Brys, Mia Hubert and Anja Struyf (2004)
  A Robust Measure of Skewness;
  \emph{JCGS} \bold{13} (4), 996--1017.

  Hubert, M. and Vandervieren, E. (2008).
  An adjusted boxplot for skewed distributions,
  \emph{Computational Statistics and Data Analysis} \bold{52}, 5186--5201.
}
\author{Guy Brys; modifications by Tobias Verbeke and bug fixes and
  extensions by Manuel Koller and Martin Maechler.
}
\seealso{\code{\link{Qn}} for a robust measure of scale (aka
  \dQuote{dispersion}), ....
}
\examples{
mc(1:5) # 0 for a symmetric sample

x1 <- c(1, 2, 7, 9, 10)
mc(x1) # = -1/3

data(cushny)
mc(cushny) # 0.125

stopifnot(mc(c(-20, -5, -2:2, 5, 20)) == 0,
          mc(x1, doReflect=FALSE) ==  -mc(-x1, doReflect=FALSE),
          all.equal(mc(x1, doReflect=FALSE), -1/3, tolerance = 1e-12))
}
\keyword{robust}
\keyword{univar}
