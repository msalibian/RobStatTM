\name{adjboxStats}
\alias{adjboxStats}
\title{Statistics for Skewness-adjusted Boxplots}
\description{
  Computes the \dQuote{statistics} for producing boxplots adjusted for
  skewed distributions as proposed in Hubert and Vandervieren (2004),
  see \code{\link{adjbox}}.
}
\usage{
adjboxStats(x, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE,
            \dots)
}
\arguments{
  \item{x}{a numeric vector for which adjusted boxplot statistics are computed.}
  \item{coef}{number determining how far \sQuote{whiskers} extend out
    from the box, see \code{\link{boxplot.stats}}.}
  \item{a, b}{scaling factors multiplied by the medcouple
    \code{\link{mc}()} to determine outlyer boundaries; see the references.}
  \item{do.conf,do.out}{logicals; if \code{FALSE}, the \code{conf} or
    \code{out} component respectively will be empty in the result.}
  \item{\dots}{further optional arguments to be passed to
    \code{\link{mc}()}, such as \code{doReflect}.}
}
\details{
  Given the quartiles \eqn{Q_1}{Q1}, \eqn{Q_3}{Q3}, the interquartile
  range \eqn{\Delta Q := Q_3 - Q_1}{IQR := Q3-Q1}, and the medcouple
  \eqn{M :=}\code{mc(x)}, \eqn{c =}\code{coef},
  the \dQuote{fence} is defined,
  for \eqn{M \ge 0} as
  \deqn{[Q_1 - c e^{a \cdot M}\Delta Q, Q_3 + c e^{b \cdot M}\Delta Q],%
  }{[Q1 - c*exp(a * M)*IQR, Q3 + c*exp(b * M)*IQR],}
  and for \eqn{M < 0} as
  \deqn{[Q_1 - c e^{-b \cdot M}\Delta Q, Q_3 + c e^{-a \cdot M}\Delta Q],%
  }{[Q1 - c*exp(-b * M)*IQR, Q3 + c*exp(-a * M)*IQR],}
  and all observations \code{x} outside the fence, the \dQuote{potential
  outliers}, are returned in \code{out}.

  Note that a typo in robustbase version up to 0.7-8,
  for the (rare left-skewed) case where \link{mc}(x) < 0, lead to a
  \dQuote{fence} not wide enough in the upper part, and hence
  \emph{less} outliers there.
}
\value{
  A \code{\link{list}} with the components

  \item{stats}{a vector of length 5, containing the extreme of the lower
    whisker, the lower hinge, the median, the upper hinge and the extreme of
    the upper whisker.}
  \item{n}{the number of observations}
  \item{conf}{the lower and upper extremes of the \sQuote{notch}
    (\code{if(do.conf)}). See \code{\link{boxplot.stats}}.}
  \item{fence}{length 2 vector of interval boundaries which
    define the non-outliers, and hence the whiskers of the plot.}
  \item{out}{the values of any data points which lie beyond the fence,
    and hence beyond the extremes of the whiskers.}
}
\author{R Core Development Team (\code{\link{boxplot.stats}}); adapted
  by Tobias Verbeke and Martin Maechler.}
\note{The code only slightly modifies the code of \R's
  \code{\link{boxplot.stats}}.
}
\seealso{\code{\link{adjbox}()}, also for references,
  the function which mainly uses this one;
  further \code{\link{boxplot.stats}}.
}
\examples{
data(condroz)
adjboxStats(ccA <- condroz[,"Ca"])
adjboxStats(ccA, doReflect = TRUE)# small difference in fence

## Test reflection invariance [was not ok, up to and including robustbase_0.7-8]
a1 <- adjboxStats( ccA, doReflect = TRUE)
a2 <- adjboxStats(-ccA, doReflect = TRUE)

nm1 <- c("stats", "conf", "fence")
stopifnot(all.equal(       a1[nm1],
                    lapply(a2[nm1], function(u) rev(-u))),
          all.equal(a1[["out"]], -a2[["out"]]))
}
\keyword{robust}
\keyword{univar}
