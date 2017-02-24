\name{SiegelsEx}
\alias{SiegelsEx}
\docType{data}
\title{Siegel's Exact Fit Example Data}
\description{
  A small counterexample data set devised by Andrew Siegel.
  Six (out of nine) data points lie on the line \eqn{y = 0} such that
  some robust regression estimators exhibit the \dQuote{\emph{exact fit}}
  property.
}
\usage{data(SiegelsEx)}
\format{
  A data frame with 9 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\source{
  Emerson and Hoaglin (1983, p.139)
}
\references{
  Peter J. Rousseeuw and Annick M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection} Wiley, p.60--61
}
\examples{
data(SiegelsEx)
plot(SiegelsEx, main = "Siegel's example for 'exact fit'")
abline(          lm(y ~ x, data = SiegelsEx))
abline(MASS::lqs(y ~ x, data = SiegelsEx, method = "lms"), col = 2)
legend("topright", leg = c("lm", "LMS"), col=1:2, lwd=1, inset = 1/20)
}
\keyword{datasets}
