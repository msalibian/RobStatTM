\name{print.lmrob}
\alias{print.lmrob}
\title{Print Method for Objects of Class "lmrob"}
\description{
  Print method for elements of class \code{"lmrob"}.
}
\usage{
\method{print}{lmrob}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{x}{an \R object of class \code{lmrob}, typically created by
    \code{\link{lmrob}}.}
  \item{digits}{number of digits for printing, see \code{digits} in
    \code{\link{options}}.}
  \item{\dots}{potentially more arguments passed to methods.}
}
\seealso{\code{\link{lmrob}}, \code{\link{summary.lmrob}},
    \code{\link{print}} and \code{\link{summary}}.
}
\examples{
data(coleman)
( m1 <- lmrob(Y ~ ., data=coleman) ) # -> print.lmrob() method
}
\keyword{robust}
\keyword{regression}
