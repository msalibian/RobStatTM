\name{sigma}
\title{Extract 'Sigma' - Standard Deviation of Errors for Robust Models}
\alias{sigma}
\alias{sigma.lmrob}
\description{
  Extract the estimated standard deviation of the errors, the
  \dQuote{residual standard deviation} (misnomed also
  \dQuote{residual standard error}) from a fitted model.
}
\usage{
% Declare S3 method as the generic is no longer in lme4 for R (>= 2015-06-01):
\S3method{sigma}{lmrob}(object, \dots)
}
\arguments{
  \item{object}{a fitted model.}

  \item{\dots}{additional, optional arguments.  (None are
  used in our methods)}
}
\value{
  the residual standard error as a scalar
}
\details{
  For \R \code{<= 3.2.x}, we provide an (S3) generic function (as e.g.,
  package \pkg{lme4}) and methods for \code{\link{lmrob}},
  \code{\link{nlrob}}, and \code{\link{nls}}.

  From \R \code{>= 3.3.0}, we provide methods for our
  \code{\link{lmrob}} and \code{\link{nlrob}} models.
}
\examples{
      m.cl <-   lm (Y ~ ., data=coleman)
if(getRversion() >= "3.3.0")  sigma(m.cl)  else  summary(m.cl)$sigma
sigma( m1  <- lmrob(Y ~ ., data=coleman) )
sigma( m2  <- lmrob(Y ~ ., data=coleman, setting = "KS2014") )
}
\keyword{models}
