\name{estimethod}
\alias{estimethod}
\title{Extract the Estimation Method 'Estimethod' from a Fitted Model}
\description{
  Extract the estimation method as a \code{\link{character}} string from a fitted model.
}
\usage{
  estimethod(object, ...)
}
\arguments{
  \item{object}{a fitted model.}
  \item{...}{additional, optional arguments.  (None are
  used in our methods)}
}
\value{
  a \code{\link{character}} string, the estimation method used.
}
\details{
  This is a (S3) generic function for which we provide methods,
  currently for \code{\link{nlrob}} only.% TODO
%%   \code{\link{lmrob}},
%%   \code{\link{glmrob}},
%% , and \code{\link{nls}}.
}
\seealso{
  \code{\link{nlrob}}, and \code{\link{nlrob.MM}}, notably for examples.
}
\keyword{models}
