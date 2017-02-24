\name{nlrob.control}
\alias{nlrob.control}
\title{Control Nonlinear Robust Regression Algorithms}
\description{
  Allow the user to specify details for the different nonlinear robust
  regression algorithms in \code{\link{nlrob}}.
}
\usage{
nlrob.control(method,
              psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
              init = c("S", "lts"), 
              optimizer = "JDEoptim", optArgs  = list(),
              ...)
}
\arguments{
  \item{method}{\code{\link{character}} string specifying the method}
  \item{psi}{string specifying the psi-function which defines the estimator.}
  \item{init}{for some methods, currently, \code{"MM"} only, a string
    specifying the initial estimator.
  }
  \item{optimizer}{currently only \code{"JDEoptim"} from package \pkg{DEoptimR}.}
  \item{optArgs}{
    a \code{\link{list}} of optional arguments to the optimizer.
    Currently, that is \code{\link[DEoptimR]{JDEoptim}} from package
    \pkg{DEoptimR}.
  }
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}
%% \details{
%%
%% }
\value{
  a \code{\link{list}} with several named components.
  The contents depend quite a bit on the \code{method}.
}
%\author{Martin Maechler}
%% \note{
%% }
\seealso{
  \code{\link{nlrob}},   \code{\link{nlrob}},

}
\examples{
str(nlrob.control("MM"))
str(nlrob.control("tau"))
str(nlrob.control("CM"))
str(nlrob.control("mtl"))
}
\keyword{utilities}
