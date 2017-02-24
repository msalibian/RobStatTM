\name{pulpfiber}
\alias{pulpfiber}
\encoding{utf-8}
\docType{data}
\title{Pulp Fiber and Paper Data}
\description{
  Measurements of aspects pulp fibers and the paper produced from them.
  Four properties of each are measured in sixty-two samples.
}
\usage{data(pulpfiber)}
\format{
  A data frame with 62 observations on the following 8 variables.
  \describe{
    \item{\code{X1}}{numeric vector of arithmetic fiber length}
    \item{\code{X2}}{numeric vector of long fiber fraction}
    \item{\code{X3}}{numeric vector of fine fiber fraction}
    \item{\code{X4}}{numeric vector of zero span tensile}
%    \item{\ }{\ }% ------------- white space --
    \item{\code{Y1}}{numeric vector of breaking length}
    \item{\code{Y2}}{numeric vector of elastic modulus}
    \item{\code{Y3}}{numeric vector of stress at failure}
    \item{\code{Y4}}{numeric vector of burst strength}
  }
}
\details{
  Cited from the reference article:
  \emph{The dataset contains measurements of properties of pulp fibers and the
    paper made from them.  The aim is to investigate relations between pulp
    fiber properties and the resulting paper properties.  The dataset
    contains \eqn{n = 62} measurements of the following four pulp fiber
    characteristics: arithmetic fiber length, long fiber fraction, fine
    fiber fraction, and zero span tensile.  The four paper properties that
    have been measured are breaking length, elastic modulus, stress at
    failure, and burst strength.}

  The goal is to predict the \eqn{q = 4} paper properties from the
  \eqn{p = 4} fiber characteristics.
}
\source{
  Rousseeuw, P. J., Van Aelst, S., Van Driessen, K., and Agulló, J. (2004)
  Robust multivariate regression;
  \emph{Technometrics} \bold{46}, 293--305.

  \url{http://users.ugent.be/~svaelst/data/pulpfiber.txt}
}
\author{port to \R and this help page: Martin Maechler
}
\references{
  Lee, J. (1992)
  \emph{Relationships Between Properties of Pulp-Fibre and Paper},
  unpublished doctoral thesis, U. Toronto, Faculty of Forestry.
}
\examples{
data(pulpfiber)
str(pulpfiber)

pairs(pulpfiber, gap=.1)
## 2 blocks of 4 ..
c1 <- cov(pulpfiber)
cR <- covMcd(pulpfiber)
## how different are they: The robust estimate has more clear high correlations:
symnum(cov2cor(c1))
symnum(cov2cor(cR$cov))
}
\keyword{datasets}
