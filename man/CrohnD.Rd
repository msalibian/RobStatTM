\name{CrohnD}
\encoding{utf8}
\alias{CrohnD}
\docType{data}
\title{Crohn's Disease Adverse Events Data}
\description{
  Data set issued from a study of the adverse events of a drug on
  117 patients affected by Crohn's disease (a chronic inflammatory
  disease of the intestines).
}
\usage{data(CrohnD)}
\format{
  A data frame with 117 observations on the following 9 variables.
  \describe{
%% FIXME: leave these away -- or explain: -- they code patient sub-groups
    \item{\code{ID}}{the numeric patient IDs}
    \item{\code{nrAdvE}}{the number of adverse events}
    \item{\code{BMI}}{Body MASS Index, i.e., \eqn{weight[kg] / (height[m])^2}.}
    \item{\code{height}}{in cm}
    \item{\code{country}}{a factor with levels \code{0} and \code{1}}
    \item{\code{sex}}{the person's gender, a binary factor with levels
      \code{M} \code{F}}
    \item{\code{age}}{in years, a numeric vector}
    \item{\code{weight}}{in kilograms, a numeric vector}
    \item{\code{treat}}{how CD was treated: a factor with levels
      \code{0}, \code{1} and \code{2}, meaning placebo, drug 1 and drug 2.}
  }
}
% \details{
%   ~~ If necessary, more details than the __description__ above ~~
% }
\source{
  form the authors of the reference, with permission by the original
  data collecting agency.
}
\references{
  Serigne N. Lô and Elvezio Ronchetti (2006).
  Robust Second Order Accurate Inference for Generalized Linear Models.
  Technical report, University of Geneva, Switzerland.
}
\examples{
data(CrohnD)
str(CrohnD)
with(CrohnD, ftable(table(sex,country, treat)))
}
\keyword{datasets}
