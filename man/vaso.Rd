\name{vaso}
\alias{vaso}
\docType{data}
\title{Vaso Constriction Skin Data Set}
\description{
  Finney's data on vaso constriction in the skin of the digits.
}
\usage{data(vaso)}
\format{
  A data frame with 39 observations on the following 3 variables.
  \describe{
    \item{\code{Volume}}{Inhaled volume of air}
    \item{\code{Rate}}{Rate of inhalation}
    \item{\code{Y}}{vector of 0 or 1 values.}
  }
}
\details{The data taken from Finney (1947) were obtained in a carefully
  controlled study in human physiology where a reflex
  \dQuote{vaso constriction} may occur in the skin of the digits after taking a
  single deep breath.  The response y is the occurence (y = 1) or
  non-occurence (y = 0) of vaso constriction in the skin of the digits
  of a subject after he or she inhaled a certain volume of air at a certain
  rate.  The responses of three subjects are available.  The first
  contributed 9 responses, the second contributed 8 responses, and the
  third contributed 22 responses.

  Although the data represent repeated measurements, an analysis that
  assumes independent observations may be applied, as claimed by Pregibon
  (1981).
}
\source{
  Finney, D.J. (1947)
  The estimation from individual records of the relationship between
  dose and quantal response.
  \emph{Biometrika} \bold{34}, 320--334
}
\references{
  Atkinson, A.C. and Riani, M. (2000)
  \emph{Robust Diagnostic Regression Analysis},
  First Edition. New York: Springer, Table A.23.

  Fahrmeir, L. and Tutz, G. (2001)
  \emph{Multivariate Statistical Modelling Based on Generalized Linear Models},
  Springer, Table 4.2.

  Kuensch, H.R., Stefanski, A. and Carrol, R.J. (1989)
  Conditionally unbiased bounded influence estimation in general
  regression models, with applications to generalized linear models,
  \emph{JASA} \bold{84}, 460--466.

  Pregibon, D. (1981)
  Logistic regression diagnostics,
  \emph{Annals of Statistics} \bold{9}, 705--724.
}
\examples{
data(vaso)
str(vaso)
pairs(vaso)

glmV <- glm(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso)
summary(glmV)
## -->  example(glmrob)  showing classical & robust GLM
}
\keyword{datasets}
