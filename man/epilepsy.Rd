\name{epilepsy}
\alias{epilepsy}
\docType{data}
\title{Epilepsy Attacks Data Set}
\description{Data from a clinical trial of 59 patients with epilepsy
  (Breslow, 1996) in order to illustrate diagnostic techniques in
  Poisson regression.
}
\usage{data(epilepsy)}
\format{
  A data frame with 59 observations on the following 11 variables.
  \describe{
    \item{\code{ID}}{Patient identification number}
    \item{\code{Y1}}{Number of epilepsy attacks patients have during the
      first follow-up period}
    \item{\code{Y2}}{Number of epilepsy attacks patients have during the
      second follow-up period}
    \item{\code{Y3}}{Number of epilepsy attacks patients have during the
      third follow-up period}
    \item{\code{Y4}}{Number of epilepsy attacks patients have during the
      forth follow-up period}
    \item{\code{Base}}{Number of epileptic attacks
      recorded during 8 week period prior to randomization}
    \item{\code{Age}}{Age of the patients}
    \item{\code{Trt}}{a factor with levels \code{placebo}
      \code{progabide} indicating whether the anti-epilepsy
      drug Progabide has been applied or not}
    \item{\code{Ysum}}{Total number of epilepsy attacks patients have
      during the four follow-up periods  }
    \item{\code{Age10}}{Age of the patients devided by 10}
    \item{\code{Base4}}{Variable \code{Base} devided by 4}
  }
}
\details{Thall and Vail reported data from  a clinical trial of 59 patients
  with epilepsy, 31 of whom were randomized to receive the anti-epilepsy
  drug Progabide and 28 of whom received a placebo. Baseline data
  consisted of the patient's age and the number of epileptic seizures
  recorded during 8 week period prior to randomization. The response
  consisted of counts of seizures occuring during the four consecutive
  follow-up periods of two weeks each.
}
\source{
  Thall, P.F. and Vail S.C. (1990)
  Some covariance models for longitudinal count data with overdispersion.
  \emph{Biometrics} \bold{46}, 657--671.
}
\references{
  Diggle, P.J., Liang, K.Y., and Zeger, S.L. (1994)
  \emph{Analysis of Longitudinal Data}; Clarendon Press.

  Breslow N. E. (1996)
  Generalized linear models: Checking assumptions and strengthening
  conclusions.
  \emph{Statistica Applicata} \bold{8}, 23--41.
}
\examples{
data(epilepsy)
str(epilepsy)
pairs(epilepsy[,c("Ysum","Base4","Trt","Age10")])

Efit1 <- glm(Ysum ~ Age10 + Base4*Trt, family=poisson, data=epilepsy)
summary(Efit1)

## Robust Fit :
Efit2 <- glmrob(Ysum ~ Age10 + Base4*Trt, family=poisson, data=epilepsy,
                method = "Mqle",
                tcc=1.2, maxit=100)
summary(Efit2)
}
\keyword{datasets}
