\name{coleman}
\alias{coleman}
\docType{data}
\title{Coleman Data Set}
\description{
  Contains information on 20 Schools from the Mid-Atlantic and New
  England States, drawn from a population studied by Coleman et
  al. (1966). Mosteller and Tukey (1977) analyze this sample consisting
  of measurements on six different variables, one of which will be
  treated as a responce.
}
\usage{data(coleman)}
\format{
  A data frame with 20 observations on the following 6 variables.
  \describe{
    \item{\code{salaryP}}{staff salaries per pupil}
    \item{\code{fatherWc}}{percent of white-collar fathers}
    \item{\code{sstatus}}{socioeconomic status composite deviation: means for
      family size, family intactness, father's education, mother's
      education, and home items}
    \item{\code{teacherSc}}{mean teacher's verbal test score}
    \item{\code{motherLev}}{mean mother's educational level, one unit is equal
      to two school years}
    \item{\code{Y}}{verbal mean test score (y, all sixth graders)}
  }
}
\author{Valentin Todorov}
\source{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection} Wiley, p.79, table 2.
}
\examples{
data(coleman)
pairs(coleman)
summary( lm.coleman <-     lm(Y ~ . , data = coleman))
summary(lts.coleman <- ltsReg(Y ~ . , data = coleman))

coleman.x <- data.matrix(coleman[, 1:6])
(Cc <- covMcd(coleman.x))
}
\keyword{datasets}
