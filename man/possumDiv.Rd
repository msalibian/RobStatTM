\name{possumDiv}
\alias{possumDiv}
\alias{possum.mat}
\docType{data}
\title{Possum Diversity Data}
\description{
  Possum diversity data: As issued from a study of the diversity
  of possum (arboreal marsupials) in the Montane ash forest (Australia),
  this dataset was collected in view of the management of hardwood
  forest to take conservation and recreation values, as well as wood
  production, into account.

  The study is fully described in the two references.
  The number of different species of arboreal marsupials (possum) was
  observed on 151 different 3ha sites with uniform vegetation. For each
  site the nine variable measures (see below) were recorded.
  The problem is to model the relationship between \code{diversity} and these
  other variables.
}
\usage{data(possumDiv)}
\format{
  Two different representations of the same data are available:

  \code{possumDiv} is a data frame of 151 observations
  of 9 variables, where the last two are factors, \code{eucalyptus} with
  3 levels and \code{aspect} with 4 levels.

  \code{possum.mat} is a numeric (integer) matrix of 151 rows
  (observations) and 14 columns (variables) where the last seven ones
  are 0-1 dummy variables, three (\code{E.*}) are coding for the kind of
  \code{eucalyptus} and the last four are 0-1 coding for the
  \code{aspect} factor.

  The variables have the following meaning:
  \describe{
    \item{Diversity}{main variable of interest is the number of
      different species of arboreal marsupial (possum) observed, with
      values in 0:5.}
    \item{Shrubs}{the number of shrubs.}
    \item{Stumps}{the number of cut stumps from past logging operations.}
    \item{Stags}{the number of stags (hollow-bearing trees).}
    \item{Bark}{bark index (integer) vector reflecting the quantity of
      decorticating bark.}
    \item{Habitat}{an integer score indicating the suitability of
      nesting and foraging habitat for Leadbeater's possum.}
    \item{BAcacia}{a numeric vector giving the basal area of acacia species.}
    \cr

    \item{eucalyptus}{a 3-level \code{\link{factor}}
      specifying the species of eucalypt with the greatest stand basal
      area.  This has the same information as the following three variables}
    \item{E.regnans}{0-1 indicator for Eucalyptus regnans}
    \item{E.delegatensis}{0-1 indicator for Eucalyptus deleg.}
    \item{E.nitens}{0-1 indicator for Eucalyptus nitens}
    \cr

    \item{aspect}{a 4-level \code{\link{factor}} specifying the aspect
      of the site.  It is the same information as the following four
      variables.}
    \item{NW-NE}{0-1 indicator}
    \item{NW-SE}{0-1 indicator}
    \item{SE-SW}{0-1 indicator}
    \item{SW-NW}{0-1 indicator}
  }
}
\source{
  Eva Cantoni (2004)
  Analysis of Robust Quasi-deviances for Generalized Linear Models.
  \emph{Journal of Statistical Software} \bold{10}, 04,
  \url{http://www.jstatsoft.org/v10/i04}
}
\references{
  Lindenmayer, D. B., Cunningham, R. B., Tanton, M. T., Nix, H. A. and
  Smith, A. P. (1991)
  The conservation of arboreal marsupials in the montane ash forests of
  the central highlands of victoria, south-east australia: III. The habitat
  requirements of leadbeater's possum \emph{gymnobelideus leadbeateri} and
  models of the diversity and abundance of arboreal marsupials.
  \emph{Biological Conservation} \bold{56}, 295--315.

  Lindenmayer, D. B., Cunningham, R. B., Tanton, M. T., Smith, A. P. and
  Nix, H. A. (1990)
  The conservation of arboreal marsupials in the montane ash forests of
  the victoria, south-east australia, I. Factors influencing the occupancy of
  trees with hollows, \emph{Biological Conservation} \bold{54}, 111--131.

  See also the references in \code{\link{glmrob}}.
}
\examples{
data(possumDiv)
head(possum.mat)

str(possumDiv)
## summarize all variables as multilevel factors:
summary(as.data.frame(lapply(possumDiv, function(v)
                             if(is.integer(v)) factor(v) else v)))

## Following Cantoni & Ronchetti (2001), JASA, p.1026 f.:% cf. ../tests/poisson-ex.R
pdFit <- glmrob(Diversity ~ . , data = possumDiv,
                family=poisson, tcc = 1.6, weights.on.x = "hat", acc = 1e-15)
summary(pdFit)
summary(pdF2 <- update(pdFit, ~ . -Shrubs))
summary(pdF3 <- update(pdF2,  ~ . -eucalyptus))
summary(pdF4 <- update(pdF3,  ~ . -Stumps))
summary(pdF5 <- update(pdF4,  ~ . -BAcacia))
summary(pdF6 <- update(pdF5,  ~ . -aspect))# too much ..
anova(pdFit, pdF3, pdF4, pdF5, pdF6, test = "QD") # indeed,
## indeed, the last simplification is too much
possumD.2 <- within(possumDiv, levels(aspect)[1:3] <- rep("other", 3))
## and use this binary 'aspect' instead of the 4-level one:
summary(pdF5.1 <- update(pdF5, data = possumD.2))

if(FALSE) # not ok, as formually not nested.
anova(pdF5, pdF5.1)

summarizeRobWeights(weights(pdF5.1, type="rob"), eps = 0.73)
##->  "outliers"  (1, 59, 110)
wrob <- setNames(weights(pdF5.1, type="rob"), rownames(possumDiv))
head(sort(wrob))
}
\keyword{datasets}
