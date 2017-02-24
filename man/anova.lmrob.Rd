\name{anova.lmrob}
\alias{anova.lmrob}
\title{Analysis of Robust Deviances ('anova') for "lmrob" Objects}
\description{
  Compute an analysis of robust Wald-type or deviance-type test tables
  for one or more linear regression models fitted by \code{\link{lmrob}}.
}
\usage{
\method{anova}{lmrob}(object, \dots, test = c("Wald", "Deviance"),
      verbose = getOption("verbose"))
}
\arguments{
  \item{object, \dots}{objects of class \code{"lmrob"}, typically
    the result of a call to \code{\link{lmrob}}. \code{\dots} arguments
    may also be symbolic descriptions of the reduced models
    (cf. argument \code{formula} in \code{\link{lm}}).
  }
  \item{test}{a character string specifying the test statistic to be used.
    Can be one of \code{"Wald"} or \code{"Deviance"},
    with partial matching allowed, for specifying a \code{"Wald"}-type
    test or \code{"Deviance"}-type test.}
  \item{verbose}{logical; if true some informative messages are printed.}
}
\details{
  Specifying a single object gives a sequential analysis of a robust
  quasi-deviance table for that fit.  That is, the reductions in the
  robust residual deviance as each term of the formula is added in turn
  are given in as the rows of a table. (Currently not yet implemented.)

  If more than one object is specified, the table has a row for the
  residual quasi-degrees of freedom (however, this information is never
  used in the asymptotic tests).  For all but the first model, the
  change in degrees of freedom and robust deviance is also given.  (This
  only makes statistical sense if the models are nested.)  As opposed to
  the convention, the models are forced to be listed from largest to
  smallest due to computational reasons.

  In addition, the table will contain test statistics and P values
  comparing the reduction in robust deviances for the model on the row
  to that on top of it.  There are two different robust tests available:
  The "Wald"-type test (\code{test = "Wald"}) and the Deviance-type test
  (\code{test = "Deviance"}).  When using formula description of the
  nested models in the dot arguments and \code{test = "Deviance"}, you
  may be urged to supply a \code{\link{lmrob}} fit for these models by
  an error message.  This happens when the coefficients of the largest
  model reduced to the nested models result in invalid initial estimates
  for the nested models (indicated by robustness weights which are all 0).

  The comparison between two or more models by \code{\link{anova.lmrob}}
  will only be valid if they are fitted to the same dataset.
}

\value{
  Basically, an object of class \code{\link{anova}} inheriting from class
  \code{\link{data.frame}}.
}

%%\references{ ~put references to the literature/web site here ~ }
\author{Andreas Ruckstuhl}

\seealso{\code{\link{lmrob}}, \code{\link{anova}}.
}
\examples{
data(salinity)
summary(m0.sali  <- lmrob(Y ~ . , data = salinity))
anova(m0.sali, Y ~ X1 + X3)
## -> X2 is not needed
(m1.sali  <- lmrob(Y ~ X1 + X3, data = salinity))
anova(m0.sali, m1.sali) # the same as before
anova(m0.sali, m1.sali, test = "Deviance")
## whereas 'X3' is highly significant:
m2 <- update(m0.sali, ~ . -X3)
anova(m0.sali, m2)
anova(m0.sali, m2,  test = "Deviance")
## Global test [often not interesting]:
anova(m0.sali, update(m0.sali, . ~ 1), test = "Wald")
anova(m0.sali, update(m0.sali, . ~ 1), test = "Deviance")

if(require("MPV")) { ## Montgomery, Peck & Vining  datasets
  Jet <- table.b13
  Jet.rflm1 <- lmrob(y ~ ., data=Jet,
                     control = lmrob.control(max.it = 500))
  summary(Jet.rflm1)

  anova(Jet.rflm1, y ~ x1 + x5 + x6, test="Wald")

  try( anova(Jet.rflm1, y ~ x1 + x5 + x6, test="Deviance") )
  ## -> Error in anovaLm....  Please fit the nested models by lmrob

  ## {{ since  all robustness weights become 0 in the nested model ! }}

  ## Ok: Do as the error message told us:
  ##    test by comparing the two *fitted* models:

  Jet.rflm2 <- lmrob(y ~ x1 + x5 + x6, data=Jet,
                     control=lmrob.control(max.it=100))
  anova(Jet.rflm1, Jet.rflm2, test="Deviance")

} # end{"MPV" data}

}
\keyword{robust}
\keyword{models}
\keyword{regression}

