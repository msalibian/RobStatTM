\name{anova.glmrob}
\alias{anova.glmrob}
\title{Analysis of Robust Quasi-Deviance for "glmrob" Objects}
\description{
  Compute an analysis of robust quasi-deviance table for one or more
  generalized linear models fitted by \code{\link{glmrob}}.
}
\usage{
\method{anova}{glmrob}(object, ..., test = c("Wald", "QD", "QDapprox"))
}
\arguments{
  \item{object, \dots}{objects of class \code{glmrob}, typically
    the result of a call to \code{\link{glmrob}}.}
  \item{test}{a character string specifying the test statistic to be
    used.  (Partially) matching one of \code{"Wald"}, \code{"QD"} or
    \code{"QDapprox"}. See Details.}
}
\details{
  Specifying a single object gives a sequential analysis of robust
  quasi-deviance table for that fit.  That is, the reductions in the
  robust residual quasi-deviance as each term of the formula is added in
  turn are given in as the rows of a table. \emph{(Currently not yet
  implemented.)}

  If more than one object is specified, the table has a row for the
  residual quasi-degrees of freedom (However, this information is never used in
  the asymptotic tests).  For all but the first model, the
  change in degrees of freedom and robust quasi-deviance is also
  given.  (This only makes statistical sense if the models are nested.)
  It is conventional to list the models from smallest to largest, but
  this is up to the user.

  In addition, the table will contain test statistics and P values
  comparing the reduction in robust quasi-deviance for the model on the
  row to that on top of it.  For all robust fitting methods, the
  \dQuote{Wald}-type test between two models can be applied (\code{test
    = "Wald"}).

  When using Mallows or Huber type robust estimators
  (\code{method="Mqle"} in \code{\link{glmrob}}), then there are
  additional test methods.  One is the robust quasi-deviance test
  (\code{test = "QD"}), as described by Cantoni and Ronchetti (2001).
  The asymptotic distribution is approximated by a chi-square
  distibution.  Another test (\code{test = "QDapprox"}) is based on a
  quadratic approximation of the robust quasi-deviance test
  statistic.  Its asymptotic distribution is chi-square (see the reference).

  The comparison between two or more models by \code{anova.glmrob}
  will only be valid if they are fitted to the same dataset and by the same
  robust fitting method using the same tuning constant \eqn{c} (\code{tcc} in
  \code{\link{glmrob}}).
}
\value{
  Basically, an object of class \code{\link{anova}} inheriting from class
  \code{\link{data.frame}}.
}
\references{
  E. Cantoni and E. Ronchetti (2001)
  Robust Inference for Generalized Linear Models.
  \emph{JASA} \bold{96} (455), 1022--1030.

  E.Cantoni (2004)
  Analysis of Robust Quasi-deviances for Generalized Linear Models.
  \emph{Journal of Statistical Software} \bold{10},
  \url{http://www.jstatsoft.org/v10/i04}
}

\author{ Andreas Ruckstuhl }

\seealso{ \code{\link{glmrob}}, \code{\link{anova}}.
%%
%%   \code{\link{drop1}} for
%%   so-called \sQuote{type II} anova where each term is dropped one at a
%%   time respecting their hierarchy.
}
\examples{
## Binomial response -----------
data(carrots)
Cfit2 <- glmrob(cbind(success, total-success) ~ logdose + block,
                family=binomial, data=carrots, method="Mqle",
                control=glmrobMqle.control(tcc=1.2))
summary(Cfit2)

Cfit4 <- glmrob(cbind(success, total-success) ~ logdose * block,
                family=binomial, data=carrots, method="Mqle",
                control=glmrobMqle.control(tcc=1.2))

anova(Cfit2, Cfit4, test="Wald")

anova(Cfit2, Cfit4, test="QD")

anova(Cfit2, Cfit4, test="QDapprox")

## Poisson response ------------
data(epilepsy)

Efit2 <- glmrob(Ysum ~ Age10 + Base4*Trt, family=poisson, data=epilepsy,
               method="Mqle", control=glmrobMqle.control(tcc=1.2,maxit=100))
summary(Efit2)

Efit3 <- glmrob(Ysum ~ Age10 + Base4 + Trt, family=poisson, data=epilepsy,
               method="Mqle", control=glmrobMqle.control(tcc=1.2,maxit=100))

anova(Efit3, Efit2, test = "Wald")

anova(Efit3, Efit2, test = "QD")

## trivial intercept-only-model:
E0 <- update(Efit3, . ~ 1)
anova(E0, Efit3, Efit2, test = "QDapprox")
%% failed in robustbase <= 2013-11-27
}
\keyword{robust}
\keyword{models}
\keyword{regression}
