% Origin: src/library/stats/man/glm.summaries.Rd (as of 2011-10-23)
\name{residuals.glmrob}
\alias{residuals.glmrob}
\title{Residuals of Robust Generalized Linear Model Fits}
\usage{
\method{residuals}{glmrob}(object,
          type = c("deviance", "pearson", "working",
                   "response", "partial"),
          \dots)
}
\arguments{
  \item{object}{an object of class \code{glmrob}, typically the result of
    a call to \code{\link{glmrob}}.}
  \item{type}{the type of residuals which should be returned.
    The alternatives are: \code{"deviance"} (default), \code{"pearson"},
    \code{"working"}, \code{"response"}, and \code{"partial"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  Compute residuals of a fitted \code{\link{glmrob}} model, i.e., robust
  generalized linear model fit.
}
\details{
  The references in \code{\link{glm}} define the types of residuals:
  Davison & Snell is a good reference for the usages of each.

  The partial residuals are a matrix of working residuals, with each
  column formed by omitting a term from the model.

  The \code{residuals} (S3) method (see \code{\link{methods}}) for
  \code{\link{glmrob}} models has been modeled to follow closely the
  method for classical (non-robust) \code{\link{glm}} fitted models.
  Possibly, see its documentation, i.e., \link{residuals.glm}, for
  further details.
}
\seealso{
  \code{\link{glmrob}} for computing \code{object}, \code{\link{anova.glmrob}};
  the corresponding \emph{generic} functions, \code{\link{summary.glmrob}},
  \code{\link{coef}},
  % \code{\link{deviance}}, \code{\link{effects}},
  \code{\link{fitted}},
  \code{\link{residuals}}.

}
\references{
  See those for the classical GLM's, \code{\link{glm}}.
}
\examples{
### -------- Gamma family -- data from example(glm) ---
clotting <- data.frame(
            u = c(5,10,15,20,30,40,60,80,100),
         lot1 = c(118,58,42,35,27,25,21,19,18),
         lot2 = c(69,35,26,21,18,16,13,12,12))
summary(cl <- glm   (lot1 ~ log(u), data=clotting, family=Gamma))
summary(ro <- glmrob(lot1 ~ log(u), data=clotting, family=Gamma))
clotM5.high <- within(clotting, { lot1[5] <- 60 })
cl5.high <- glm   (lot1 ~ log(u), data=clotM5.high, family=Gamma)
ro5.high <- glmrob(lot1 ~ log(u), data=clotM5.high, family=Gamma)
rr <- range(residuals(ro), residuals(cl), residuals(ro5.high))
plot(residuals(ro5.high) ~ residuals(cl5.high), xlim = rr, ylim = rr, asp = 1)
abline(0,1, col=2, lty=3)
points(residuals(ro) ~ residuals(cl), col = "gray", pch=3)

## Show all kinds of residuals:
r.types <- c("deviance", "pearson", "working", "response")
sapply(r.types, residuals, object = ro5.high)
}
\keyword{models}
\keyword{regression}
